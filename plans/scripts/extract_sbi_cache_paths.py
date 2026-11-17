#!/usr/bin/env python3
"""Print unique serve-time cache paths from the Streamlit SBI registry.

The page module is deliberately not imported: importing it executes Streamlit
page setup.  ``ast.literal_eval`` accepts only the registry's literal data and
cannot execute source code.
"""
from __future__ import annotations

import argparse
import ast
import importlib.util
import json
import pickle
import sys
from pathlib import Path


REGISTRY_NAME = '_SBI_ARTIFACT_SLOTS'


def _load_registry(source: Path) -> dict:
    tree = ast.parse(source.read_text(encoding='utf-8'), filename=str(source))
    registry_node = None
    for node in tree.body:
        if not isinstance(node, ast.Assign):
            continue
        if any(isinstance(target, ast.Name) and target.id == REGISTRY_NAME
               for target in node.targets):
            registry_node = node.value
            break
    if registry_node is None:
        raise ValueError(f'{REGISTRY_NAME} assignment not found in {source}')

    registry = ast.literal_eval(registry_node)
    if not isinstance(registry, dict):
        raise TypeError(f'{REGISTRY_NAME} must be a dict literal')
    return registry


def _repo_relative(value, what: str) -> str:
    if not isinstance(value, str) or not value:
        raise TypeError(f'{what} is invalid')
    path = Path(value)
    if path.is_absolute() or '..' in path.parts:
        raise ValueError(f'{what} is not repo-relative')
    return path.as_posix()


def extract_cache_paths(source: Path) -> list[str]:
    paths = set()
    for slot_id, slot in _load_registry(source).items():
        if not isinstance(slot, dict):
            raise TypeError(f'slot {slot_id!r} must be a dict literal')
        cache_path = slot.get('cache_path')
        if cache_path is None:  # intentional placeholder slot
            continue
        paths.add(_repo_relative(cache_path, f'slot {slot_id!r} cache_path'))
    return sorted(paths)


def _config_induction_labels(raw: dict) -> bool:
    """True when the runner's Ae precompute would run for this config (an
    induction observable or an induction_bounds label; mirrors
    MCMCRunner._precompute_ae_grid's gate)."""
    obs = raw.get('observables') or {}
    if raw.get('induction_bounds'):
        return True
    return any(k.startswith(('Ae_', 'BiAmp_', 'BiPhase_', 'Bind_'))
               for k in obs)


def build_manifest(source: Path, repo_root: Path) -> dict:
    """Per-slot required/optional companion files, derived from the registry
    and each slot's config (nothing hard-coded except the Ae-sidecar and
    ``_offsets.json`` naming conventions the runner itself uses).

    Returns ``{slot_id: [{'path', 'kind', 'required'}, ...]}`` (placeholder
    slots with ``cache_path`` None are omitted).  Kinds:

    - ``cache``: the structure-grid pkl (required).
    - ``ae_sidecar``: ``<cache>.ae_sidecar.pkl``; required when the slot's
      config triggers the Ae precompute (induction observable/bounds) or
      declares ``metadata.induction_provenance`` (Enceladus display reads it).
    - ``cmr2_offsets``: ``<cache stem>_offsets.json`` (optional, shipped if
      present; the runner treats absence as 'no offset').
    - ``lib_factor_sidecar``: ``<cache> + metadata.libration_model_correction.
      elastic_factor_sidecar.suffix`` (required when declared; the runner
      raises on every servable draw if it is missing).
    - ``is_fiducial_manifest`` / ``is_fiducial_npz``: slot ``is_correction``
      bundle (required; ships via the ``PlanetProfile/`` rsync, not Test/).
    """
    out = {}
    for slot_id, slot in _load_registry(source).items():
        cache = slot.get('cache_path')
        if cache is None:
            continue
        cache = _repo_relative(cache, f'slot {slot_id!r} cache_path')
        files = [{'path': cache, 'kind': 'cache', 'required': True}]
        cfg_rel = slot.get('config_path')
        raw = {}
        if cfg_rel:
            raw = json.loads((repo_root / _repo_relative(
                cfg_rel, f'slot {slot_id!r} config_path')
            ).read_text(encoding='utf-8'))
        meta = raw.get('metadata') or {}
        ae_req = _config_induction_labels(raw) or 'induction_provenance' in meta
        files.append({'path': cache + '.ae_sidecar.pkl',
                      'kind': 'ae_sidecar', 'required': bool(ae_req)})
        files.append({'path': cache[:-len('.pkl')] + '_offsets.json'
                      if cache.endswith('.pkl') else cache + '_offsets.json',
                      'kind': 'cmr2_offsets', 'required': False})
        lib = (meta.get('libration_model_correction') or {}).get(
            'elastic_factor_sidecar')
        if lib:
            files.append({'path': cache + lib.get(
                'suffix', '.lib_factor_sidecar.pkl'),
                'kind': 'lib_factor_sidecar', 'required': True})
        isc = slot.get('is_correction')
        if isc:
            mpath = _repo_relative(isc['manifest'],
                                   f'slot {slot_id!r} is_correction.manifest')
            files.append({'path': mpath, 'kind': 'is_fiducial_manifest',
                          'required': True})
            mfile = repo_root / mpath
            if mfile.is_file():
                npz = json.loads(mfile.read_text(encoding='utf-8')).get(
                    'npz_file')
                if npz:
                    files.append({
                        'path': (Path(mpath).parent / npz).as_posix(),
                        'kind': 'is_fiducial_npz', 'required': True})
        out[slot_id] = files
    return out


def _load_key_module(repo_root: Path):
    """Load the stdlib-only Ae key definition by file path (no package
    import, so the deploy host needs neither numpy nor PlanetProfile)."""
    path = repo_root / 'PlanetProfile/Inference/ae_sidecar_key.py'
    spec = importlib.util.spec_from_file_location('_pp_ae_sidecar_key', path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def validate_ae_sidecars(source: Path, repo_root: Path) -> list[str]:
    """For each slot that REQUIRES an Ae sidecar, check the shipped legacy
    ``<cache>.ae_sidecar.pkl`` key equals the key that slot's runner would
    compute (labels from the config's observables/induction_bounds, node
    count and size of the cache file). Returns human-readable errors; empty
    means all shipped sidecars will be reused (not recomputed) at serve time.
    """
    keymod = _load_key_module(repo_root)
    errors = []
    seen = set()
    for slot_id, slot in _load_registry(source).items():
        cache = slot.get('cache_path')
        cfg_rel = slot.get('config_path')
        if cache is None or not cfg_rel:
            continue
        cache = _repo_relative(cache, f'slot {slot_id!r} cache_path')
        raw = json.loads((repo_root / _repo_relative(
            cfg_rel, f'slot {slot_id!r} config_path')).read_text('utf-8'))
        if not _config_induction_labels(raw):
            continue                      # runner would not build a sidecar
        sidecar = repo_root / (cache + '.ae_sidecar.pkl')
        cache_file = repo_root / cache
        if not sidecar.is_file() or not cache_file.is_file():
            continue                      # presence is the manifest's check
        labels = keymod.requested_ae_labels(
            raw.get('observables') or {}, raw.get('induction_bounds') or {})
        ident = (cache, tuple(sorted(labels)))
        if ident in seen:
            continue
        seen.add(ident)
        with open(cache_file, 'rb') as fh:
            data = pickle.load(fh)
        if not (isinstance(data, dict) and 'Tb_K_grid' in data
                and 'structures' in data):
            continue                      # runner skips the Ae precompute
        ours = keymod.sidecar_key(labels, len(data['structures']),
                                  cache_file.stat().st_size)
        try:
            with open(sidecar, 'rb') as fh:
                stored = pickle.load(fh).get('key')
        except Exception as exc:
            errors.append(f'slot {slot_id!r}: Ae sidecar {sidecar} is '
                          f'unreadable ({exc})')
            continue
        if stored != ours:
            errors.append(
                f'slot {slot_id!r}: shipped Ae sidecar {cache}.ae_sidecar.pkl '
                f'key {stored} != the key the slot runner computes {ours}; '
                'the public app would recompute the Ae grid on every runner '
                'build. Rebuild the sidecar (construct an MCMCRunner for the '
                "slot's config on the current cache) and commit it.")
    return errors


def flatten_manifest(manifest: dict) -> list[dict]:
    """Unique files across slots (required if any slot requires it), sorted."""
    merged = {}
    for files in manifest.values():
        for f in files:
            cur = merged.get(f['path'])
            if cur is None:
                merged[f['path']] = dict(f)
            else:
                cur['required'] = cur['required'] or f['required']
    return [merged[k] for k in sorted(merged)]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('source', type=Path,
                        help='Inference.py containing _SBI_ARTIFACT_SLOTS')
    parser.add_argument('--manifest', action='store_true',
                        help='emit the full companion-file manifest as '
                             'TSV lines: path<TAB>kind<TAB>required(0/1)')
    parser.add_argument('--json', action='store_true',
                        help='emit the per-slot manifest as JSON')
    parser.add_argument('--validate-ae', action='store_true',
                        help='check shipped Ae sidecar keys against the keys '
                             'the slot runners compute; exit 1 on mismatch')
    parser.add_argument('--repo-root', type=Path, default=None,
                        help='repo root for config lookup (default: '
                             'two levels above the source file)')
    args = parser.parse_args()
    if args.validate_ae:
        root = args.repo_root or args.source.resolve().parents[2]
        errs = validate_ae_sidecars(args.source, root)
        for e in errs:
            print('ERROR: ' + e, file=sys.stderr)
        sys.exit(1 if errs else 0)
    if args.manifest or args.json:
        root = args.repo_root or args.source.resolve().parents[2]
        man = build_manifest(args.source, root)
        if args.json:
            print(json.dumps(man, indent=1))
        else:
            for f in flatten_manifest(man):
                print(f"{f['path']}\t{f['kind']}\t{int(f['required'])}")
        return
    print(*extract_cache_paths(args.source), sep='\n')


if __name__ == '__main__':
    main()
