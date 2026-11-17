"""Shared runtime guard for the Enceladus ocean-branch SBI scripts.

The deployed artifact (sbi_artifacts/enceladus_ocean6D_rhorock_posterior_1m.pt)
and its training dataset were produced under the FROZEN libration convention.
The live config (enceladus_cassini_isostasy_ocean6D_rhorock.json) now carries
the REPAIRED libration treatment.  Gating / training the existing artifact or
dataset against the live config would pair the frozen-convention flow with the
repaired simulator.  Torch-free on purpose (libomp ordering in the gate script).
"""
import hashlib
import json
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
CFG_DIR = REPO / 'PlanetProfile/Inference/configs'
SNAPSHOT = CFG_DIR / 'enceladus_cassini_isostasy_ocean6D_rhorock_sbi_frozenconv.json'
ARTIFACT = REPO / 'PlanetProfile/Inference/sbi_artifacts/enceladus_ocean6D_rhorock_posterior_1m.pt'
COMMITTED_G0 = REPO / 'validation_reports/enceladus_isostasy/sbi/G0_report.json'


def _lib(path):
    return (json.loads(Path(path).read_text()).get('metadata') or {}).get(
        'libration_model_correction')


def _sha256(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for c in iter(lambda: f.read(1 << 20), b''):
            h.update(c)
    return h.hexdigest()


def is_existing_artifact(path):
    p = Path(path)
    return (ARTIFACT.exists() and p.exists()
            and (p.resolve() == ARTIFACT.resolve()
                 or _sha256(p) == _sha256(ARTIFACT)))


def is_existing_dataset_dir(d):
    """True if <d>/G0_report.json is the committed G0 record of the deployed
    artifact's training set (same cache sha and n_kept)."""
    g = Path(d) / 'G0_report.json'
    if not (g.exists() and COMMITTED_G0.exists()):
        return False
    a, b = json.loads(g.read_text()), json.loads(COMMITTED_G0.read_text())
    return (a.get('single_cache_sha') == b.get('single_cache_sha')
            and a.get('n_kept') == b.get('n_kept'))


def refuse_if_mispaired(cfg_path, what):
    """Exit unless ``cfg_path`` carries the frozen snapshot's
    libration_model_correction (``what`` names the existing artifact/dataset)."""
    if _lib(cfg_path) != _lib(SNAPSHOT):
        raise SystemExit(
            f'REFUSED: {what} is the existing FROZEN-libration-convention '
            f'artifact/dataset, but {cfg_path} has a different '
            'metadata.libration_model_correction (repaired convention); that '
            'would pair the frozen flow with the repaired simulator. To '
            f're-gate the current artifact use --config {SNAPSHOT}. A retrain '
            'under the repaired config is a NEW artifact and needs a new '
            'preregistration (see validation_reports/enceladus_isostasy/sbi/'
            'RESULT.md).')
