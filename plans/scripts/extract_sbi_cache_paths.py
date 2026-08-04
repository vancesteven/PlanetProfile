#!/usr/bin/env python3
"""Print unique serve-time cache paths from the Streamlit SBI registry.

The page module is deliberately not imported: importing it executes Streamlit
page setup.  ``ast.literal_eval`` accepts only the registry's literal data and
cannot execute source code.
"""
from __future__ import annotations

import argparse
import ast
from pathlib import Path


REGISTRY_NAME = '_SBI_ARTIFACT_SLOTS'


def extract_cache_paths(source: Path) -> list[str]:
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

    paths = set()
    for slot_id, slot in registry.items():
        if not isinstance(slot, dict):
            raise TypeError(f'slot {slot_id!r} must be a dict literal')
        cache_path = slot.get('cache_path')
        if cache_path is None:  # intentional placeholder slot
            continue
        if not isinstance(cache_path, str) or not cache_path:
            raise TypeError(f'slot {slot_id!r} has invalid cache_path')
        path = Path(cache_path)
        if path.is_absolute() or '..' in path.parts:
            raise ValueError(f'slot {slot_id!r} cache_path is not repo-relative')
        paths.add(path.as_posix())
    return sorted(paths)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument('source', type=Path,
                        help='Inference.py containing _SBI_ARTIFACT_SLOTS')
    args = parser.parse_args()
    print(*extract_cache_paths(args.source), sep='\n')


if __name__ == '__main__':
    main()
