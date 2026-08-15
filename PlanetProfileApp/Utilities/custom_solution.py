"""Builders for CustomSolution ocean composition strings.

Kept as pure functions so the GUI table logic is unit-testable: the output
format `"CustomSolution<Name> = species: conc, ..."` is consumed verbatim
by Planet.Ocean.comp -> SpeciesParser (PlanetProfile HydroEOS) and must
not drift.
"""
import logging

log = logging.getLogger('PlanetProfile')


def build_custom_solution_comp(solution_name, rows):
    """Build the CustomSolution composition string from table rows.

    Args:
        solution_name: user-chosen name (whitespace-stripped; '' -> None).
        rows: iterable of (species, concentration). Rows with an empty
            species or a null/NaN/zero concentration are treated as the
            ion being ABSENT and dropped. Duplicate species keep the LAST
            occurrence.

    Returns:
        (comp_string_or_None, kept_rows, warnings): comp string is None
        when no valid rows or no name; kept_rows is the deduped
        [(species, conc)] actually used; warnings is a list of
        user-facing strings (e.g. duplicates dropped).
    """
    warnings = []
    kept = {}
    for species, conc in rows:
        name = (species or '').strip() if isinstance(species, str) else ''
        if not name:
            continue
        try:
            c = float(conc)
        except (TypeError, ValueError):
            continue
        if not (c > 0):  # None/NaN/0/negative -> ion absent
            continue
        if name in kept:
            warnings.append(f"Duplicate species {name}: keeping the last entry.")
        kept[name] = c

    solution_name = (solution_name or '').strip()
    if not kept or not solution_name:
        return None, list(kept.items()), warnings

    salt_string = ", ".join(f"{sp}: {c}" for sp, c in kept.items())
    return f"CustomSolution{solution_name} = {salt_string}", list(kept.items()), warnings
