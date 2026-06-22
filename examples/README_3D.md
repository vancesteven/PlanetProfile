# 3D Lateral Structure Examples

This directory contains example scripts demonstrating 3D lateral structure capabilities in PlanetProfile.

## Quick Start

### Simplest Example: Equilibrium Ice Thickness
```bash
python examples/europa_3d_simple.py
```

This runs Europa with the **recommended** equilibrium ice thickness mode:
- Physics-based ice thickness from steady-state heat balance
- Self-consistent thermal-tidal coupling
- ~5 minutes on 4 cores (nSide=4, 192 pixels)

### Full Comparison: 1D vs 3D
```bash
python compare_europa_3d.py
```

Runs both traditional 1D and 3D calculations with identical parameters. Generates:
- Comparison plots
- Quantitative metrics (JSON)
- Performance benchmarks

### Test Equilibrium Solver
```bash
python test_equilibrium_ice.py
```

Focused test of the equilibrium ice thickness solver with convergence diagnostics.

## Scripts Overview

| Script | Purpose | Run Time | Output |
|--------|---------|----------|--------|
| `europa_3d_simple.py` | Minimal equilibrium example | ~5 min | Lateral plots, NPZ data |
| `compare_europa_3d.py` | 1D vs 3D comparison | ~10 min | Comparison plots + JSON |
| `test_equilibrium_ice.py` | Equilibrium solver test | ~5 min | Convergence diagnostics |

## Ice Thickness Modes

PlanetProfile supports three modes for specifying ice thickness in 3D:

### 1. Equilibrium (RECOMMENDED for science)
```python
Planet.Lateral.DO_EQUILIBRIUM_ICE = True
Planet.Lateral.DO_TIDAL_3D = True
# Do NOT set dIce_Cpq_km — solver computes thickness
```

**Use when**: You need physically self-consistent ice thickness for publications.

**Physics**: Solves heat balance equation iteratively:
```
k * (Tb - Tsurf) / d_ice = q_basal + H_tidal(pixel) * d_ice
```

See [docs/3D_ICE_THICKNESS_MODES.md](../docs/3D_ICE_THICKNESS_MODES.md) for full details.

### 2. Prescribed (for exploratory studies)
```python
Planet.Lateral.DO_EQUILIBRIUM_ICE = False
Planet.Lateral.dIce_Cpq_km = np.array([[29, 0, 0], [0, 0, 0], [-3.5, 0, -1.5]])
Planet.Lateral.dIce_Spq_km = np.array([[0, 0, 0], [0, 0, 0], [0, 0, -0.7]])
```

**Use when**: Testing sensitivity to ice thickness patterns or matching observations.

### 3. Uniform (testing only)
```python
Planet.Lateral.DO_3D = True
# Neither DO_EQUILIBRIUM_ICE nor dIce_Cpq_km set
```

**Use when**: Debugging or baseline comparison.

## Command-Line Options

### Switch Between Modes
```bash
# Equilibrium mode (default)
python compare_europa_3d.py --mode equilibrium

# Prescribed mode
python compare_europa_3d.py --mode prescribed
```

### Change Grid Resolution
```bash
# Low resolution (fast, 192 pixels)
python compare_europa_3d.py --grid-size 4

# Medium resolution (moderate, 768 pixels)
python compare_europa_3d.py --grid-size 8

# High resolution (slow, 3072 pixels)
python compare_europa_3d.py --grid-size 16
```

## Expected Results

### Equilibrium Ice Thickness (Europa)
Based on Tobie et al. (2003, 2005):
- **Mean thickness**: ~25-30 km (depending on Tb_K)
- **Peak-to-peak variation**: ~5 km
- **Thickest**: Sub/anti-Jovian points (0°, 180° longitude)
- **Thinnest**: Mid-latitudes, poles
- **Surface heat flux**: ~35-40 mW/m², nearly uniform

### Performance
| Grid Resolution | Pixels | Equilibrium Time | Prescribed Time |
|----------------|--------|------------------|-----------------|
| nSide=4 | 192 | ~5 min | ~1 min |
| nSide=8 | 768 | ~20 min | ~4 min |
| nSide=16 | 3072 | ~80 min | ~16 min |

(4 cores, parallelized. Equilibrium mode takes 5-10 iterations.)

## Outputs

All scripts save results to `<BodyName>_<Config>/figures/`:
- `*_lateral3D.npz`: 3D field data (ice thickness, tidal heating, etc.)
- `*_lateral_ice_thickness.png`: Ice thickness map
- `*_lateral_tidal_heating.png`: Tidal heating map
- `*_lateral_surface_heat_flux.png`: Surface heat flux map

JSON files contain quantitative comparison metrics.

## Common Issues

### "DO_EQUILIBRIUM_ICE requires DO_TIDAL_3D"
**Solution**: Equilibrium mode needs tidal heating. Set:
```python
Planet.Lateral.DO_TIDAL_3D = True
```

### "Equilibrium mode requires orbital eccentricity"
**Solution**: Set orbital parameters:
```python
Planet.Bulk.eccentricity = 0.0094  # Europa value
Planet.Bulk.meanMotion_radps = 2 * np.pi / (3.551181 * 86400)  # Europa
```

### Equilibrium doesn't converge
**Solutions**:
1. Increase tolerance: `Planet.Lateral.equilibriumTol_m = 200.0`
2. More iterations: `Planet.Lateral.equilibriumMaxIter = 20`
3. Check thermal conductivity: `Planet.Lateral.kThermIce_WmK = 2.3` (reasonable?)
4. Check basal heat flux is realistic (~5-10 mW/m² for Europa)

### Script is too slow
**Solutions**:
1. Lower grid resolution: `nSide = 4` (192 pixels)
2. Use prescribed mode instead (no iteration)
3. Enable parallelization: `Params.DO_PARALLEL = True`
4. Disable plots during testing: `Params.PLOT_LATERAL = False`

## Modifying for Other Bodies

These scripts are easily adapted for Enceladus, Titan, Ganymede, etc.:

```python
# Example: Enceladus equilibrium ice thickness
Planet = PlanetStruct('Enceladus')
Planet.Bulk.R_m = 252.1e3
Planet.Bulk.M_kg = 1.08e20
Planet.Bulk.Tsurf_K = 75
Planet.Bulk.Tb_K = 273.0
Planet.Bulk.eccentricity = 0.0047
Planet.Bulk.meanMotion_radps = 2 * np.pi / (1.370218 * 86400)

# Same 3D config as Europa
Planet.Lateral.DO_3D = True
Planet.Lateral.DO_TIDAL_3D = True
Planet.Lateral.DO_EQUILIBRIUM_ICE = True
# ... etc
```

## References

1. **Tobie et al. (2003)**: Equilibrium ice thickness physics, *JGR* doi:10.1029/2003JE002099
2. **Tobie et al. (2005)**: Geographic tidal heating, *Icarus* 196, 642-652
3. **Ojakangas & Stevenson (1989)**: Tidal strain pattern, *Icarus* 81, 220-241
4. **Roberts & Nimmo (2008)**: Ice-ocean coupling, *Icarus* 194, 675-689

## Documentation

- **Full mode guide**: [docs/3D_ICE_THICKNESS_MODES.md](../docs/3D_ICE_THICKNESS_MODES.md)
- **Architecture**: [.claude/ARCHITECTURE.md](../.claude/ARCHITECTURE.md)
- **CLAUDE.md**: [CLAUDE.md](../CLAUDE.md)

---

**Last updated**: 2026-06-18  
**Author**: Emma Vellard  
**Questions**: emma.vellard@outlook.fr
