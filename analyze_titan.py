#!/usr/bin/env python3
"""Analyze Titan model outputs for results section."""

import pandas as pd

def analyze_profile(filename):
    """Extract phase transitions and key parameters from profile."""
    print(f"\n{'='*60}")
    print(f"Analyzing: {filename}")
    print(f"{'='*60}")

    with open(filename, 'r') as f:
        lines = f.readlines()

    # Extract header info
    for line in lines[:82]:
        if 'Tb_K' in line:
            print(line.strip())
        if 'zb_km' in line:
            print(line.strip())
        if 'D_km' in line:
            print(line.strip())
        if 'Pb_MPa' in line:
            print(line.strip())
        if 'RaConvect' in line and 'RaConvectIII' not in line and 'RaConvectV' not in line:
            print(line.strip())
        if 'etaConv_Pas' in line:
            print(line.strip())
        if 'eLid_m' in line and 'eLidIII' not in line and 'eLidV' not in line:
            print(line.strip())
        if 'Dconv_m' in line and 'DconvIII' not in line and 'DconvV' not in line:
            print(line.strip())

    # Parse data
    data_lines = []
    for line in lines[82:]:
        if line.strip() and not line.startswith('P (MPa)'):
            parts = line.split()
            if len(parts) >= 4:
                try:
                    data_lines.append({
                        'P_MPa': float(parts[0]),
                        'T_K': float(parts[1]),
                        'r_m': float(parts[2]),
                        'phase': int(parts[3])
                    })
                except:
                    pass

    df = pd.DataFrame(data_lines)

    # Find phase transitions
    print("\nPhase transitions:")
    prev_phase = None
    for idx, row in df.iterrows():
        if row['phase'] != prev_phase:
            r_km = row['r_m'] / 1e3
            depth_km = 2574.73 - r_km  # Titan radius
            print(f"  Phase {int(row['phase'])}: P={row['P_MPa']:.1f} MPa, T={row['T_K']:.1f} K, "
                  f"r={r_km:.1f} km (depth={depth_km:.1f} km)")
            prev_phase = row['phase']

    print(f"\nTotal layers: {len(df)}")
    print(f"Phases present: {sorted(df['phase'].unique())}")

    return df

# Analyze both cases
df_ocean = analyze_profile('Titan/TitanProfile_NaCl_100.0ppt_Tb244.5K_ConstantInnerRho.txt')
df_no_ocean = analyze_profile('Titan/TitanProfile_NaCl_100.0ppt_Tb243.165K_450000.0_ConstantInnerRho.txt')

print("\n" + "="*60)
print("COMPARISON SUMMARY")
print("="*60)
print("Ocean case: Has liquid ocean layer")
print("No ocean case: Ice Ih transitions directly to HP ices")
