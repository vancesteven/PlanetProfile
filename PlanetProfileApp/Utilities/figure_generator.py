"""
Figure generation utilities for PlanetProfileApp.
Generates interactive Plotly figures from saved model data files.
"""

import os
import pandas as pd
import streamlit as st
from io import StringIO
import re

# Check if plotly is available
try:
    import plotly.graph_objects as go
    import plotly.express as px
    from plotly.subplots import make_subplots
    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False
    st.warning("⚠️ Plotly not installed. Install with 'pip install plotly' for interactive figures.")


def load_profile_data(file_path):
    """
    Load profile data from a PlanetProfile output text file.

    Args:
        file_path: Path to the profile text file

    Returns:
        pd.DataFrame with profile data
    """
    with open(file_path, 'r', encoding='utf-8') as f:
        lines = f.readlines()

    # Profile files have header info in lines 1-81, data starts at line 83
    data_str = "".join(lines[83:])
    df = pd.read_csv(StringIO(data_str), sep=r'\s+', header=None, na_values=['nan'])

    # Standard column names for profile files
    column_names = [
        "P (MPa)", "T (K)", "r (m)", "phase ID", "rho (kg/m3)",
        "Cp (J/kg/K)", "alpha (1/K)", "g (m/s2)", "phi (void/solid frac)",
        "sigma (S/m)", "k (W/m/K)", "VP (km/s)", "VS (km/s)", "QS",
        "KS (GPa)", "GS (GPa)", "Ppore (MPa)", "rhoMatrix (kg/m3)",
        "rhoPore (kg/m3)", "MLayer (kg)", "VLayer (m3)", "Htidal (W/m3)", "eta (Pa s)"
    ]

    if df.shape[1] == len(column_names):
        df.columns = column_names
    else:
        df.columns = [f"Col_{i+1}" for i in range(df.shape[1])]

    return df


def load_ocean_data(file_path):
    """
    Load ocean data from a liquidOceanProps.txt file.

    Args:
        file_path: Path to the ocean properties text file

    Returns:
        pd.DataFrame with ocean properties data
    """
    with open(file_path, "r", encoding="utf-8") as f:
        lines = f.readlines()

    header_line = lines[4].strip()
    column_names = re.split(r'\s{2,}|\t+', header_line)
    data_lines = "".join(lines[5:])

    df = pd.read_csv(StringIO(data_lines), sep=r'\s+', header=None,
                     engine='python', on_bad_lines='warn')

    if len(column_names) == df.shape[1]:
        df.columns = column_names
    else:
        df.columns = [f"Col_{i+1}" for i in range(df.shape[1])]

    return df


def create_hydrosphere_plot(df, planet_name=""):
    """
    Create interactive hydrosphere properties plot.

    Args:
        df: DataFrame with profile data
        planet_name: Name of the planet for title

    Returns:
        plotly Figure object
    """
    if not PLOTLY_AVAILABLE:
        return None

    fig = make_subplots(
        rows=3, cols=2,
        subplot_titles=('Temperature vs Depth', 'Density vs Pressure',
                       'Sound Speeds vs Depth', 'Electrical Conductivity vs Depth',
                       'Pressure vs Depth', 'Gravity vs Radius'),
        specs=[[{}, {}], [{}, {}], [{}, {}]]
    )

    # Convert radius to depth (assuming data is from center outward)
    if 'r (m)' in df.columns:
        max_r = df['r (m)'].max()
        df['depth (km)'] = (max_r - df['r (m)']) / 1000

    # Temperature vs Depth
    if 'T (K)' in df.columns and 'depth (km)' in df.columns:
        fig.add_trace(
            go.Scatter(x=df['T (K)'], y=df['depth (km)'], mode='lines', name='Temperature'),
            row=1, col=1
        )
        fig.update_xaxes(title_text="Temperature (K)", row=1, col=1)
        fig.update_yaxes(title_text="Depth (km)", autorange='reversed', row=1, col=1)

    # Density vs Pressure
    if 'rho (kg/m3)' in df.columns and 'P (MPa)' in df.columns:
        fig.add_trace(
            go.Scatter(x=df['P (MPa)'], y=df['rho (kg/m3)'], mode='lines', name='Density'),
            row=1, col=2
        )
        fig.update_xaxes(title_text="Pressure (MPa)", row=1, col=2)
        fig.update_yaxes(title_text="Density (kg/m³)", row=1, col=2)

    # Sound Speeds vs Depth
    if 'VP (km/s)' in df.columns and 'depth (km)' in df.columns:
        fig.add_trace(
            go.Scatter(x=df['VP (km/s)'], y=df['depth (km)'], mode='lines', name='Vp'),
            row=2, col=1
        )
    if 'VS (km/s)' in df.columns and 'depth (km)' in df.columns:
        fig.add_trace(
            go.Scatter(x=df['VS (km/s)'], y=df['depth (km)'], mode='lines', name='Vs'),
            row=2, col=1
        )
    fig.update_xaxes(title_text="Sound Speed (km/s)", row=2, col=1)
    fig.update_yaxes(title_text="Depth (km)", autorange='reversed', row=2, col=1)

    # Electrical Conductivity vs Depth
    if 'sigma (S/m)' in df.columns and 'depth (km)' in df.columns:
        fig.add_trace(
            go.Scatter(x=df['sigma (S/m)'], y=df['depth (km)'], mode='lines', name='Conductivity'),
            row=2, col=2
        )
        fig.update_xaxes(title_text="Conductivity (S/m)", type="log", row=2, col=2)
        fig.update_yaxes(title_text="Depth (km)", autorange='reversed', row=2, col=2)

    # Pressure vs Depth
    if 'P (MPa)' in df.columns and 'depth (km)' in df.columns:
        fig.add_trace(
            go.Scatter(x=df['P (MPa)'], y=df['depth (km)'], mode='lines', name='Pressure'),
            row=3, col=1
        )
        fig.update_xaxes(title_text="Pressure (MPa)", row=3, col=1)
        fig.update_yaxes(title_text="Depth (km)", autorange='reversed', row=3, col=1)

    # Gravity vs Radius
    if 'g (m/s2)' in df.columns and 'r (m)' in df.columns:
        fig.add_trace(
            go.Scatter(x=df['r (m)']/1000, y=df['g (m/s2)'], mode='lines', name='Gravity'),
            row=3, col=2
        )
        fig.update_xaxes(title_text="Radius (km)", row=3, col=2)
        fig.update_yaxes(title_text="Gravity (m/s²)", row=3, col=2)

    fig.update_layout(
        height=900,
        showlegend=True,
        title_text=f"{planet_name} Hydrosphere Properties" if planet_name else "Hydrosphere Properties"
    )

    return fig


def create_seismic_plot(df, planet_name=""):
    """
    Create interactive seismic properties plot.

    Args:
        df: DataFrame with profile data
        planet_name: Name of the planet for title

    Returns:
        plotly Figure object
    """
    if not PLOTLY_AVAILABLE:
        return None

    fig = make_subplots(
        rows=2, cols=2,
        subplot_titles=('Bulk & Shear Moduli', 'Temperature, Pressure, Density',
                       'Sound Speeds', 'Seismic Quality Factor'),
    )

    r_km = df['r (m)'] / 1000 if 'r (m)' in df.columns else None

    # Bulk & Shear Moduli
    if r_km is not None and 'KS (GPa)' in df.columns:
        fig.add_trace(
            go.Scatter(x=r_km, y=df['KS (GPa)'], mode='lines', name='Ks'),
            row=1, col=1
        )
    if r_km is not None and 'GS (GPa)' in df.columns:
        fig.add_trace(
            go.Scatter(x=r_km, y=df['GS (GPa)'], mode='lines', name='Gs'),
            row=1, col=1
        )
    fig.update_xaxes(title_text="Radius (km)", row=1, col=1)
    fig.update_yaxes(title_text="Moduli (GPa)", row=1, col=1)

    # Temperature, Pressure, Density
    if r_km is not None and 'T (K)' in df.columns:
        fig.add_trace(
            go.Scatter(x=r_km, y=df['T (K)'], mode='lines', name='T'),
            row=1, col=2
        )
    fig.update_xaxes(title_text="Radius (km)", row=1, col=2)
    fig.update_yaxes(title_text="Temperature (K)", row=1, col=2)

    # Sound Speeds
    if r_km is not None and 'VP (km/s)' in df.columns:
        fig.add_trace(
            go.Scatter(x=r_km, y=df['VP (km/s)'], mode='lines', name='Vp'),
            row=2, col=1
        )
    if r_km is not None and 'VS (km/s)' in df.columns:
        fig.add_trace(
            go.Scatter(x=r_km, y=df['VS (km/s)'], mode='lines', name='Vs'),
            row=2, col=1
        )
    fig.update_xaxes(title_text="Radius (km)", row=2, col=1)
    fig.update_yaxes(title_text="Sound Speed (km/s)", row=2, col=1)

    # Seismic Quality Factor
    if r_km is not None and 'QS' in df.columns:
        fig.add_trace(
            go.Scatter(x=r_km, y=df['QS'], mode='lines', name='Qs'),
            row=2, col=2
        )
    fig.update_xaxes(title_text="Radius (km)", row=2, col=2)
    fig.update_yaxes(title_text="Quality Factor Qs", type="log", row=2, col=2)

    fig.update_layout(
        height=800,
        showlegend=True,
        title_text=f"{planet_name} Seismic Properties" if planet_name else "Seismic Properties"
    )

    return fig


def create_custom_plot(df, x_col, y_cols, title="", x_label="", y_label="", log_x=False, log_y=False):
    """
    Create a custom interactive plot from data.

    Args:
        df: DataFrame with data
        x_col: Column name for x-axis
        y_cols: List of column names for y-axis
        title: Plot title
        x_label: X-axis label
        y_label: Y-axis label
        log_x: Use log scale for x-axis
        log_y: Use log scale for y-axis

    Returns:
        plotly Figure object
    """
    if not PLOTLY_AVAILABLE:
        return None

    fig = go.Figure()

    for y_col in y_cols:
        if y_col in df.columns and x_col in df.columns:
            fig.add_trace(
                go.Scatter(x=df[x_col], y=df[y_col], mode='lines+markers', name=y_col)
            )

    fig.update_layout(
        title=title,
        xaxis_title=x_label or x_col,
        yaxis_title=y_label or ", ".join(y_cols),
        hovermode='x unified',
        height=600
    )

    if log_x:
        fig.update_xaxes(type="log")
    if log_y:
        fig.update_yaxes(type="log")

    return fig


def save_plotly_figure(fig, filename, format='html'):
    """
    Save a Plotly figure to file.

    Args:
        fig: Plotly Figure object
        filename: Output filename
        format: Output format ('html', 'png', 'pdf', 'svg')

    Returns:
        Success boolean
    """
    if not PLOTLY_AVAILABLE or fig is None:
        return False

    try:
        if format == 'html':
            fig.write_html(filename)
        elif format == 'png':
            fig.write_image(filename, format='png')
        elif format == 'pdf':
            fig.write_image(filename, format='pdf')
        elif format == 'svg':
            fig.write_image(filename, format='svg')
        return True
    except Exception as e:
        st.error(f"Failed to save figure: {e}")
        return False
