"""
Amortized Inference (SBI) Tab for PlanetProfileApp

Instant Bayesian inference using pre-trained Simulation-Based Inference (SBI)
models (Normalizing Flows).

Author: PlanetProfile Team
Date: 2026-06-04
"""
import streamlit as st
import os
import sys
import numpy as np

# Path setup
this_file = os.path.abspath(__file__)
pages_directory = os.path.dirname(this_file)
app_directory = os.path.dirname(pages_directory)
parent_directory = os.path.dirname(app_directory)

if app_directory not in sys.path:
    sys.path.insert(0, app_directory)
if parent_directory not in sys.path:
    sys.path.insert(0, parent_directory)

from Utilities.get_planet import get_planet
from Utilities.planet_sidebar import show_planet_status

# ============================================================
# Artifact registry
#
# Maps the UI "Select pre-trained model" labels to the trained SBI
# posterior estimator file expected under
# PlanetProfile/Inference/sbi_artifacts/ (see INDEX.md in that directory
# for the naming convention). Resolved relative to the repo root
# (`parent_directory`, set up above the same way every other
# PlanetProfileApp page locates PlanetProfile/ paths).
# ============================================================
_SBI_ARTIFACTS_DIR = os.path.join(
    parent_directory, "PlanetProfile", "Inference", "sbi_artifacts"
)

_ARTIFACTS = {
    "Titan (Andrade, No-Ocean)": "titan_andrade_noocean_posterior.pt",
    "Titan (Maxwell, Ocean)": "titan_maxwell_ocean_posterior.pt",
    "Europa (Andrade, Thin Shell)": "europa_andrade_thinshell_posterior.pt",
}

# Fixed seed for UI reproducibility: repeated clicks with the same inputs
# reproduce the same posterior draw (not used for any scientific claim).
_SBI_UI_SEED = 20260704


@st.cache_resource(show_spinner="Loading SBI posterior artifact...")
def _load_sbi_artifact(artifact_path: str):
    """Load a trained SBIRunner artifact, cached per artifact path.

    Streamlit's cache_resource keys on the function arguments, so switching
    between model slots (different `artifact_path`) loads/caches each
    artifact independently within the session.
    """
    from PlanetProfile.Inference.sbi_runner import SBIRunner
    return SBIRunner.load_artifact(artifact_path)


def main():
    show_planet_status()
    st.title("🚀 Amortized Inference (SBI)")

    st.markdown("""
    **Simulation-Based Inference (SBI)** allows for near-instant posterior estimation
    by leveraging deep learning models (Normalizing Flows) pre-trained on millions
    of interior structure simulations.
    """)

    st.subheader("⚙️ Model Selection")
    model_choice = st.selectbox(
        "Select pre-trained model:",
        options=list(_ARTIFACTS.keys()),
        index=0
    )

    artifact_path = os.path.join(_SBI_ARTIFACTS_DIR, _ARTIFACTS[model_choice])
    artifact_ready = os.path.exists(artifact_path)

    if not artifact_ready:
        st.info("💡 **Status:** This tab is currently under development. Pre-trained models for Titan and Europa are being generated using the MCMC pipeline.")

    st.subheader("🎯 Observables")
    st.write("Enter the observed data to get an instant posterior estimate.")

    col1, col2 = st.columns(2)
    with col1:
        re_k2 = st.number_input("Re(k₂) — Real part:", value=0.608, format="%.4f")
    with col2:
        im_k2 = st.number_input("Im(k₂) — Imaginary part:", value=0.135, format="%.4f")

    run_clicked = st.button(
        "⚡ Generate Instant Posterior",
        disabled=not artifact_ready,
        help=None if artifact_ready else "No trained artifact found for this model slot yet.",
    )

    st.divider()

    if run_clicked and artifact_ready:
        try:
            artifact = _load_sbi_artifact(artifact_path)

            # imag_convention for all shipped artifacts is 'abs' (|Im k2|),
            # matching the MCMC likelihood/GUI convention — see
            # plans/monte-carlo-sbi-implementation-plan.md. Fail loudly on any
            # artifact trained under a different convention rather than
            # silently mis-conditioning the posterior.
            convention = getattr(artifact, 'imag_convention', 'abs')
            if convention != 'abs':
                raise ValueError(
                    f"Artifact was trained with imag_convention='{convention}', "
                    f"but this page conditions on |Im(k2)| ('abs'). "
                    f"Retrain the artifact with the 'abs' convention."
                )
            x_obs = {'Re_k2': float(re_k2), 'abs_Im_k2': float(abs(im_k2))}

            with st.spinner("Sampling posterior..."):
                samples = artifact.sample_posterior(
                    x_obs, n_samples=10000, seed=_SBI_UI_SEED
                )
            samples = np.asarray(samples)

            param_names = list(getattr(artifact, 'param_names', []))
            param_labels = list(getattr(artifact, 'param_labels', param_names))
            param_units = list(getattr(artifact, 'param_units', ['' for _ in param_names]))

            if not param_names:
                # Fall back to generic column names if the artifact does not
                # expose parameter metadata.
                param_names = [f'theta_{i}' for i in range(samples.shape[1])]
                param_labels = param_names
                param_units = ['' for _ in param_names]

            st.success(f"Sampled {samples.shape[0]} posterior draws from '{model_choice}'.")

            # Summary table: median + 16/84 percentiles per parameter, with units.
            with st.expander("📋 Parameter Summary", expanded=True):
                summary_data = []
                for i, name in enumerate(param_names):
                    col = samples[:, i]
                    unit = param_units[i] if i < len(param_units) else ''
                    summary_data.append({
                        'Parameter': param_labels[i] if i < len(param_labels) else name,
                        'Units': unit,
                        'Median': f"{np.median(col):.4f}",
                        '16th %ile': f"{np.percentile(col, 16):.4f}",
                        '84th %ile': f"{np.percentile(col, 84):.4f}",
                    })
                st.table(summary_data)

            # Corner plot, following the same inline `corner`-library pattern
            # used by the MCMC Inference tab (PlanetProfileApp/pages/Inference.py)
            # rather than mcmc_plots.plot_corner, which is designed to save a
            # PNG to disk for the batch/diagnostics pipeline. Reusing the GUI's
            # existing inline rendering convention keeps this tab consistent
            # with the rest of the app without adding new plotting
            # infrastructure.
            with st.expander("🔺 Corner Plot", expanded=True):
                try:
                    import corner
                    import matplotlib.pyplot as plt

                    valid_cols = []
                    for i in range(samples.shape[1]):
                        col_data = samples[:, i]
                        col_data = col_data[np.isfinite(col_data)]
                        if len(col_data) > 0 and np.ptp(col_data) > 1e-12:
                            valid_cols.append(i)

                    if not valid_cols:
                        st.warning("No parameters with dynamic range found for corner plot.")
                    else:
                        corner_samples_plot = samples[:, valid_cols]
                        corner_labels_plot = [param_labels[i] for i in valid_cols]

                        n_dim = len(valid_cols)
                        fig_size = max(10, 2.5 * n_dim)
                        fig = plt.figure(figsize=(fig_size, fig_size))
                        corner.corner(
                            corner_samples_plot,
                            labels=corner_labels_plot,
                            quantiles=[0.16, 0.5, 0.84],
                            show_titles=True,
                            title_fmt='.2f',
                            title_kwargs={'fontsize': 10},
                            color='steelblue',
                            fig=fig,
                        )
                        st.pyplot(fig)
                        plt.close(fig)
                except ImportError:
                    st.info("Install the `corner` library to view corner plots: `pip install corner`")
                except Exception as e:
                    st.warning(f"Corner plot unavailable: {e}")

        except Exception as e:
            st.error(f"SBI posterior sampling failed: {e}")

    with st.expander("📚 About SBI"):
        st.markdown("""
        **What is SBI?**

        Traditional MCMC requires running the forward model thousands of times for
        every new observation. SBI "amortizes" this cost by training a neural network
        to learn the mapping from observables back to parameters.

        **Advantages:**
        - **Instant:** Posteriors are generated in milliseconds.
        - **Interactive:** Perfect for "what-if" scenarios in the GUI.
        - **Robust:** Can handle complex likelihoods that are difficult for MCMC.
        """)

if __name__ == "__main__":
    main()
