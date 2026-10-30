FROM condaforge/miniforge3:latest

# poppler: pdf2image previews. texlive set: the PlanetProfile usetex
# preamble (UserConfigs/configPPplots.py) needs siunitx + mhchem
# (texlive-science) AND the stix font + upgreek (texlive-fonts-extra —
# no smaller Debian package carries stix.sty); cm-super + dvipng for
# matplotlib's usetex pipeline.
RUN apt-get update \
 && apt-get install -y --no-install-recommends poppler-utils \
      texlive-latex-base texlive-latex-extra texlive-science \
      texlive-fonts-recommended texlive-fonts-extra cm-super dvipng \
 && rm -rf /var/lib/apt/lists/*

# HF runs Spaces as uid 1000. Ubuntu-Noble bases (miniforge) already ship
# an 'ubuntu' user at uid 1000 — creating another fails the build
# ("UID 1000 is not unique"); ensure one exists and give it a writable HOME.
RUN (id -u 1000 >/dev/null 2>&1 || useradd -m -u 1000 user) \
 && mkdir -p /home/appuser && chown 1000 /home/appuser

# reaktoro from conda-forge (pip cannot install it). A DEDICATED env —
# never mutate miniforge's base (its pinned newer python makes a
# python=3.11 downgrade unresolvable, failing the build).
RUN mamba create -y -n pp python=3.11 reaktoro -c conda-forge \
 && mamba clean -afy
ENV PATH=/opt/conda/envs/pp/bin:$PATH

WORKDIR /app
COPY requirements.txt /app/requirements.txt
RUN pip install --no-cache-dir -r requirements.txt

COPY --chown=1000 . /app
# COPY --chown owns the contents, but the /app directory node itself is
# root's (created by WORKDIR): the app must be able to create dirs in its
# CWD (UserConfigs/, run outputs <Body>/, output/exploreograms/).
RUN chown 1000 /app
USER 1000
ENV HOME=/home/appuser

# Public demo mode (see PlanetProfileApp/Utilities/app_mode.py):
# amortized inference + single capped forward runs; MCMC and new
# exploreogram grids stay disabled.
ENV PP_PUBLIC_MODE=1

EXPOSE 7860
CMD ["streamlit", "run", "PlanetProfileApp/PlanetProfileApp.py", \
     "--server.port=7860", "--server.address=0.0.0.0", "--server.headless=true"]
