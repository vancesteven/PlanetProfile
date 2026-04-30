# LaTeX Integration Checklist

**Target:** `/Users/svance/Library/CloudStorage/Dropbox/planetprofile-genai-manuscript/PlanetProfile-genai/main.tex`  
**Format:** AASTeX 6.31 for Planetary Science Journal  
**Status:** Template ready, awaiting content

---

## Pre-Integration Tasks

### 1. Update Metadata (Lines 9-23)

```latex
\submitjournal{PSJ} % Already correct
\shorttitle{HP Ice Convection in Ocean Worlds}
\shortauthors{Vance et al.}

\title{Implementation and Application of High-Pressure Ice Convection with Partial Melt Prediction in PlanetProfile: Implications for Titan's Interior}

\correspondingauthor{Steven D. Vance}
\email{svance@jpl.nasa.gov}

\author[ORCID]{Steven D. Vance}
\affiliation{Jet Propulsion Laboratory, California Institute of Technology, Pasadena, CA, USA}

% Additional authors from claude-flow team if appropriate
```

### 2. Abstract (Lines 30-32)

Draft abstract summarizing:
- Kalousova & Sotin (2018) model implementation
- Application to Titan with/without ocean
- Key finding: Ocean presence critical for HP ice melting (680× Ra*c difference)
- Clathrate effect: 93% stagnant lid reduction
- Implications for ocean world characterization

Target: 200-250 words

### 3. Keywords (Line 35)

```latex
\keywords{Ocean worlds --- High-pressure ice --- Titan --- Convection --- Phase transitions --- Tidal heating}
```

---

## Section Integration

### §1 Introduction (Lines 37-38)

**Source:** `MANUSCRIPT_DRAFT_HP_ICE_CONVECTION.md` lines 7-36  
**Length:** ~1,200 words  
**Tasks:**
- Convert markdown to LaTeX
- Add citations: \citep{} for parenthetical, \citet{} for textual
- Key citations: Petricca et al. (2025), Kalousova & Sotin (2018), Vance et al. (2018)
- Ensure equation references use \ref{} if equations numbered

### §2 Methods (Lines 40-41)

**Source:** `MANUSCRIPT_METHODS_SECTION.md`  
**Length:** ~2,500 words  
**Tasks:**
- Convert LaTeX equations (already in $$ format)
- Add equation labels: \label{eq:RaStar}, \label{eq:RaCrit}, \label{eq:Ht}, \label{eq:delta}
- Table for HP ice parameters (optional):
  ```latex
  \begin{deluxetable}{lccc}
  \tablehead{\colhead{Phase} & \colhead{$\mu_0$ [Pa s]} & \colhead{$E_a$ [kJ mol$^{-1}$]} & \colhead{Source}}
  \startdata
  Ice Ih & $10^{14}$ & 60 & Goldsby \& Kohlstedt (2001) \\
  ...
  \enddata
  \end{deluxetable}
  ```
- Add figure references for Test34/35 configurations

### §3 Results (Lines 42-43)

**Source:** Results-analyst integration (in progress)  
**Data:** `TITAN_HP_ICE_CONVECTION_RESULTS.md`  
**Length:** ~1,500-2,000 words  
**Tasks:**
- Create comparison table (Ocean vs No-Ocean):
  ```latex
  \begin{deluxetable*}{lcccccc}
  \tablehead{\colhead{Phase} & \colhead{$T_{conv}$ [K]} & \colhead{Ra$^*$} & \colhead{Ra$_c^*$} & \colhead{$e_{lid}$ [km]} & \colhead{$D_{conv}$ [km]} & \colhead{Status}}
  \startdata
  \multicolumn{7}{c}{Test34: Ocean Present (T_b = 255 K)} \\
  Ice Ih & 247.6 & $6.74 \times 10^7$ & $2.59 \times 10^6$ & 48.3 & 88.1 & Supercritical \\
  ...
  \enddata
  \end{deluxetable*}
  ```
- Reference figures: \ref{fig:test34_hydro}, \ref{fig:test35_convdiag}
- Highlight key finding: Only Ice VI supercritical in ocean case

### §4 Discussion (Lines 44-45)

**Source:** `MANUSCRIPT_DRAFT_HP_ICE_CONVECTION.md` lines 63-298  
**Length:** ~3,000 words  
**Tasks:**
- Convert markdown subsections to \subsection{}
  - §4.1: \subsection{Paradigm Shift: Ocean-Free Titan}
  - §4.2: \subsection{Viscosity Interpretations}
  - §4.3: \subsection{Habitability Implications}
  - §4.4: \subsection{Other Ocean Worlds}
- Ensure citations formatted correctly
- Cross-reference Results section: \ref{sec:results}

### §5 Conclusions (Lines 46-47)

**Source:** `MANUSCRIPT_CONCLUSIONS_DRAFT.md` (after reviews)  
**Length:** ~750 words  
**Tasks:**
- Convert 5 numbered findings to itemize or enumerate
- Add future work as separate paragraph
- Strong closing statement

---

## Supplementary Content

### Acknowledgments (Line 48)

```latex
\acknowledgments
This research was carried out at the Jet Propulsion Laboratory, California Institute of Technology, under a contract with the National Aeronautics and Space Administration (NASA). We thank [reviewers/colleagues] for valuable discussions. Computational resources were provided by [institution]. This work benefited from the open-source PlanetProfile framework (github.com/vancesteven/PlanetProfile).
```

### Software (Line 50)

```latex
\software{
    \texttt{Python} \citep{vanrossum2009},
    \texttt{PlanetProfile} \citep{vance2018},
    \texttt{SeaFreeze} \citep{journaux2020},
    \texttt{PerpleX} \citep{connolly2009},
    \texttt{NumPy} \citep{harris2020},
    \texttt{SciPy} \citep{virtanen2020},
    \texttt{Matplotlib} \citep{hunter2007}
}
```

### Bibliography (Line 52)

**File:** `references.bib` (needs to be created)  
**Essential references:**
- Kalousová & Sotin (2018) - Main model
- Petricca et al. (2025) - Motivation
- Vance et al. (2018) - PlanetProfile
- Journaux et al. (2020) - SeaFreeze
- Deschamps & Sotin (2001) - Comparison
- Durham et al. (1996, 1997) - Viscosity data
- Goldsby & Kohlstedt (2001) - Ice Ih rheology
- Choukroun & Sotin (2012) - Clathrate effects

---

## Figure Integration

### Required Figures

**From Test34:**
1. `Test34Profile_*Hydrosphere.pdf` - Ice phases and ocean structure
2. `Test34Profile_*ConvDiag.pdf` - Convection layer breakdown
3. `Test34Profile_*Viscosity.pdf` - Effective viscosity profile

**From Test35:**
4. `Test35Profile_*Hydrosphere.pdf` - No-ocean structure with clathrate
5. `Test35Profile_*ConvDiag.pdf` - Convection comparison

**Figure captions template:**
```latex
\begin{figure}
\centering
\includegraphics[width=\columnwidth]{Test34Profile_PureH2O_Tb255.0K_ConstantInnerRhoHydrosphere.pdf}
\caption{Interior structure for Titan with subsurface ocean (Test34). 
Upper panel shows ice phases (colors) and temperature profile (black line). 
Ice VI layer (purple) at base of hydrosphere forms 21.8 km temperate layer 
(Ra$^*$ = $1.32 \times 10^{10}$ > Ra$_c^*$ = $1.18 \times 10^9$), 
indicated by dashed region. Ocean thickness ~28 km at P = 178 MPa.}
\label{fig:test34_hydro}
\end{figure}
```

---

## Post-Integration Checks

### Compilation
```bash
cd /Users/svance/Library/CloudStorage/Dropbox/planetprofile-genai-manuscript/PlanetProfile-genai
pdflatex main.tex
bibtex main
pdflatex main.tex
pdflatex main.tex
```

### Quality Checks
- [ ] All citations resolved (no [?])
- [ ] All references resolved (no ??)
- [ ] Figures appear in correct locations
- [ ] Tables formatted properly
- [ ] Equations numbered correctly
- [ ] No overfull hboxes (or acceptable)
- [ ] Two-column format correct
- [ ] Page count reasonable (~15-20 pages target)

### Editorial Review
- [ ] Consistent notation throughout
- [ ] No orphan headers
- [ ] Figure/table placement near references
- [ ] Abstract reflects complete findings
- [ ] Keywords comprehensive

---

## Timeline

**Phase 1:** Metadata + Abstract (15 min)  
**Phase 2:** Introduction + Methods (30 min)  
**Phase 3:** Results (20 min, after results-analyst)  
**Phase 4:** Discussion + Conclusions (30 min)  
**Phase 5:** Bibliography + Figures (20 min)  
**Phase 6:** Compilation + fixes (15 min)  
**Phase 7:** Final editorial review (30 min)

**Total:** ~2.5 hours after all dependencies met

---

## Status

- ✅ Template structure ready
- ✅ Section source files identified
- 🔄 Awaiting results-analyst completion
- 🔄 Awaiting science/editor reviews of Conclusions
- 🔲 Ready to begin integration on completion

**Next action:** Start with metadata + abstract while waiting for final sections.
