from PlanetProfile.Utilities.ResultsStructs import ExplorationResultsStruct
from PlanetProfile.Utilities.defineStructs import ParamsStruct, FigureFilesSubstruct
from PlanetProfile.Main import WriteProfile, ReloadProfile, ExploreOgram
from scipy.interpolate import RegularGridInterpolator
import numpy as np
from PlanetProfile.GetConfig import Color, Style, FigLbl, FigSize, FigMisc
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import os
from copy import deepcopy
import matplotlib.patches as patches
import numpy as np
from scipy.io import loadmat
import pickle
from matplotlib.gridspec import GridSpec
from PlanetProfile.Utilities.DataManip import smoothGrid
from PlanetProfile.Plotting.EssentialHelpers import *
from PlanetProfile.Utilities.SetupInit import SetupFilenames
import logging
logger = logging.getLogger('PlanetProfile')
def InvertBestPlanet(TruePlanet, Params, ExplorationGrid, otherExplorationList=[], interpolateGrid = False, saveGrid = False):
    """ExplorationList = []
    if len(fNames) == 1:
        Exploration, Params = ExploreOgram(TruePlanet.bodyname, Params, fNameOverride=fNames[0])
        ExplorationList.append(Exploration)
    else:
        for fName in fNames:
            Exploration, Params = ExploreOgram(TruePlanet.bodyname, Params, fNameOverride=fName)
            ExplorationList.append(Exploration)"""
    _, DataFiles, _ = SetupFilenames(TruePlanet, Params)
    Params.Inversion.InductionResponseUncertainty_nT = 3.0
    #Params.Inversion.kLoveAmpUncertainity = 0.035
    baseName = f'{DataFiles.saveBase}'
    comparePath = os.path.join(TruePlanet.bodyname, 'figures')
    FigureFiles = FigureFilesSubstruct(comparePath, baseName, FigMisc.xtn)
    Params.FigureFiles = FigureFiles
    Params = setTruePlanetData(TruePlanet, Params)
    explorationListPkl = os.path.join(TruePlanet.bodyname, 'explorationList.pkl')
    if os.path.exists(explorationListPkl) and not saveGrid:
        ExplorationList = pickle.load(open(explorationListPkl, 'rb'))
        Params.Inversion.GRID_RESOLUTION_SCALE = 1
    else:
        if interpolateGrid:
            ExplorationGrid = InterpolateExplorationData(TruePlanet, ExplorationGrid, Params)
        ExplorationList = list(ExplorationGrid.flatten()) + otherExplorationList
        Params.Inversion.GRID_RESOLUTION_SCALE = 5
        
    ExplorationList, Params.Inversion = CalcGridWithinUncertainty(TruePlanet, ExplorationList, Params, interpolateGrid=interpolateGrid)
    if saveGrid:
        pickle.dump(ExplorationList, open(explorationListPkl, 'wb'))
    Params.Inversion.GRID_RESOLUTION_SCALE = 1
    InvertBestPlanetPlot(TruePlanet, ExplorationList, Params)
    return ExplorationGrid

def InterpolateExplorationData(TruePlanet, ExplorationGrid, Params: ParamsStruct):
    nA, nB = ExplorationGrid.shape
    nx = ExplorationGrid[0, 0].nx
    ny = ExplorationGrid[0, 0].ny

    # If you truly just want to interpolate across the model-index axes:
    a_grid = np.arange(nA, dtype=float)
    b_grid = np.arange(nB, dtype=float)

    attrs = ['base.kLoveAmp', 'base.hLoveAmp', 'induction.Bi1Tot_nT']

    interpolators = {}

    for attr in attrs:
        if 'induction' in attr:
            nExc, _ = countPlottableExcitations(TruePlanet.Magnetic.calcedExc, Params.Induct)
            data = np.full((nA, nB, nExc, nx, ny), np.nan, dtype=np.complex64)
        else:
            data = np.full((nA, nB, nx, ny), np.nan, dtype=np.float32)

        for i in range(nA):
            for j in range(nB):
                Exploration = ExplorationGrid[i, j]
                # Extract the attribute value using nested getattr for dot notation
                attr_parts = attr.split('.')
                value = Exploration
                for part in attr_parts:
                    value = getattr(value, part)
                if 'induction' in attr:
                    explorationInductionIndices, planetInductionIndices = inductionIndices(TruePlanet, Exploration, Params)
                    data[i, j, explorationInductionIndices, :, :] = value
                else:
                    data[i, j, :, :] = value

        # If complex interpolation gives trouble, split real/imag (see below).
        interpolators[attr] = RegularGridInterpolator(
            points=(a_grid, b_grid),
            values=data,
            method="linear",
            bounds_error=True,
            fill_value=None
        )
    
    newnA = 101
    newnB = 41
    newa_grid = np.linspace(0, nA - 1, newnA)
    newb_grid = np.linspace(0, nB - 1, newnB)
    ExplorationGridNew = np.empty((newnA, newnB), dtype=object)
    for i in range(newnA):
        aInterpolator = newa_grid[i]
        for j in range(newnB):
            bInterpolator = newb_grid[j]
            ExplorationGridNew[i, j] = ExplorationResultsStruct()
            ExplorationGridNew[i, j].nx = nx
            ExplorationGridNew[i, j].ny = ny
            ExplorationGridNew[i, j].xName = ExplorationGrid[0, 0].xName
            ExplorationGridNew[i, j].yName = ExplorationGrid[0, 0].yName
            ExplorationGridNew[i, j].base.__dict__[ExplorationGrid[0, 0].xName] = ExplorationGrid[0, 0].base.__dict__[ExplorationGrid[0, 0].xName]
            ExplorationGridNew[i, j].base.__dict__[ExplorationGrid[0, 0].yName] = ExplorationGrid[0, 0].base.__dict__[ExplorationGrid[0, 0].yName]
            for attr in attrs:
                if 'base' in attr:
                    ExplorationGridNew[i, j].base.__setattr__(attr.split('.')[-1], interpolators[attr]((aInterpolator, bInterpolator)))
                elif 'induction' in attr:
                    nExc, nExcNames = countPlottableExcitations(TruePlanet.Magnetic.calcedExc, Params.Induct)
                    ExplorationGridNew[i, j].induction.calcedExc = nExcNames
                    ExplorationGridNew[i, j].induction.__setattr__(attr.split('.')[-1], interpolators[attr]((aInterpolator, bInterpolator)))
    return ExplorationGridNew


def InvertBestPlanetPlot(BestPlanet, ExplorationList, Params: ParamsStruct):
    """
    Create uncertainty plots for a single best-fit planet in a single figure.
    """
    FigLbl.SetInversion(BestPlanet.bodyname, ExplorationList[0].xName, ExplorationList[0].yName)
    # Plot settings
    lineStyles = {'kLove': 'solid', 'hLove': 'dashed', 'orbital': 'dashed', 'synodic': 'solid', 
                  'gravity': 'dashed', 'induction': 'solid'}
    colors = {
                'kLove':    '#1f78b4',   # strong blue
        'hLove':    '#33a0c2',   # blue–cyan (same brightness)

        'synodic':  '#e31a1c',   # strong red
        'orbital':  '#ff3a3c',   # red–orange (same brightness)

        'gravity':  '#2e9ab0',   # blend of blue pair (equal strength)
        'induction':'#f4461f'    # blend of red pair (equal strength)
    }

    gravityLabel = f'Gravity Inversion'
    inductionLabel = f'Magnetic Inversion'
    linewidth = 3

    
    # Number of rows depends on the observations we have
    n_rows = 0
    n_rows += 1 if Params.Inversion.INVERT_GRAVITY else 0
    n_rows += 1 if Params.Inversion.INVERT_INDUCTION else 0
    n_rows += 1 if Params.Inversion.INVERT_JOINT else 0
    
    # Calculate figure size with scaling (same as PlotExploreOgramMultiSubplot)
    legendFontSize = 20
    base_size = (12, 9)  # TODO Implement figure size
    fig, axes = plt.subplots(n_rows, 1, figsize=base_size, constrained_layout=True, squeeze=False)
    plot_data = extract_and_validate_plot_data(ExplorationList[0], ExplorationList[0].xName, ExplorationList[0].yName, 
                                               custom_x_axis=FigLbl.xCustomAxis, custom_y_axis=FigLbl.yCustomAxis)
    xData = plot_data['x'].reshape(plot_data['original_shape'])
    yData = plot_data['y'].reshape(plot_data['original_shape'])
    xData, yData, _ = smoothGrid(xData, yData, [xData, yData], Params.Inversion.GRID_RESOLUTION_SCALE)
    allLegendHandles = []
    currentAxIndex = 0
    if Params.Inversion.INVERT_GRAVITY:
        gravityLegendHandles = []
        ax = axes[currentAxIndex, 0]
         # plot h number
        lbl = FigLbl.axisLabelsExplore['hLoveAmp']
        ax.contourf(xData, yData, Params.Inversion.gridWithinhLoveAmpUncertainty, levels=[0.5, 1.5], colors=[colors['hLove']], alpha=0.5, linestyles=lineStyles['hLove'])
        ax.contour(xData, yData, Params.Inversion.gridWithinhLoveAmpUncertainty, levels=[0.5], colors=[colors['hLove']], alpha=0.5, linestyles=lineStyles['hLove'], linewidths=linewidth)
        # Proxy line for legend
        gravityLegendHandles.append(
            mlines.Line2D([], [], color=colors['hLove'], linestyle=lineStyles['hLove'], label=lbl)
        )
        # Plot k love
        lbl = FigLbl.axisLabelsExplore['kLoveAmp']
        ax.contourf(xData, yData, Params.Inversion.gridWithinkLoveAmpUncertainty, levels=[0.5, 1.5], colors=[colors['kLove']], alpha=0.5, linestyles=lineStyles['kLove'])
        ax.contour(xData, yData, Params.Inversion.gridWithinkLoveAmpUncertainty, levels=[0.5], colors=[colors['kLove']], alpha=0.5, linestyles=lineStyles['kLove'], linewidths=linewidth)
        
        # Proxy line for legend
        gravityLegendHandles.append(
            mlines.Line2D([], [], color=colors['kLove'], linestyle=lineStyles['kLove'], label=lbl)
        )
       
        
        allLegendHandles.append(gravityLegendHandles)
        currentAxIndex += 1

    if Params.Inversion.INVERT_INDUCTION:
        inductionLegendHandles = []
        ax = axes[currentAxIndex, 0]
        nExc, excNames = countPlottableExcitations(ExplorationList[0].induction.calcedExc, Params.Induct)
        for iExc, excName in zip(range(nExc), excNames):
            lbl = f'{FigLbl.Bi1Tot_nTLabel} ({excName})'
            ax.contourf(xData, yData, Params.Inversion.gridWithinBi1Tot_nTUncertainty[iExc, :, :], levels=[0.5, 1.5], colors=[colors[excName]], alpha=0.5, linestyles=lineStyles[excName])
            ax.contour(xData, yData, Params.Inversion.gridWithinBi1Tot_nTUncertainty[iExc, :, :], levels=[0.5], colors=[colors[excName]], alpha=0.5, linestyles=lineStyles[excName], linewidths=linewidth)
            # Proxy line for legend
            inductionLegendHandles.append(
                mlines.Line2D([], [], color=colors[excName], linestyle=lineStyles[excName], label=lbl)
            )
        allLegendHandles.append(inductionLegendHandles)
        currentAxIndex += 1
    
    if Params.Inversion.INVERT_JOINT:
        ax = axes[currentAxIndex, 0]
        jointInversionLegendHandles = []
        ax.contourf(xData, yData, Params.Inversion.gridWithinGravityUncertainty, levels=[0.5, 1.5], colors=[colors['gravity']], alpha=0.5, linestyles=lineStyles['gravity'])
        ax.contour(xData, yData, Params.Inversion.gridWithinGravityUncertainty, levels=[0.5], colors=[colors['gravity']], alpha=0.5, linestyles=lineStyles['gravity'], linewidths=linewidth)
        # Proxy line for legend
        jointInversionLegendHandles.append(
            mlines.Line2D([], [], color=colors['gravity'], linestyle=lineStyles['gravity'], label=gravityLabel)
        )
        
        ax.contourf(xData, yData, Params.Inversion.gridWithinInductionUncertainty, levels=[0.5, 1.5], colors=[colors['induction']], alpha=0.5, linestyles=lineStyles['induction'])
        ax.contour(xData, yData, Params.Inversion.gridWithinInductionUncertainty, levels=[0.5], colors=[colors['induction']], alpha=0.5, linestyles=lineStyles['induction'])
        # Proxy line for legend
        jointInversionLegendHandles.append(
            mlines.Line2D([], [], color=colors['induction'], linestyle=lineStyles['induction'], label=inductionLabel)
        )
        allLegendHandles.append(jointInversionLegendHandles)

    # Set axis labels and scales and plot best planet point
    for i, ax in enumerate(axes.flatten()):
        ax.legend(handles=allLegendHandles[i], loc='upper right', fontsize=legendFontSize)
        # Set tick label size
        ax.tick_params(axis='both', labelsize=Style.TS_ticks)
        # Set y axis
        ax.set_yscale(FigLbl.yScaleExplore)
        ax.set_ylim(FigLbl.yCustomAxis)
        # Plot best planet point
        ax.scatter(Params.Inversion.xTrueFit, Params.Inversion.yTrueFit, color='black', marker='x', s=100)
        ax.set_xlabel(FigLbl.xLabelExplore, fontsize=18)
        ax.set_ylabel(FigLbl.yLabelExplore, fontsize=18)
        """if i < len(axes.flatten()) - 1:
            ax.tick_params(axis='x', labelbottom=False)"""
    
    #fig.suptitle(FigLbl.inversionTitle, verticalalignment='top', horizontalalignment='center', fontsize = Style.TS_super)
    # Set figure y label
    #fig.supxlabel(FigLbl.xLabelExplore, verticalalignment='bottom', horizontalalignment='center', fontsize = Style.TS_axis)
    # Save to uncertainty multiplot file
    fig.savefig(Params.FigureFiles.inversion, format=FigMisc.figFormat, dpi=FigMisc.dpi, bbox_inches='tight', transparent=True)
    log.debug(f'Uncertainty multiplot saved to file: {Params.FigureFiles.inversion}')
    
    plt.close()

def setTruePlanetData(TruePlanet, Params: ParamsStruct):
    """
    Set the true planet data for the inversion.
    """
    if Params.Explore.xName in Params.Explore.provideExploreRange:
        xDataList = loadmat(Params.DataFiles.xRangeData)['Data'].flatten().tolist()
        xDataList = [s.strip() if isinstance(s, str) else s for s in xDataList]
        xInputList = np.linspace(Params.Explore.xRange[0], Params.Explore.xRange[1], Params.Explore.nx)
    else:
        xDataList, xInputList = [], []
    if Params.Explore.yName in Params.Explore.provideExploreRange:
        yDataList = loadmat(Params.DataFiles.yRangeData)['Data'].flatten().tolist()
        yDataList = [s.strip() if isinstance(s, str) else s for s in yDataList]
        yInputList = np.linspace(Params.Explore.yRange[0], Params.Explore.yRange[1], Params.Explore.ny)
    else:
        yDataList, yInputList = [], []
    Params.Inversion.xTrueFit = getTruePlanetData(TruePlanet, Params.Inversion.invertXName, dataList = xDataList, inputList = xInputList)
    Params.Inversion.yTrueFit = getTruePlanetData(TruePlanet, Params.Inversion.invertYName, dataList = yDataList, inputList = yInputList)
    return Params

def inductionIndices(bestPlanet, Exploration, Params: ParamsStruct):
    nExcPlanet, nExcNamesPlanet = countPlottableExcitations(bestPlanet.Magnetic.calcedExc, Params.Induct)
    explorationInductionIndices = []
    planetInductionIndices = []
    for i, planet_exc_name in enumerate(nExcNamesPlanet):
        if planet_exc_name in Exploration.induction.calcedExc:
            exploration_index = Exploration.induction.calcedExc.index(planet_exc_name)
            explorationInductionIndices.append(exploration_index)
            planetInductionIndices.append(i)
    return explorationInductionIndices, planetInductionIndices
def CalcGridWithinUncertainty(bestPlanet, ExplorationList, Params: ParamsStruct, interpolateGrid = False):
    """
    Calculate the grid of models within the uncertainty of the best-fit model.
    """ 
    grid_nx = ExplorationList[0].nx * Params.Inversion.GRID_RESOLUTION_SCALE
    grid_ny = ExplorationList[0].ny * Params.Inversion.GRID_RESOLUTION_SCALE
    # Get 
    nExcPlanet, nExcNamesPlanet = countPlottableExcitations(bestPlanet.Magnetic.calcedExc, Params.Induct)
    
    # Create grid
    Params.Inversion.gridWithinkLoveAmpUncertainty = np.zeros((grid_nx, grid_ny), dtype=bool)
    Params.Inversion.gridWithinhLoveAmpUncertainty = np.zeros((grid_nx, grid_ny), dtype=bool)
    Params.Inversion.gridWithinBi1Tot_nTUncertainty = np.zeros((nExcPlanet, grid_nx, grid_ny), dtype=bool)
    for i, Exploration in enumerate(ExplorationList): 
        print(i)
        ExplorationList[i].nx = grid_nx
        ExplorationList[i].ny = grid_ny
        plot_data = extract_and_validate_plot_data(Exploration, Exploration.xName, Exploration.yName, 
                                               custom_x_axis=FigLbl.xCustomAxis, custom_y_axis=FigLbl.yCustomAxis)
        xData = plot_data['x'].reshape(plot_data['original_shape'])
        yData = plot_data['y'].reshape(plot_data['original_shape'])
        if interpolateGrid:
            xFine, yFine, [kLoveAmp, hLoveAmp] = smoothGrid(xData, yData, [Exploration.base.kLoveAmp, Exploration.base.hLoveAmp], Params.Inversion.GRID_RESOLUTION_SCALE)
            ExplorationList[i].base.kLoveAmp = kLoveAmp
            ExplorationList[i].base.hLoveAmp = hLoveAmp
            ExplorationList[i].xData = xFine
            ExplorationList[i].yData = yFine
            ExplorationList[i].base.__setattr__(Exploration.xName, xFine)
            ExplorationList[i].base.__setattr__(Exploration.yName, yFine)
        else:
            kLoveAmp = Exploration.base.kLoveAmp
            hLoveAmp = Exploration.base.hLoveAmp
        Params.Inversion.gridWithinkLoveAmpUncertainty = np.logical_or(Params.Inversion.gridWithinkLoveAmpUncertainty, withinUncertainty(bestPlanet.Gravity.kAmp, kLoveAmp, Params.Inversion.kLoveAmpUncertainity))
        Params.Inversion.gridWithinhLoveAmpUncertainty = np.logical_or(Params.Inversion.gridWithinhLoveAmpUncertainty, withinUncertainty(bestPlanet.Gravity.hAmp, hLoveAmp, Params.Inversion.hLoveAmpUncertainity))
        
        # Find matching excitation names and their indices
        explorationInductionIndices, planetInductionIndices = inductionIndices(bestPlanet, Exploration, Params)
        
        # Get the corresponding induction data from the exploration
        smoothedInductionData = np.zeros((len(explorationInductionIndices), grid_nx, grid_ny), dtype=np.complex128) * np.nan
        for j, exploration_index in enumerate(explorationInductionIndices):
            if interpolateGrid:
                _, _, [smoothed_induction] = smoothGrid(xData, yData, [Exploration.induction.Bi1Tot_nT[exploration_index, :, :]], Params.Inversion.GRID_RESOLUTION_SCALE)
            else:
                smoothed_induction = Exploration.induction.Bi1Tot_nT[exploration_index, :, :]
            smoothedInductionData[j, :, :] = smoothed_induction
        ExplorationList[i].induction.Bi1Tot_nT = smoothedInductionData
        ExplorationList[i].induction.calcedExc = bestPlanet.Magnetic.calcedExc
        inductionWithinUncertaintyGrid = withinUncertainty(bestPlanet.Magnetic.Bi1Tot_nT[planetInductionIndices], smoothedInductionData, Params.Inversion.InductionResponseUncertainty_nT)
        
        Params.Inversion.gridWithinBi1Tot_nTUncertainty[planetInductionIndices, :, :] = np.logical_or(Params.Inversion.gridWithinBi1Tot_nTUncertainty[planetInductionIndices, :, :], inductionWithinUncertaintyGrid)
    
    # Get summary results
    Params.Inversion.gridWithinGravityUncertainty = np.logical_and(Params.Inversion.gridWithinkLoveAmpUncertainty, Params.Inversion.gridWithinhLoveAmpUncertainty)
    Params.Inversion.gridWithinInductionUncertainty = np.all(Params.Inversion.gridWithinBi1Tot_nTUncertainty, axis=0)
    return ExplorationList, Params.Inversion

  
def withinUncertainty(truePlanetData, GridData, UncertaintyData):
    """
    Calculate the grid of models within the uncertainty of the best-fit model.
    For complex data, uses magnitude of difference (circular uncertainty region).
    For real data, uses direct comparison (linear uncertainty bounds).
    """
    # Check if data is complex
    if np.iscomplexobj(truePlanetData):
        UncertaintyData = np.array([UncertaintyData, UncertaintyData*2])
        # For complex data, check if magnitude of difference is within uncertainty
        # This creates a circular uncertainty region in the complex plane
        # GridData shape: (nexc, planetGridWidth, planetGridHeight)
        # bestPlanetData shape: (nexc,)
        # Need to broadcast bestPlanetData to match GridData dimensions
        realGridData = np.real(GridData)
        imagGridData = np.imag(GridData)
        rTruePlanetData = np.real(truePlanetData)
        iTruePlanetData = np.imag(truePlanetData)

        if realGridData.ndim == 3:
            # Expand true planet data to match grid dimensions
            rTruePlanetData = rTruePlanetData[:, np.newaxis, np.newaxis]
            iTruePlanetData = iTruePlanetData[:, np.newaxis, np.newaxis]

        realDifference = np.abs(realGridData - rTruePlanetData)
        imagDifference = np.abs(imagGridData - iTruePlanetData)
        realWithinUncertainty = realDifference <= UncertaintyData[:, np.newaxis, np.newaxis]
        imagWithinUncertainty = imagDifference <= UncertaintyData[:, np.newaxis, np.newaxis]
        withinUncertainty = np.logical_and(realWithinUncertainty, imagWithinUncertainty)
    else:
        # For real data, use direct comparison (linear bounds)
        upper_bound = truePlanetData + UncertaintyData
        lower_bound = truePlanetData - UncertaintyData
        withinUncertainty = (GridData >= lower_bound) & (GridData <= upper_bound)
    return withinUncertainty

def getTruePlanetData(TruePlanet, name, dataList = None, inputList = None):
    """
    Get the true planet data for a given variable name.
    """
    if name == 'oceanComp':
        return inputList[dataList.index(TruePlanet.Ocean.comp)]
    elif name == 'zb_approximate_km':
        return TruePlanet.Bulk.zb_approximate_km