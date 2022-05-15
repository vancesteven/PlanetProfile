function cfg = config
    % Set CALC_NEW options to 0 to re-use profile data when possible.
    % It is recommended to keep CALC_NEW=1 except when intermediate parameters
    % such as layer thicknesses will not change between runs.
    
    % Runtime options
    cfg.CALC_NEW =       0;
    cfg.TESTING =        0; % For comparing between branches in python conversion
    cfg.CALC_NEW_REF =   1;
    cfg.CALC_NEW_SOUND = 1;
    cfg.CALC_NEW_INDUC = 1;
    cfg.SKIP_PROFILES =  0; % Whether to skip past all PlanetProfile.m plotting and calculations
    cfg.NO_PLOTS =       0;
    cfg.HOLD =           1; % Whether to overlay runs when possible
    cfg.CONDUCT =        1; % Calculate electrical conductivity
    cfg.REDUCED =        1; % Whether to limit number of ocean layers for faster computation of layered induction
    cfg.DISP_LAYERS =    1; % Whether to display layer depths and heat fluxes for user
    cfg.DISP_TABLES =    0; % Whether to print latex-formatted tables to Matlab command line
    cfg.DEPRECATED =     0; % Whether to allow deprecated code to run. Will often cause errors.
    
    % Magnetic induction calculation settings
    cfg.DO_EUR = 1;
    cfg.DO_GAN = 1;
    cfg.DO_CAL = 1;
    cfg.DO_ENC = 1;
    cfg.DO_MIR = 1;
    cfg.DO_ARI = 1;
    cfg.DO_PER = 1; % Convert frequency axes to periods
    cfg.DO_LEGEND = 0;
    cfg.PLOT_PERMEABILITY = 0; % Permeability is not currently implemented, so this should be 0 until it is.
    cfg.PLOT_FFT = 1; % Whether to show plots of fourier space
    cfg.PLOT_CONTOURS = 1; % Contours or surfaces
    cfg.PLOT_V2020S = 1; % Mark the selected ocean/conductivity combos used in Vance et al. 2020
    cfg.intMethod = 'makima'; % Interpolation method. Certain ones can cause wiggles, notably 'linear'.
    cfg.npts_k = 50;
    cfg.npts_D = 60;
    cfg.np_intp = 200;
    cfg.npts_w = 100;
    cfg.np_wfine = 1000;
    cfg.nIntL = 3; % Number of ocean layers to use when cfg.REDUCED = 1
    cfg.opts_odeParams = odeset('RelTol',1e-10,'AbsTol',1e-10,'MaxStep', 2e3,'InitialStep',1e-2);
    cfg.opts_odeLayers = odeset('RelTol',1e-8, 'AbsTol',1e-10,'MaxStep',10e3,'InitialStep',1e-2);
    
    % General figure options
    cfg.dft_font = 'Arial';
    cfg.dft_math = 'STIX Math';
    cfg.interpreter = 'tex';
    cfg.fig_fmt = '-depsc';
    cfg.xtn = '_pp.eps';
    
    % Color selection
    cfg.b000ff = [176,  0,255]/255; % b000ff.
    % Note the second colons are required to retrieve all 3 RGB elements
    % due to Matlab idiosyncracies.
    cfg.cmap = parula(100);
    cc = 100/256;
    cfg.col_contSyn = cfg.cmap(floor(100*cc),:);
    cfg.col_contOrb = cfg.cmap(floor( 10*cc),:);
    cfg.col_contHrm = cfg.cmap(floor(200*cc),:);
    
    cfg.col_Sw = summer(200); cfg.col_Sw = cfg.col_Sw(51:150,:); % Take bottom half of 'summer' colormap
    %cfg.col_Sw = winter(200); cfg.col_Sw = cfg.col_Sw(51:150,:); % Take middle half of 'winter' colormap
    cfg.col_coldestSw = 'c'; %cfg.col_Sw(1, :);
    cfg.col_midColdSw = cfg.col_Sw(25,:);
    cfg.col_middestSw = cfg.col_Sw(50,:);
    cfg.col_midWarmSw = cfg.col_Sw(75,:);
    cfg.col_warmestSw = cfg.b000ff; %cfg.col_Sw(99,:);
    cfg.Sw_alt = [  0,175,238]/255;
    
    cfg.col_MgSO4 = cool(133); cfg.col_MgSO4 = cfg.col_MgSO4(34:133,:); % Take top 3/4 of 'cool' colormap
    cfg.col_coldestMgSO4 = 'b'; %cfg.col_MgSO4(1, :);
    cfg.col_midColdMgSO4 = cfg.col_MgSO4(25,:);
    cfg.col_middestMgSO4 = cfg.col_MgSO4(50,:);
    cfg.col_midWarmMgSO4 = cfg.col_MgSO4(75,:);
    cfg.col_warmestMgSO4 = cfg.col_MgSO4(99,:);
    
    % Linestyle options
    cfg.LS_syn = '-';
    cfg.LS_orb = ':';
    cfg.LS_hrm = '-.';
    cfg.LS_Sw  = '-';
    cfg.LS_Mg = '--';
    cfg.LS_sp =  ':';
    cfg.LW_syn = 2;
    cfg.LW_orb = 2;
    cfg.LW_hrm = 2;
    cfg.LW_sal = 3;
    cfg.LW_dil = 1;
    cfg.LW_std = 2;
    cfg.LW_sound = 1.5;
    cfg.LW_seism = 1;
end
