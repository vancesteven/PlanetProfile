function cfg = config
    % Set CALC_NEW options to 0 to re-use profile data when possible.
    % It is recommended to keep CALC_NEW=1 except when intermediate parameters
    % such as layer thicknesses will not change between runs.
    
    % Runtime options
    cfg.calc_new =       1;
    cfg.calc_new_ref =   1;
    cfg.calc_new_sound = 1;
    cfg.calc_new_induc = 1;
    cfg.skip_profiles =  0; % Whether to skip past all PlanetProfile.m plotting and calculations
    cfg.no_plots =       0;
    cfg.hold =           1; % Whether to overlay runs when possible
    cfg.conduct =        1; % Calculate electrical conductivity
    cfg.disp_layers =    1; % Whether to display layer depths and heat fluxes for user
    cfg.disp_tables =    1; % Whether to print latex-formatted tables to Matlab command line
    cfg.deprecated =     0; % Whether to allow deprecated code to run. Will often cause errors.
    
    % Magnetic induction calculation settings
    cfg.do_eur = 1;
    cfg.do_gan = 1;
    cfg.do_cal = 1;
    cfg.do_per = 1; % Convert frequency axes to periods
    cfg.do_legend = 1;
    cfg.plot_fft = 1; % Whether to show plots of fourier space
    cfg.plot_contours = 1; % Contours or surfaces
    cfg.plot_V2020s = 1; % Mark the selected ocean/conductivity combos used in Vance et al. 2020
    cfg.intMethod = 'makina'; % Interpolation method. Certain ones can cause wiggles, notably 'linear'.
    cfg.npts_k = 50;
    cfg.npts_D = 60;
    cfg.np_intp = 200;
    cfg.npts_w = 100;
    cfg.np_wfine = 1000;
    cfg.opts_odeParams = odeset('RelTol',1e-10,'AbsTol',1e-10,'MaxStep', 2e3,'InitialStep',1e-2);
    cfg.opts_odeLayers = odeset('RelTol',1e-8, 'AbsTol',1e-10,'MaxStep',10e3,'InitialStep',1e-2);
    
    % General figure options
    cfg.dft_font = 'Arial';
    cfg.dft_math = 'STIX Math';
    cfg.interpreter = 'tex';
    cfg.fig_fmt = '-depsc';
    cfg.xtn = '.eps';
    
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
    
    cfg.col_MgSO4 = cool(133); cfg.col_MgSO4 = cfg.col_MgSO4(34:133,:); % Take top 3/4 of 'cool' colormap
    cfg.col_coldestMgSO4 = 'b'; %cfg.col_MgSO4(1, :);
    cfg.col_midColdMgSO4 = cfg.col_MgSO4(25,:);
    cfg.col_middestMgSO4 = cfg.col_MgSO4(50,:);
    cfg.col_midWarmMgSO4 = cfg.col_MgSO4(75,:);
    cfg.col_warmestMgSO4 = cfg.col_MgSO4(99,:);
    
    % Linestyle options
    cfg.ls_syn = '-';
    cfg.ls_orb = ':';
    cfg.ls_hrm = '-.';
    cfg.ls_Sw  = '-';
    cfg.ls_Mg = '--';
    cfg.ls_sp =  ':';
    cfg.lw_syn = 2;
    cfg.lw_orb = 2;
    cfg.lw_hrm = 2;
    cfg.lw_sal = 3;
    cfg.lw_dil = 1;
    cfg.lw_std = 2;
    cfg.lw_sound = 1.5;
    cfg.lw_seism = 1;
end