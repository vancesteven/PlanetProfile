function cfg = config
    % Set CALC_NEW options to 0 to re-use profile data when possible.
    % It is recommended to keep CALC_NEW=1 except when intermediate parameters
    % such as layer thicknesses will not change between runs.
    
    % Runtime options
    cfg.calc_new = 1;
    cfg.calc_new_ref = 1;
    cfg.calc_new_sound = 1;
    cfg.no_plots = 0;
    cfg.conduct = 1;
    
    % General figure options
    cfg.dft_font = 'Arial';
    cfg.dft_math = 'STIX Math';
    cfg.interpreter = 'tex';
    cfg.fig_fmt = '-depsc';
    cfg.xtn = '.eps';
    
    % Magnetic induction calculation settings
    cfg.do_eur = 1;
    cfg.do_gan = 1;
    cfg.do_cal = 1;
    cfg.do_per = 1; % Convert frequency axes to periods
    cfg.do_legend = 1;
    cfg.tables = 1; % Whether to print .tex formatting of results
    cfg.plot_fft = 1; % Whether to show plots of fourier space
    cfg.plot_contours = 1; % Contours or surfaces
    cfg.plot_V2020s = 0; % Mark the selected ocean/conductivity combos used in Vance et al. 2020
end