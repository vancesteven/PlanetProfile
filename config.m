function cfg = config
    % Set CALC_NEW options to 0 to re-use profile data when possible.
    % It is recommended to keep CALC_NEW=1 except when intermediate parameters
    % such as layer thicknesses will not change between runs.
    
    cfg.calc_new = 1;
    cfg.calc_new_ref = 1;
    cfg.calc_new_sound = 1;
    cfg.no_plots = 0;
    cfg.conduct = 1;
end