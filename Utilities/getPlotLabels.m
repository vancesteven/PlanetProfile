function lbl = getPlotLabels(dft_font, dft_math, interpreter)
    % Apply default settings'
    if ~exist('dft_font','var'); dft_font = 'Arial'; end
    set(0,'defaultTextFontName',dft_font)
    set(0,'defaultAxesFontName',dft_font)
    set(0,'defaultLegendFontName',dft_font)
    set(0,'defaultColorbarFontName',dft_font)
    
    if ~exist('interpreter','var'); interpreter = 'tex'; end
    set(0,'defaultTextInterpreter',interpreter);
    set(0,'defaultLegendInterpreter',interpreter);
    set(0, 'defaultAxesTickLabelInterpreter',interpreter);
    set(0, 'defaultColorbarTickLabelInterpreter',interpreter);
    
    % Create font toggles
    if ~exist('dft_lbl.math','var'); dft_math = 'STIX Math'; end
    lbl.nm   = ['\rm\fontname{' dft_font '}'];
    lbl.bnm  = ['\rm\bf\fontname{' dft_font '}'];
    lbl.math = ['\it\fontname{' dft_math '}'];
    
    % Set default font sizes
    lbl.smtext = 11;
    lbl.mdtext = 14;
    lbl.mltext = 16;
    lbl.lgtext = 18;

    % Figure and axis labels
    lbl.inter = 'P-T dependences';
    lbl.mantl = 'Mantle size vs. density';
    lbl.corsz = 'Core size vs. mantle size';
    lbl.seism = 'Seismic data and attenuation';
    lbl.panl4 = 'Conductivity with interior properties';
    lbl.wedge = 'Interior wedge diagrams';
    lbl.gravg = 'Gravity and pressure';
    lbl.porvP = 'Porosity vs. pressure';
    lbl.porvR = 'Porosity vs radius';
    lbl.perme = 'Permeability';
    lbl.gs_ks = 'G_S and K_S';
    
    lbl.conts = 'Induced B contour plot';
    lbl.xcond = [lbl.math '\sigma_{' lbl.nm 'ocean} ' lbl.nm ' (S/m)'];
    lbl.ythic = [lbl.math 'D_{' lbl.nm 'ocean} ' lbl.nm ' (km)'];
    lbl.ttamp = 'Amplitude (nT)';
    lbl.ttphs = 'Phase Delay (\circ)';
    lbl.fftsp = 'Fourier spectrum';
    
    lbl.sigsr = 'Radial conductivity function';
    lbl.sigma = [lbl.math '\sigma' lbl.nm ' (S/m)'];
    lbl.r_kms = [lbl.math 'r' lbl.nm ' (km)'];
    
    lbl.amphs = 'Amplitude and phase';
    lbl.Aecpx = 'Ae complex plane';
    lbl.T_hrs = 'Period (hr)';
    lbl.Bcomp = 'Component Amplitude (nT)';
    lbl.normA = 'Normalized Amplitude';
    
end %getPlotLabels