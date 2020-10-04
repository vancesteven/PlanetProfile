function figs = getLayeredFigRefs(lbl, Tb_K, name, PLOT_SIGS, overlay, hideplots)

    set(0,'defaultfigurecolor',[1 1 1])
    
    nTbs = length(Tb_K);
    if hideplots; vis = 'off'; else; vis = 'on'; end

    if PLOT_SIGS
        sigsName = [lbl.sigsr ' ' name];
        figs.sigs = findobj('type','figure','Name', sigsName); % sigma(r) plot
        if isempty(figs.sigs); figs.sigs = figure('visible',vis); set(gcf,'Name',sigsName); end % sigma(r) plot
    end
    
    if ~overlay
        warning('on','all')
        warning('WARNING: cfg.hold=0 is not currently implemented in LayeredInductionResponse. Only one overlayed plot will be printed for each body.');
        overlay = 1;
    end
    
    if overlay
        amphName = [lbl.amphs ' ' name];
        BxAeName = ['Bx' lbl.Aecpx ' ' name];
        ByAeName = ['By' lbl.Aecpx ' ' name];
        BzAeName = ['Bz' lbl.Aecpx ' ' name];
        figs.amph = findobj('type','figure','Name',amphName); % mag and phase
        figs.BxAe = findobj('type','figure','Name',BxAeName); % BxAe
        figs.ByAe = findobj('type','figure','Name',ByAeName); % ByAe
        figs.BzAe = findobj('type','figure','Name',BzAeName); % BzAe
        if isempty(figs.amph); figs.amph = figure('visible',vis); set(gcf,'Name',amphName); end % mag and phase
        if isempty(figs.BxAe); figs.BxAe = figure('visible',vis); set(gcf,'Name',BxAeName); end % BxAe
        if isempty(figs.ByAe); figs.ByAe = figure('visible',vis); set(gcf,'Name',ByAeName); end % ByAe
        if isempty(figs.BzAe); figs.BzAe = figure('visible',vis); set(gcf,'Name',BzAeName); end % BzAe
    else
        [figs.amph, figs.BxAe, figs.ByAe, figs.BzAe] = deal(gobjects(1,nTbs));
        for iT=1:nTbs
            amphName = [lbl.amphs ' ' name ' Tb = ' num2str(Tb_K(iT))];
            BxAeName = ['Bx' lbl.Aecpx ' ' name ' Tb = ' num2str(Tb_K(iT))];
            ByAeName = ['By' lbl.Aecpx ' ' name ' Tb = ' num2str(Tb_K(iT))];
            BzAeName = ['Bz' lbl.Aecpx ' ' name ' Tb = ' num2str(Tb_K(iT))];
            appl = findobj('type','figure','Name',amphName); % mag and phase
            bxpl = findobj('type','figure','Name',BxAeName); % BxAe
            bypl = findobj('type','figure','Name',ByAeName); % ByAe
            bzpl = findobj('type','figure','Name',BzAeName); % BzAe
            if isempty(appl); figs.amph(iT) = figure('visible',vis); set(gcf,'Name',amphName); else; figs.amph(iT) = appl; end % mag and phase
            if isempty(bxpl); figs.BxAe(iT) = figure('visible',vis); set(gcf,'Name',BxAeName); else; figs.BxAe(iT) = bxpl; end % BxAe
            if isempty(bypl); figs.ByAe(iT) = figure('visible',vis); set(gcf,'Name',ByAeName); else; figs.ByAe(iT) = bypl; end % ByAe
            if isempty(bzpl); figs.BzAe(iT) = figure('visible',vis); set(gcf,'Name',BzAeName); else; figs.BzAe(iT) = bzpl; end % BzAe
        end
    end
end