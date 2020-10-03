function clrAllLayered(name)
    % Assumes cfg.hold = 1.
    
    lbl = getPlotLabels;
    
    figs.sigs = findobj('type','figure','Name',[lbl.sigsr ' ' name]); % sigma(r) plot
    figs.amph = findobj('type','figure','Name',[lbl.amphs ' ' name]); % mag and phase
    figs.BxAe = findobj('type','figure','Name',['Bx' lbl.Aecpx ' ' name]); % BxAe
    figs.ByAe = findobj('type','figure','Name',['By' lbl.Aecpx ' ' name]); % ByAe
    figs.BzAe = findobj('type','figure','Name',['Bz' lbl.Aecpx ' ' name]); % BzAe
    
    figList = [figs.sigs, figs.amph, figs.BxAe, figs.ByAe, figs.BzAe];
    for i=1:length(figList)
        if ~isempty(figList(i)); clf(figList(i)); end
    end
end