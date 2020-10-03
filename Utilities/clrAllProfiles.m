function clrAllProfiles()
    % Assumes cfg.hold = 1 and also clears all porosity-related and
    % core/mantle figures.
    
    lbl = getPlotLabels;
    
    figs.grav = findobj('type','figure','Name',lbl.gravg); % gravity
    figs.cond = findobj('type','figure','Name',lbl.panl4); % conductivity
    figs.pvt6 = findobj('type','figure','Name',[lbl.inter ' x 6']); % P-T x 6
    figs.pvt4 = findobj('type','figure','Name',[lbl.inter ' x 4']); % P-T x 4
    figs.wedg = findobj('type','figure','Name',lbl.wedge); % wedge diagrams
    figs.mant = findobj('type','figure','Name',lbl.mantl); % mantle density vs r
    figs.core = findobj('type','figure','Name',lbl.corsz); % core size vs mantle size
    figs.seis = findobj('type','figure','Name',lbl.seism); % seismic
    figs.gsks = findobj('type','figure','Name',lbl.gs_ks); % G_S and K_S
    figs.porP = findobj('type','figure','Name',lbl.porvP); % porosity vs P
    figs.porR = findobj('type','figure','Name',lbl.porvR); % porosity vs R
    figs.perm = findobj('type','figure','Name',lbl.perme); % permeability
    
    figList = [figs.grav, figs.cond, figs.pvt6, figs.pvt4, figs.wedg, figs.mant, figs.core, figs.seis, figs.gsks, figs.porP, figs.porR, figs.perm];
    for i=1:length(figList)
        if ~isempty(figList(i)); clf(figList(i)); end
    end
end