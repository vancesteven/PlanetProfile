function figs = getProfileFigRefs(lbl, Tbs, thereIsCore, porous, overlay, hideplots)
    
    nTbs = length(Tbs);
    if hideplots; vis = 'off'; else; vis = 'on'; end
    
    figs.grav = findobj('type','figure','Name',lbl.gravg); % gravity
    figs.cond = findobj('type','figure','Name',lbl.panl4); % conductivity
    figs.pvt6 = findobj('type','figure','Name',[lbl.inter ' x 6']); % P-T x 6
    figs.pvt4 = findobj('type','figure','Name',[lbl.inter ' x 4']); % P-T x 4
    figs.wedg = findobj('type','figure','Name',lbl.wedge); % wedge diagrams

    if isempty(figs.grav); figs.grav = figure('visible',vis); set(gcf,'Name',lbl.gravg); end % gravity
    if isempty(figs.cond); figs.cond = figure('visible',vis); set(gcf,'Name',lbl.panl4); end % conductivity
    if isempty(figs.pvt6); figs.pvt6 = figure('visible',vis); set(gcf,'Name',[lbl.inter ' x 6']); end % P-T x 6
    if isempty(figs.pvt4); figs.pvt4 = figure('visible',vis); set(gcf,'Name',[lbl.inter ' x 4']); end % P-T x 4
    if isempty(figs.wedg); figs.wedg = figure('visible',vis); set(gcf,'Name',lbl.wedge); end % wedge diagrams

    figs.mant = findobj('type','figure','Name',lbl.mantl); % mantle density vs r
    figs.core = findobj('type','figure','Name',lbl.corsz); % core size vs mantle size
    if ~thereIsCore
        if isempty(figs.mant); figs.mant = figure('visible',vis); set(gcf,'Name',lbl.mantl); end % mantle density vs r
    else
        if isempty(figs.core); figs.core = figure('visible',vis); set(gcf,'Name',lbl.corsz); end % core size vs mantle size
    end

    if overlay
        figs.seis = findobj('type','figure','Name',lbl.seism); % seismic
        figs.gsks = findobj('type','figure','Name',lbl.gs_ks); % G_S and K_S
        if isempty(figs.seis); figs.seis = figure('visible',vis); set(gcf,'Name',lbl.seism); end % seismic
        if isempty(figs.gsks); figs.gsks = figure('visible',vis); set(gcf,'Name',lbl.gs_ks); end % G_S and K_S
    else
        [figs.seis, figs.gsks] = deal(gobjects(1,nTbs));
        for iT=1:nTbs
            seisName = [lbl.seism ' Tb = ' num2str(Tbs(iT))];
            gsksName = [lbl.gs_ks ' Tb = ' num2str(Tbs(iT))];
            szpl = findobj('type','figure','Name',seisName); % seismic
            gkpl = findobj('type','figure','Name',gsksName); % G_S and K_S
            if isempty(szpl); figs.seis(iT) = figure('visible',vis); set(gcf,'Name',seisName); else; figs.seis(iT) = szpl; end % seismic
            if isempty(gkpl); figs.gsks(iT) = figure('visible',vis); set(gcf,'Name',gsksName); else; figs.gsks(iT) = gkpl; end % G_S and K_S
        end
    end

    if porous
        if overlay
            figs.porP = findobj('type','figure','Name',lbl.porvP); % porosity vs P
            figs.porR = findobj('type','figure','Name',lbl.porvR); % porosity vs R
            figs.perm = findobj('type','figure','Name',lbl.perme); % permeability
            if isempty(figs.porP); figs.porP = figure('visible',vis); set(gcf,'Name',lbl.porvP); end % porosity vs P
            if isempty(figs.porR); figs.porR = figure('visible',vis); set(gcf,'Name',lbl.porvR); end % porosity vs R
            if isempty(figs.perm); figs.perm = figure('visible',vis); set(gcf,'Name',lbl.perme); end % permeability
        else
            [figs.porP, figs.porR, figs.perm] = deal(gobjects(1,nTbs));
            for iT=1:nTbs
                porPName = [lbl.porvP ' Tb = ' num2str(Tbs(iT))];
                porRName = [lbl.porvR ' Tb = ' num2str(Tbs(iT))];
                permName = [lbl.perme ' Tb = ' num2str(Tbs(iT))];
                pPpl = findobj('type','figure','Name',porPName); % porosity vs P
                pRpl = findobj('type','figure','Name',porRName); % porosity vs R
                pmpl = findobj('type','figure','Name',permName); % permeability
                if isempty(pPpl); figs.porP(iT) = figure('visible',vis); set(gcf,'Name',porPName); else; figs.porP(iT) = pPpl; end % porosity vs P
                if isempty(pRpl); figs.porR(iT) = figure('visible',vis); set(gcf,'Name',porRName); else; figs.porR(iT) = pRpl; end % porosity vs R
                if isempty(pmpl); figs.perm(iT) = figure('visible',vis); set(gcf,'Name',permName); else; figs.perm(iT) = pmpl; end % permeability
            end
        end
    end

end