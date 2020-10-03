function figs = getMagResFigRefs(lbl, bnames, hideplots)

    set(0,'defaultfigurecolor',[1 1 1])
    set(0,'defaultAxesFontSize',15)

    nbodies = length(bnames);
    if hideplots; vis = 'off'; else; vis = 'on'; end
    
    [figs.cont, figs.ffts] = deal(gobjects(1,nbodies));
    for ibody=1:nbodies
        contName = [lbl.conts ' ' char(bnames(ibody))];
        fftsName = [lbl.fftsp ' ' char(bnames(ibody))];
        cthold = findobj('type','figure','Name',contName); % induced B contour
        fthold = findobj('type','figure','Name',fftsName); % FFT plot
    
        if isempty(cthold)
            figs.cont(ibody) = figure('visible',vis);
            set(figs.cont(ibody),'Name',contName);
        else
            figs.cont(ibody) = cthold;
        end
        if isempty(fthold)
            figs.ffts(ibody) = figure('visible',vis);
            set(figs.ffts(ibody),'Name',fftsName);
        else
            figs.ffts(ibody) = fthold;
        end
    end
    
end