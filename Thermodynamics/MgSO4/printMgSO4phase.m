function printMgSO4phase
    nPs = 200;
    nTs = 100;
    nws = 40;
    P_MPa = linspace(0, 2200, nPs);
    T_K = linspace(246.5, 400, nTs);
    w_ppt = linspace(0, 10, nws);
    phase = zeros(nPs, nTs, nws);
    parfor i=1:nPs
        disp(['i = ' num2str(i)])
        for j=1:nTs
            for k=1:nws
                phase(i,j,k) = getIcePhaseMgSO4(P_MPa(i),T_K(j),w_ppt(k)/10);
            end
        end
    end
    save('phaseFinderMgSO4', 'P_MPa', 'T_K', 'w_ppt', 'phase');
end