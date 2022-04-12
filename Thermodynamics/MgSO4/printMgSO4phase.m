% For creating a lookup table for MgSO4 phase calculations under Margules
% formulation

function printMgSO4phase
    Plim = [0, 3000];
    Tlim = [246.5, 500];
    wlim = [0, 150];
    Pstep = 6.0;
    Tstep = 0.5;
    wstep = 0.5;
    P_MPa = Plim(1):Pstep:Plim(2);
    T_K =   Tlim(1):Tstep:Tlim(2);
    w_ppt = wlim(1):wstep:wlim(2);
    nPs = length(P_MPa);
    nTs = length(T_K);
    nws = length(w_ppt);
    disp(['Calculating MgSO4(aq) phase for ' num2str(nPs) ' x ' ...
           num2str(nTs) ' x ' num2str(nws) ' P x T x w grid with resolution ' ...
           num2str(Pstep) ' MPa x ' num2str(Tstep) ' K x ' num2str(wstep) ' ppt.'])
    parfor i=1:nPs
        for j=1:nTs
            for k=1:nws
                phase(i,j,k) = getIcePhaseMgSO4(P_MPa(i),T_K(j),w_ppt(k)/10);
            end
        end
        disp(['Finished i = ' num2str(i) '/' num2str(nPs)])
    end
    save('phaseLookupMgSO4', 'P_MPa', 'T_K', 'w_ppt', 'phase');
end