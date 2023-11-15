% For creating a lookup table for MgSO4 phase calculations under Margules
% formulation

function printMgSO4phaseChunks
    outDir = fullfile('.');
    Plim = [0, 3000];
    Tlim = [246.5, 500];
    wlim = [0, 150];
    Pstep = 1.0;
    Tstep = 0.1;
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
    tStart = tic;
    for k=1:nws
        thisw = w_ppt(k);
        parfor i=1:nPs
            for j=1:nTs
                    phase(i,j) = getIcePhaseMgSO4(P_MPa(i),T_K(j),thisw/10);
            end
            disp(['Finished i = ' num2str(i) '/' num2str(nPs) ' after ' num2str(toc(tStart)) ' s elapsed.'])
        end
    
        fName = fullfile(outDir, ['phaseLookupMgSO4_' num2str(k)]);
        save(fName, 'P_MPa', 'T_K', 'thisw', 'phase');
        disp(['Phase lookup table for ' num2str(thisw) ' ppt printed to ' [fName '.mat'] '.' newline 'Approx. ' tRemStr(toc(tStart), k, nws) 'remaining.'])
    end
end

function tStr = tRemStr(t, i, N)
    tTot = t * N / i;
    tRem_s = tTot - t;
    tRem_d = floor(tRem_s / 86400);
    tRem_h = floor((tRem_s - tRem_d * 86400) / 3600);
    tRem_min = floor((tRem_s - tRem_d * 86400 - tRem_h * 3600) / 60);
    tRem = floor(mod(tRem_s, 60));
    if tRem_d
        tStr = [num2str(tRem_d) ' d ' num2str(tRem_h) ' h '];
    elseif tRem_h
        tStr = [num2str(tRem_h) ' h ' num2str(tRem_min) ' min ' ];
    else
        tStr = [num2str(tRem_min) ' min ' num2str(tRem) ' s ' ];
    end

end
