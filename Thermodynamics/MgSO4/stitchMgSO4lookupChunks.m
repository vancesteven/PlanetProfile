% Combines output files from printMgSO4phaseChunks into a single lookup
% table compatible with PlanetProfile.

function stitchMgSO4lookupChunks
    datDir = fullfile('Thermodynamics', 'MgSO4', 'phaseData');
    
    fNameBase = 'phaseLookupMgSO4';
    files = dir(fullfile(datDir, [fNameBase '_*']));
    nws = length(files);

    load(fullfile(datDir, files(1).name), 'P_MPa', 'T_K', 'phase', 'thisw');
    nPs = length(P_MPa);
    nTs = length(T_K);
    allPhase = zeros(nPs, nTs, nws);
    wUnsort_ppt = zeros(1,nws);
    allPhase(:,:,1) = phase;
    wUnsort_ppt(1) = thisw;
    for k=2:nws
        load(fullfile(datDir, files(k).name), 'phase', 'thisw');
        allPhase(:,:,k) = phase;
        wUnsort_ppt(k) = thisw;
    end
    [w_ppt, iSortw] = sort(wUnsort_ppt);
    phase = cast(allPhase(:,:,iSortw), 'uint8');
    fName = fullfile(datDir, fNameBase);
    save(fName, 'P_MPa', 'T_K', 'w_ppt', 'phase', '-v7.3');
    disp(['Phase lookup table printed to ' [fName '.mat'] '.'])

end