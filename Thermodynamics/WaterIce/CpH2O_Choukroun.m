function Cp_J_kg_K = CpH2O_Choukroun(P_MPa,T_K,ind)
% from Choukroun and Grasset 2010
inds = getInds(P_MPa,T_K,ind);
Cp_J_kg_K = getCp2010(T_K,inds);

function inds = getInds(P_MPa,T_K,ind)
% if USE2010
    if ~ind
    iceI_Pinds = (P_MPa<=209.5 &P_MPa>=6.11657e-4);
    iceIII_Pinds = (P_MPa<=355.0 & P_MPa>209.5);
    iceV_Pinds = (P_MPa<=618.4 & P_MPa>355.0);
    iceVI_Pinds = (P_MPa<=2050 & P_MPa>618.4);

    iceI_Tinds =(T_K<=273.15 & T_K>=0);
    iceIII_Tinds = (T_K<=258  & T_K>=0);
    iceV_Tinds = (T_K<=276 & T_K>=0);
    iceVI_Tinds = (T_K<=355 & T_K>=0);

    inds = zeros(1,length(P_MPa));
    inds = 2*(iceI_Pinds & iceI_Tinds);
    inds = inds + 3*(iceIII_Pinds & iceIII_Tinds);
    inds = inds + 4*(iceV_Pinds & iceV_Tinds);
    inds = inds + 5*(iceVI_Pinds & iceVI_Tinds);
    inds(~inds) = 1;
    else
    inds = ind*ones(1,length(P_MPa));
    end


function Cp_J_kg_K = getCp2010(T_K,inds)
%  Ice Ih, Ice II, Ice III, Ice V, Ice VI
c0 = [0 74.11 2200 820  700 940]; 
c1 = [0 7.56 0 7  7.56 5.5];

Cp_J_kg_K = polyval([c1(inds) c0(inds)],T_K);

