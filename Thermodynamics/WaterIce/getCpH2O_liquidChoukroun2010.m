function Cp_J_kg_K = getCpH2O_liquidChoukroun2010(T_K)
ci = [4190 9 -0.11];
Tref_K = 281.6;
Cp_J_kg_K = ones(size(T_K));
Cp_J_kg_K(T_K>=231) = ci(1) + ci(2).*exp(ci(3).*(T_K(T_K>=231)-Tref_K));
ci = [2142511.11 -35312.772 53.606 2.025691667 -0.012166 2.33191e-6 1.36462e-7 -2.68836e-10];
Cp_J_kg_K(T_K<231) = polyval(ci(end:-1:1),T_K(T_K<231));