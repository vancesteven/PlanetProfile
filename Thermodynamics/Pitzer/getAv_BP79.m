function Av_BP79_mL_mol = getAv_BP79(P_bar,T_oC)

if find(P_bar >8000 | P_bar < 0 )
    error('P out of bounds [0 to 8000 bar]')
end
if find(T_oC >110 | T_oC < -22)
    error('T out of bounds [-22 to 110 oC]')
end
Pnew_bar = linspace(0,8000);
Tnew_oC = linspace(-22,100,101);

f = load('Av.mat');
Av_mL_mol = f.mat_new;

Av_BP79_mL_mol = interp2(Pnew_bar,Tnew_oC,Av_mL_mol,P_bar,T_oC);

%This is how I generated the *smooth* extrapolation of Av
%%

% % %
% Av_mat = [
%     NaN	0	100	200	400	600	800	1000
% 0	1.51	1.48	1.46	1.41	1.37	1.34	1.31
% 10	1.63	1.61	1.58	1.53	1.48	1.43	1.38
% 20	1.78	1.75	1.72	1.66	1.59	1.53	1.47
% 25	1.87	1.83	1.8	1.73	1.66	1.59	1.52
% 30	1.96	1.92	1.88	1.8	1.72	1.65	1.57
% 40	2.15	2.1	2.06	1.96	1.87	1.78	1.69
% 50	2.37	2.31	2.26	2.15	2.04	1.93	1.83
% 60	2.62	2.55	2.49	2.36	2.23	2.11	2
% 70	2.91	2.83	2.75	2.59	2.45	2.31	2.18
% 80	3.23	3.14	3.04	2.86	2.69	2.53	2.39
% 90	3.61	3.49	3.38	3.17	2.97	2.79	2.62
% 100	4.04	3.9	3.77	3.52	3.28	3.07	2.88
% 110   4.53    4.37    4.21    3.91    3.64    3.4 3.17];
% 
% Pnew_bar = linspace(0,8000);
% Tnew_oC = linspace(-22,110,101);
% mat_new = ones(101,100);
% T_mat_oC = Av_mat(2:end,1);
% P_mat_bar = Av_mat(1,2:end);
% mat_intermed = ones(101,length(P_mat_bar));
% for ij = 2:length(P_mat_bar)+1
%     mat_intermed(:,ij) = interp1(T_mat_oC,Av_mat(2:end,ij),Tnew_oC,'linear','extrap');
% end
% for jk = 1:101
%     mat_new(jk,:) = interp1(P_mat_bar,mat_intermed(jk,2:end),Pnew_bar,'linear','extrap');
% end
% % figure(123435); clf; 
% dimcell = {Tnew_oC,Pnew_bar};
% pp = spaps(dimcell,mat_new,100);
% mat_new = fnval(pp,dimcell);
% T_K = Tnew_oC + 273.15;
% P_bar = Pnew_bar;
% for ij = 1:length(Pnew_bar)
%     for jk = 1:length(Tnew_oC)
%         Av_ananth(ij,jk) = 3.73387-0.0289662*T_K(jk)+1.29461e-3*P_bar(ij)-5.62291e-6*T_K(jk).*P_bar(ij)+7.62143e-5*T_K(jk).^2+4.09944e-8*P_bar(ij).^2;
%     end
% end
% surf(Tnew_oC,Pnew_bar,pdiff(mat_new',Av_ananth'));
% % save('Av_BP79.mat','mat_new');
%%