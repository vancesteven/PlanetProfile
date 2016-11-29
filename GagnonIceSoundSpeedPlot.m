function GagnonIceSoundSpeedPlot

P_MPa=0.1:100:1600.1; lP = length(P_MPa);
T_K = 93:10:273; lT = length(T_K);

phase = nan(lP,lT);
for iP = 1:lP
    for iT = 1:lT
        phase(iP,iT) = getIcePhaseMgSO4(P_MPa(iP),T_K(iT),10);%p = 0 for liquid, I for I, 2 for II, 3 for III, 5 for V and 6 for VI
    end
end
[Pg,Tg]=meshgrid(P_MPa,T_K);
indsI = phase==1;
V_I = iceVelsGagnon1990(Pg(indsI),Tg(indsI),2);
indsII = phase==2;
V_II = iceVelsGagnon1990(Pg(indsII),Tg(indsII),3);
indsIII = phase==3;
V_III = iceVelsGagnon1990(Pg(indsIII),Tg(indsIII),4);
indsV = phase==5;
V_V = iceVelsGagnon1990(Pg(indsV),Tg(indsV),5);
indsVI = phase==5;
V_VI = iceVelsGagnon1990(Pg(indsVI),Tg(indsVI),6);
