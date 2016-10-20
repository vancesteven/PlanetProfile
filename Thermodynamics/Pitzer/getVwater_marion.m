function Vw = getVwater_marion(P,T)
%The equilibrum extrapolation (Eq. 17 and 20) from Marion et al
%(2005).  Used to calculate the volume of 1 kg of water.
% P || bars
% T || K

Mw = 18.01528;
n_1kg = 1000/Mw;

delP = P-1.01; % this is the compression range for the model

Vl = get_Vl(T);% get molar volume of water at 1.01 bar
Kl = get_Kl(P,T);% get compressibility for equlibrum water

for ij = 1:length(delP)
    Vw(ij,:) = n_1kg*(Vl-Kl(ij,:)*delP(ij));
end
%==========================================================================
function Kl = get_Kl(P,T)
% calculate integral average compressibility, representing the compression
% of water as it is brought from 1.01 bar up to the pressure, P.
% I have integrated in P the polynomial from Marion (2005), Eq 20, and
% divided by the range of compression (P-1.01)
    Kave_fn = inline('(P-1.01).^(-1)*(8.62420e-3*P-5.06616e-5*T.*P-3.78615e-6*(P.^2)/2+2.27679e-8*P.^2*T/2+1.18481e-10*(P.^3)/3+8.20762e-8*P.*T.^2-3.63261e-11*(P.^2)/2.*T.^2-2.87099e-13*T.*(P.^3)/3)','P','T');
    for ij = 1:length(P)
        Kl(ij,:) = Kave_fn(P(ij),T);
    end
%==========================================================================
function Vl = get_Vl(T)
    TC = T-273.15;
    Vfn = inline('18.0182-1.407964e-3*TC+1.461418e-4*TC.^2','TC');
    Vl = Vfn(TC);