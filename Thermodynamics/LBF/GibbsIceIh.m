function out=GibbsIceIh(PT)
% function based on Feistal and Wagner
% out=GibbsIceIh(PT)
P=PT{1}*1e6;  % convert MPa to Pa;
T=PT{2};
nP=length(P);
nT=length(T);
pt=611.657;  %Pa
po=101325;   %Pa
Tt=273.160;
Tm=273.152519;
[pm,tm]=ndgrid(P,T);
P=pm(:);
T=tm(:);
npts=length(P);

g0k=[-632020.233449497;   %J kg21
     0.655022213658955;  % J kg21
     -1.89369929326131e-8;  % J kg21
     3.39746123271053e-15;  % J kg21
     -5.56464869058991e-22];  % J kg21
s0=189.13; % J kg21 K21
s0IAPWS=-3327.33756492168;  % J kg21 K21

r1=44.7050716285388 +i*65.6876847463481;  % J kg21 K21
tk=[3.68017112855051e-2+i*5.10878114959572e-02;
    0.337315741065416+i*0.335449415919309];

r2k=[ -72.597457432922-i*78.100842711287;  % J kg21 K21
      -5.57107698030123e-5+i*4.64578634580806e-05; % J kg21 K21
       2.34801409215913e-11-i*2.85651142904972e-11]; % J kg21 K21
rho=zeros(npts,1);
G=rho;
S=rho;
Cp=rho;
Cv=rho;
Kt=rho;
Ks=rho;
gamma=rho;
alpha=rho;



for kk=1:npts
pir=(P(kk)-po)/pt;
tau=T(kk)/Tt;

g0=sum(g0k.*pir.^(0:4)');
g0p=pt^-1*sum((1:4)'.* g0k(2:5).*pir.^(0:3)');
g0pp=pt^-2*sum((2:4)'.*(1:3)'.*g0k(3:5).*pir.^(0:2)');
r2=sum(r2k.*pir.^(0:2)');
r2p=pt^-1*sum((1:2)'.*r2k(2:3).*pir.^(0:1)');
r2pp=r2k(3)*2*pt^-2;
rk=[r1;r2];

G(kk)=g0-s0IAPWS*T(kk) + Tt*real(sum(rk.*((tk-tau).*log(tk-tau)+(tk+tau).*log(tk+tau)-2*tk.*log(tk)-tau^2.*tk.^-1)));

Gp=g0p+Tt*real(r2p*(((tk(2)-tau).*log(tk(2)-tau)+(tk(2)+tau).*log(tk(2)+tau)-2*tk(2).*log(tk(2))-tau^2.*tk(2).^-1)));
rho(kk)=1/Gp;

Gt=-s0IAPWS+real(sum(rk.*(-log(tk-tau)+log(tk+tau)-2*tau./tk)));
S(kk)=-Gt;

Gtt=Tt^-1*real(sum(rk.*((tk-tau).^-1+(tk+tau).^-1-2*tk.^-1)));
Cp(kk)=-Gtt*T(kk);

Gtp=real(r2p*(-log(tk(2)-tau)+log(tk(2)+tau)-2*tau/tk(2)));
alpha(kk)=Gtp/Gp;

Gpp=g0pp+Tt*real(r2pp*((tk(2)-tau)*log(tk(2)-tau)+(tk(2)+tau)*log(tk(2)+tau)-2*tk(2)*log(tk(2))-tau^2/tk(2)));
Kt(kk)=-(Gp/Gpp)/1e9;

Ks(kk)=(Gp*Gtt)/(Gtp^2-Gtt*Gpp)/1e9;
gamma(kk)=alpha(kk)*Ks(kk)/rho(kk)/Cp(kk)*1e9;
Cv(kk)=Cp(kk)*(1+alpha(kk)*gamma(kk)*T(kk))^-1;
end

out.rho=reshape(rho,nP,nT);
out.Kt=reshape(Kt,nP,nT);
out.Ks=reshape(Ks,nP,nT);
out.Cp=reshape(Cp,nP,nT);
out.Cv=reshape(Cv,nP,nT);
out.gamma=reshape(gamma,nP,nT);
out.alpha=reshape(alpha,nP,nT);
out.G=reshape(G,nP,nT);
out.S=reshape(S,nP,nT);



