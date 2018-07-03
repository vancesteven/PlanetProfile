function  Results=fnGval(sp,input,M)
% function to return rho,vel, G, Cp, alpha S U H K and Kp for G splines in either (P and T) or (m, P, and T)
%  ALL MKS with P in MPa
%    Usage: Results=G_sp_eval(sp,input, MW) 
%   Results=Result.rho,vel,G,Cp,alpha,S,U,H,Kt,Kp,Ks,mus,muw
%          where input is either a npts by (2 or 3) matrix of scatter points in [P  T m] or input is a cell of {P,T,m}  
%                 m in molality, P in MPa and T in K.  optional MW is molecular
%                 weight in kg/mol - needed for chemical potentials and
%                 partial molar quantities
%           rho in kg/m^3, vel in m/s, G in J/Mg Cp in J/kg/K alpha in K^(-1)
%           mu is dG/dm where m is in units determined externally
%  
% JMB 2015  
nw=1000/18.0152;  % number of moles of water in a kilogram of water

% if a 3rd input argument is given (a molecular weight) then calculate
% partial molar quantities (mu_flg=1)
mu_flg=1;
if nargin==2
    mu_flg=0;
end

%Determine whether this is P-T or P-T-m spline for gridded or scattered points 
if iscell(input) % gridded data
   nd=length(input);  
else % scattered data in columns
   [~,nd]=size(input);
end

% handle P T representations
if nd==2   % spline in P and T only
    if iscell(input)  % gridded output
        P=input{1};T=input{2};    
        [Pm,Tm]=ndgrid(P,T);
        G=fnval(sp,{P,T});
        d1T=fnval(fnder(sp,[ 0 1]),{P,T});
        d2T=fnval(fnder(sp,[ 0 2]),{P,T});
        d1P=fnval(fnder(sp,[ 1 0]),{P,T});
        dPT=fnval(fnder(sp,[ 1 1]),{P,T});
        d2P=fnval(fnder(sp,[ 2 0]),{P,T});
        d3P=fnval(fnder(sp,[ 3 0]),{P,T});
    else % scatter output
        Tm=input(:,2);
        Pm=input(:,1);
        G=fnval(sp,input')';
        d1T=fnval(fnder(sp,[ 0 1]),input')';
        d2T=fnval(fnder(sp,[0 2]),input')';
        d1P=fnval(fnder(sp,[1 0]),input')';
        dPT=fnval(fnder(sp,[1 1]),input')';
        d2P=fnval(fnder(sp,[2 0]),input')';
        d3P=fnval(fnder(sp,[3 0]),input')';
    end

else  % spline in P, T and compositions
    if iscell(input) % gridded output
        P=input{1};T=input{2};  m=input{3}; 
        m(m==0)=eps; % add eps to zero concentrations
        mflg=0;
        if mu_flg
        if (m(1)~=eps)
          m=[eps;m(:)];  % add in a zero value in order to calculate apparent quantities and remove it later
          mflg=1;
        end
        end
        [Pm,Tm,mm]=ndgrid(P,T,m);
        G=fnval(sp,{P,T,m});
        d1T=fnval(fnder(sp,[ 0 1 0]),{P,T,m});
        d2T=fnval(fnder(sp,[ 0 2 0]),{P,T,m});
        d3Tm=fnval(fnder(sp,[ 0 2 1]),{P,T,m});
        d1P=fnval(fnder(sp,[ 1 0 0]),{P,T,m});
        d2Pm=fnval(fnder(sp,[ 1 0 1]),{P,T,m});
        dPT=fnval(fnder(sp,[ 1 1 0]),{P,T,m});
        d2P=fnval(fnder(sp,[ 2 0 0]),{P,T,m});
        d3P=fnval(fnder(sp,[ 3 0 0]),{P,T,m});
        dGdm=fnval(fnder(sp,[ 0 0 1]),{P,T,m});
    else % scatter output
        Tm=input(:,2);
        Pm=input(:,1);
        mm=input(:,3); 
        m=mm;
        if mu_flg
          mm0=zeros(size(mm))+eps;
          mflg=1;
        else
            mflg=0;
        end
        G=fnval(sp,input')';
        d1T=fnval(fnder(sp,[ 0 1 0]),input')';
        d2T=fnval(fnder(sp,[0 2 0]),input')';
        d3Tm=fnval(fnder(sp,[0 2 1]),input')';
        d1P=fnval(fnder(sp,[1 0 0]),input')';
        d2Pm=fnval(fnder(sp,[1 0 1]),input')';
        dPT=fnval(fnder(sp,[1 1 0]),input')';
        d2P=fnval(fnder(sp,[2 0 0]),input')';
        d3P=fnval(fnder(sp,[3 0 0]),input')';
        dGdm=fnval(fnder(sp,[ 0 0 1]),input')';
        % calculate zero concentration derivatives to determine apparent
        % Volume and specific heat
        if mu_flg
           d2T0=fnval(fnder(sp,[0 2 0]),[Pm(:) Tm(:) mm0(:)]')';
           d1P0=fnval(fnder(sp,[1 0 0]),[Pm(:) Tm(:) mm0(:)]')';
           V0=1e-6*d1P0;  % 1e6 for MPa to Pa
           Cp0=-d2T0.*Tm(:);
        end
    end
end

Cp=-d2T.*Tm;
S=-d1T;
vel=real(sqrt(d1P.^2./(dPT.^2./d2T - d2P))); % MPa-Pa units conversion cancels
rho=1e6*d1P.^(-1);  % 1e6 for MPa to Pa
Ks=rho.*vel.^2/1e6;
alpha=1e-6*dPT.*rho; % 1e6 for MPa to Pa

U=G-1e6*Pm./rho+Tm.*S;
H=U-Tm.*S;
Kt=-d1P./d2P;
Kp=d1P.*d2P.^(-2).*d3P -1;

if iscell(input) % gridded output
  if mu_flg==1   
     f=1+M*mm;
     V=rho.^-1;
     mus=M*G + f.*dGdm;
     muw=G/nw - 1/nw*f.*mm.*dGdm;
%     muw=G/nw - (mm/nw).*dGdm;
     Vm=M*d1P +f.*d2Pm;
     Cpm=M*Cp - f.* d3Tm.*Tm;
     Cpa=(Cp.*f -repmat(squeeze(Cp(:,:,1)),1,1,length(m)))./mm;
     Va=1e6*(V.*f - repmat(squeeze(V(:,:,1)),1,1,length(m)))./mm;
     if mflg  % remove the added zero concentration ppoint
       Va=Va(:,:,2:end);
       Cpa=Cpa(:,:,2:end);
       Vm=Vm(:,:,2:end);
       Cpm=Cpm(:,:,2:end);
       G=G(:,:,2:end);
       rho=rho(:,:,2:end);
       Cp=Cp(:,:,2:end);
       S=S(:,:,2:end);
       vel=vel(:,:,2:end);
       Ks=Ks(:,:,2:end);
       alpha=alpha(:,:,2:end);
       U=U(:,:,2:end);
       H=H(:,:,2:end);
       Kt=Kt(:,:,2:end);
       Kp=Kp(:,:,2:end);
       mus=mus(:,:,2:end);
       muw=muw(:,:,2:end);
       f=f(:,:,2:end);
     end
  else
    mus=[];
    muw=[];
    Cpa=[];
    Va=[];
    f=[];
    Vm=[];
    Cpm=[];
  end
else  % Scattered data
   if mu_flg==1   
      f=1+M*mm;
      V=rho.^-1;
      mus=M*G + f.*dGdm;
      muw=G/nw - 1/nw*f.*mm.*dGdm;
      Vm=M*d1P +f.*d2Pm;
      Cpm=M*Cp - f.* d3Tm.*Tm;
      Cpa=(Cp.*f - Cp0)./mm;
      Va=1e6*(V.*f - V0)./mm;
   else
     mus=[];
     muw=[];
     Cpa=[];
     Va=[];
     f=[];
     Vm=[];
     Cpm=[];
   end
end

Results.rho=rho;
Results.Va=Va;
Results.Cpa=Cpa;
Results.Cp=Cp;
Results.G=G;
Results.vel=vel;
Results.Kt=Kt;
Results.Ks=Ks;
Results.Kp=Kp;
Results.S=S;
Results.U=U;
Results.H=H;
Results.alpha=alpha;
Results.mus=mus;
Results.muw=muw;
Results.f=f;
Results.Vm=Vm;
Results.Cpm=Cpm;

Results.d1P=d1P;
Results.d2P=d2P;
Results.d3P=d3P;
