function [rho_kgm3,v_kms,G_Jkg,Cp_JkgK,S_JkgK,alpha_Km1,Pmelt_MPa,Tmelt_K,rho_crit_kgm3,P_crit_MPa,T_crit_K] = refproppy32(PT,varspecies,x_frac,phase,flag,guess)

% PT, P in MPa and T in K can be a vector or cell, just a vector for now
% P_MPa will be converted to kPa for refprop
% x_frac is the vector of concentrations
% phase:
%  1 = liquid
%  2 = vapor
%  0 = stable phase--NOT ALLOWED (use TPFLSH)
%      (unless an initial guess is supplied for rho)
%  -1 = force the search in the liquid phase
%  -2 = force the search in the vapor phase


wd = pwd; % the working directory
pydir = '../Thermodynamics/librefprop.so-master/python-refprop-master/python3.2'; % refprop directory, relative to an input file directory.
cd(pydir);
% rp = py.importlib.import_module('refprop');
rp = py.importlib.import_module('multiRP');

species = py.list;
concentrations = py.list;
for is = 1:length(varspecies)
    species.append(varspecies{is});
    concentrations.append(x_frac(is));
end
rp.setup('def',species);

if ~exist('guess')
    guess = 0;
    flag = 0;
else
    guess = guess/w_kgmol*1e-3;
end

% incomplete
if(iscell(PT))
    rp.multirefprop(); % 'Set parent multiprocessing variables as globals'
    flg_grd=1;
    nXvar=length(PT);
    if nXvar>2
        error('First entry should be a matrix or cell of P (MPa) and T (K)');
    end
    for i=1:2
        nX(i)=length(X{i});      
    end
else
    [npts,~]=size(PT);
end

if npts==1
    try
        rhodict = rp.tprho(PT(2),PT(1)*1000,concentrations,py.int(phase),py.int(flag),guess);
        try
            tmelt = rp.meltp(PT(2),concentrations);
            Tmelt_K = tmelt{'t'};
            pmelt = rp.meltt(PT(1)*1e3,concentrations);
            Pmelt_MPa = pmelt{'p'}*1e-3;
        catch
            Tmelt_K = NaN;
            Pmelt_MPa = NaN;
        end
        
        
        rho_molL = rhodict{'D'}; % this is the input to many other functions;
        thermpy = rp.therm2(rhodict{'t'},rho_molL,rhodict{'x'});
        moldict = rp.wmol(concentrations);
        w_kgmol = 1e-3*moldict{'wmix'};

        rho_kgm3 = rho_molL.*w_kgmol*1e3; % convert D from mol/L to kg/m3
        v_kms = thermpy{'w'}*1e-3; % convert w from m/s to km/s
        S_JkgK = thermpy{'s'}./w_kgmol; % convert S from J/mol/K to J/kg/K
        Cp_JkgK = thermpy{'cp'}./w_kgmol; % convert Cp from J/mol/K to J/kg/K
        G_Jkg = thermpy{'G'}./w_kgmol; % convert G from J/mol to J/kg
        alpha_Km1 = -thermpy{'dDdt'}./rhodict{'D'};

        crit = rp.critp(concentrations);
        rho_crit_kgm3 = crit{'Dcrit'}.*w_kgmol*1e3;
        T_crit_K = crit{'tcrit'};
        P_crit_MPa = crit{'pcrit'}*1e-3;
    catch 
        [rho_kgm3,v_kms,G_Jkg,Cp_JkgK,S_JkgK,alpha_Km1,Pmelt_MPa,Tmelt_K,rho_crit_kgm3,P_crit_MPa,T_crit_K] = deal(NaN);
    end

% else
%     if iscell(PT)
%         for 
%        processlist =  
%     else
%     end
end
    
% function thermout = parseOutput(thermpy)


%return to the original working directory
cd(wd);




