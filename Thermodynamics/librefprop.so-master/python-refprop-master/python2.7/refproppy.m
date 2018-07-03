function rho_kgm3 = refproppy(P_MPa,T_K,varspecies,x_frac)

wd = pwd;
pydir = '/Users/svance/Dropbox/Developer/MatlabFiles/REFPROP/librefprop.so-master/python-refprop-master/python2.7';
cd(pydir);
rp = py.importlib.import_module('refprop');
udef = py.unicode('def');
species = py.list;
for is = 1:length(varspecies)
    species.append(py.unicode(varspecies{is}));
end
rp.setup(udef,species);
rp.wmol(x_frac)
rho_molar = rp.tprho(T_K,P_MPa*1000,x_frac);


