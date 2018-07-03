function EOS = readh5spline(filename)
% make an EOS spline from hdf5 data

coefs = h5read(filename,'/coefs');
n = size(coefs);

EOS.sp.form = 'B-';
EOS.sp.dim = 1;
EOS.sp.coefs = coefs;
EOS.sp.number = n;

tx = h5read(filename,'/tx');
switch length(n)
    case 1
    EOS.sp.knots = {tx'};
    EOS.sp.order = [length(tx)]-n; 
    case 2
    ty = h5read(filename,'/ty');
    EOS.sp.knots = {tx',ty'};
    EOS.sp.order = [length(tx) length(ty)]-n;
    case 3
    ty = h5read(filename,'/ty');
    tz = h5read(filename,'/tz');
    EOS.sp.knots = {tx',ty',tz'};
    EOS.sp.order = [length(tx) length(ty) length(tz)]-n;
end

% additional lines here to read in the metadata and molecular weight
% contained in the hdf file....
EOS.MW = [0.058443];
