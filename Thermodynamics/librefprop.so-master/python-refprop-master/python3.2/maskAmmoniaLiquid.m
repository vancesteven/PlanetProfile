function mask = maskAmmoniaLiquid(PT)
% P in MPa, T in K
if iscell(PT)
    [png,tng]=meshgrid(PT{1},PT{2});
elseif ndims(PT>2)
    png = PT(:,1);
    tng = PT(:,2);
end
Pvg_MPa= vaporAmmonia(tng);
mask = png>=Pvg_MPa;
Psg_MPa = solidAmmonia(tng);
mask = double(mask & png <= Psg_MPa)';
mask(mask==0)=nan;