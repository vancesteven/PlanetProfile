function [header,out] = read_perplex_table(mp_table)
% read_perplex_table *tab formatted table from WERAMI (perplex)
% mp_table is table name
% output are header, out is a structure with den,vp,vs, specific cp and thermal expansion

a=fopen(mp_table);
 but=fscanf(a,'%s',2);
 f=fscanf(a,'%d',1);but=fscanf(a,'%s',1);
 minT=fscanf(a,'%g',1);steT=fscanf(a,'%g',1);tstep=fscanf(a,'%g',1);
 but=fscanf(a,'%s',1);
 minP=fscanf(a,'%g',1);steP=fscanf(a,'%g',1);pstep=fscanf(a,'%g',1);
 colu=fscanf(a,'%d',1);
 type=fscanf(a,'%s',5);
fclose (a);

[den vp vs cp alpha]=textread(mp_table,'%f%f%f%f%f','headerlines',13); % takes some time if the model is large

maxt=round(minT+(tstep-1)*steT);
maxp=round(minP+(pstep-1)*steP);
xax=[minT:steT:maxt];
yax=[minP:steP:maxp];
[t p]=meshgrid(yax,xax);
%t=reshape(t,tstep,pstep);p=reshape(p,tstep,pstep);

den=reshape(den,tstep,pstep);vp=reshape(vp,tstep,pstep);vs=reshape(vs,tstep,pstep);
p=p.*1e-4; % P in GPa

header=[minT steT tstep minP steP pstep colu-2];
out.p=p;
out.t=t;
out.den=den;
out.vp=vp;
out.vs=vs;
end

