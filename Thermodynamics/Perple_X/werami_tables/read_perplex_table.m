function [header,out] = read_perplex_table(mp_table)
% read_perplex_table *tab formatted table from WERAMI (perplex)
% mp_table is table name
% output are header, out is a structure with den,vp,vs, specific cp and thermal expansion

a=fopen(mp_table);
 but=fscanf(a,'%s',2);
 f=fscanf(a,'%d',1);but=fscanf(a,'%s',1);
 if strcmp(but,'T(K)')
    minT=fscanf(a,'%g',1);steT=fscanf(a,'%g',1);tstep=fscanf(a,'%g',1);
 elseif strcmp(but,'P(bar)')
    minP=fscanf(a,'%g',1);steP=fscanf(a,'%g',1);pstep=fscanf(a,'%g',1);
 end
 but=fscanf(a,'%s',1);
 if strcmp(but,'T(K)')
    minT=fscanf(a,'%g',1);steT=fscanf(a,'%g',1);tstep=fscanf(a,'%g',1);
 elseif strcmp(but,'P(bar)')
    minP=fscanf(a,'%g',1);steP=fscanf(a,'%g',1);pstep=fscanf(a,'%g',1);
 end
 colu=fscanf(a,'%d',1);
 type=fscanf(a,'%s',5);
fclose (a);

if strfind(type,'P(bar)') % added by svance for compatibility with newer releases of werami that automatically output P and T columns
    if strcmp(but,'P(bar)')
        [Tread Pread den vp vs cp alpha Ks Gs]=textread(mp_table,'%f%f%f%f%f%f%f%f%f','headerlines',13); % takes some time if the model is large
        xax = Tread(1:tstep);
        yax = Pread(tstep:tstep:end);
    else
        [Pread Tread den vp vs cp alpha Ks Gs]=textread(mp_table,'%f%f%f%f%f%f%f%f%f','headerlines',13); % takes some time if the model is large
        xax = Tread(pstep:pstep:end);
        yax = Pread(1:pstep);
    end
    %      [Pread Tread den vp vs cp alpha]=textread(mp_table,'%f%f%f%f%f%f%f','headerlines',13); % takes some time if the model is large
else
    % older versions of perplex (pre v6?) didn't include columns of p and t
    [den vp vs cp alpha]=textread(mp_table,'%f%f%f%f%f','headerlines',13); % takes some time if the model is large
    % In theory, P and T are more correct if calculated from the long from
    % increments provided by perplex. However, numerical errors can result from
    % the scheme below, as illustrated by the occasional need to use "ceil"
    % instead of "round". For this reason, I decided to just use the less
    % precise table outputs from newer versions of perplex.
    maxt=round(minT+(tstep-1)*steT);
    maxt=ceil(minT+(tstep-1)*steT);
    maxp=round(minP+(pstep-1)*steP);
    xax=[minT:steT:maxt];
    yax=[minP:steP:maxp];
end
[t,p]=meshgrid(xax,yax);


% t=reshape(t,tstep,pstep);p=reshape(p,tstep,pstep);

den=reshape(den,pstep,tstep);vp=reshape(vp,pstep,tstep);vs=reshape(vs,pstep,tstep);cp=reshape(cp,pstep,tstep);alpha=reshape(alpha,pstep,tstep);

p=p.*1e-4;

header=[minT steT tstep minP steP pstep colu-2];
out.p=p;
out.t=t;
out.den=den;
out.vp=vp;
out.vs=vs;
out.cp=cp;
out.alpha=alpha;
if exist('Ks')
    Ks=reshape(Ks,pstep,tstep);Gs=reshape(Gs,pstep,tstep);
    ks=Ks*1e-4;Gs=Gs*1e-4; % values in in GPa
    out.ks = ks;
    out.gs = Gs;
end

end

