function [header,out] = read_perplex_table_nvars(mp_table)
% read_perplex_table_nvars *tab formatted table from WERAMI (perplex)
% mp_table is table name
% output are header, out is a structure with p-t matrices of all the
% properties in the table

% use inpaintn to fill NaN and Inf values in the constructed grids
ALLOW_INPAINT = 1;

%NEED TO IMPLEMENT COUNTING OF NUMBER OF PHASES AND AUTOMATIC NAMING
a=fopen(mp_table);
disp(['Reading Perple_X Table:' mp_table])
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
 fheader=fgetl(a);
 fheader=fgetl(a);
fclose (a);

columns = strsplit(fheader);
if isempty(columns{1}) % adjust for indented format in file sent by Mohit Melwani Daswani to SDV on 20210103
    columns = columns(2:end);
end
if isempty(columns{end}) % adjust to eliminate kluge below when creating nvars SDV on 20210103
    columns = columns(1:end-1);
end
for ic = 1:length(columns)
    if strfind(columns{ic},',')
        csplit = strsplit(columns{ic},',');
        columns{ic} = csplit{1};
    elseif strfind(columns{ic},'[') % strip of characters to create variable names for volatile outputs from perplex
        csplit = strsplit(columns{ic},'[');
        csplit = strsplit(csplit{2},']');
        columns{ic} = csplit{1};
    end
    columns{ic} = replace(columns{ic},'-','_');
end
output = dlmread(mp_table,'',13,0);

if strcmp(columns{1},'P(bar)') || strcmp(columns{1},'T(K)') % for compatibility with newer releases of werami that automatically output P and T columns
        npt = 2;
        if strcmp(columns{1},'T(K)') 
            xax = output(1:tstep,1);
            yax = output(tstep:tstep:end,2);
        elseif strcmp(columns{1},'P(bar)')
            xax = output(pstep:pstep:end,2);
            yax = output(1:pstep,1);           
        end
else
    npt = 0;
    % older versions of perplex (pre v6?) didn't include columns of p and t

    % In theory, P and T are more correct if calculated from the long from
    % increments provided by perplex. However, numerical errors can also result from
    % the scheme below, as illustrated by the occasional need to use "ceil"
    % instead of "round". For this reason, I decided to just use the less
    % precise table outputs from newer versions of perplex.
    maxt=round(minT+(tstep-1)*steT);
    maxt=ceil(minT+(tstep-1)*steT);
    maxp=round(minP+(pstep-1)*steP);
    xax=[minT:steT:maxt];
    yax=[minP:steP:maxp];
end
[t,p]=meshgrid(xax,yax); % make a matrix of p and t
p=p.*1e-4; % set units of GPa

mpsplit = strsplit(mp_table,'.');
savefile = ['Thermodynamics/Perple_X/output_data/' mpsplit{1} '.mat'];
paintfile = dir(savefile);
if ~isempty(paintfile)
    load(savefile);
else
    % nvars = length(columns)-1-npt; % subtract columns of P and T, and the erroneous empty string for the endline
    nvars = length(columns)-npt; % subtract columns of P and T, and the erroneous empty string for the endline
    FLIP = 0; % many input data files have the same dimensions of P and T and not all are the same shape. This needs to be fixed. for now, FLIP is used to correct for misshapen inputs
    for iv = 1:nvars
        out.(columns{iv+npt}) = reshape(output(:,iv+npt),tstep,pstep); %read each column of data, and reshape it into a matrix of p and t
%         out.(columns{iv+npt}) = reshape(output(:,iv+npt),pstep,tstep)'; %read each column of data, and reshape it into a matrix of p and t
        if ALLOW_INPAINT
            if find(isnan(out.(columns{iv+2})))
                out.(columns{iv+2})=inpaintn(out.(columns{iv+2}));
            end
        end
        if FLIP
            X = t;        Y = 1e3*p;
        else
            X = 1e3*p;        Y = t;
        end
        header=[minT steT tstep minP steP pstep colu-2];
        try
            out.([columns{iv+npt} '_fn']) = griddedInterpolant(X,Y,out.(columns{iv+npt})); %interpolant with P in MPa
        catch
            out.([columns{iv+npt} '_fn']) = scatteredInterpolant(X(:),Y(:),out.(columns{iv+npt})(:)); %interpolant with P in MPa
        end
    end
    out.p=p;
    out.t=t;
    save(savefile);
end

