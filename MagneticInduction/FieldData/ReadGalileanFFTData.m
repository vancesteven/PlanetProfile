function FieldData = ReadGalileanFFTData
% all data are in NAIF IAU reference frames 

% 10 year duration ...
%    data points: 524288
%    Sep 1, 2018 12:00:00 AM UTC
%    Sep 1, 2028 12:00:00 AM UTC
%    Time Resolution: 0.0016611 Hz, 601.9958 seconds, 10.0333 minutes, 0.16722 hours

files = {'FTdataEuropa.mat',...
        'FTdataGanymede.mat',...
        'FTdataCallisto.mat'};
pnames = {'Europa','Ganymede','Callisto'};
FieldData = [];
for in = 1:3
    load(files{in});
    FieldData.(['BxFFT_' pnames{in}]) = FTdata.By_fft; %deliberatly switching Bx and By consistent with IAU convention used elsewhere
    FieldData.(['ByFFT_' pnames{in}]) = FTdata.Bx_fft;
    FieldData.(['BzFFT_' pnames{in}]) = FTdata.Bz_fft;
end %getFdata
FieldData.frequency = FTdata.frequency;