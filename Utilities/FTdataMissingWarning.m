function dlNow = FTdataMissingWarning
    disp("WARNING: Fourier Transform data was not found.")
    disp("Make sure the MagneticInduction/FieldData folder is on the path.")
    disp("If so, you likely need to download the FTdata files from the following URL:")
    disp("https://zenodo.org/record/5057572")
    disp(" ")
    response = input("Would you like to download the files now? (Y/n)",'s');
    if isempty(response)
        response = "Y";
    end
    if response ~= "n" && response ~= "N" && response ~= "no"
        zenodoBase = 'https://zenodo.org/record/5057572/files/FTdata';
        downloadBase = 'MagneticInduction/FieldData/FTdata';
        disp("Downloading FTdataAriel.mat ...")
        thisBody = 'Ariel.mat';
        websave(fullfile([downloadBase thisBody]),[zenodoBase thisBody]);
        
        disp("Downloading FTdataCallisto.mat ...")
        thisBody = 'Callisto.mat';
        websave(fullfile([downloadBase thisBody]),[zenodoBase thisBody]);
        
        disp("Downloading FTdataEuropa.mat ...")
        thisBody = 'Europa.mat';
        websave(fullfile([downloadBase thisBody]),[zenodoBase thisBody]);
        
        disp("Downloading FTdataGanymede.mat ...")
        thisBody = 'Ganymede.mat';
        websave(fullfile([downloadBase thisBody]),[zenodoBase thisBody]);
        
        disp("Downloading FTdataMiranda.mat ...")
        thisBody = 'Miranda.mat';
        websave(fullfile([downloadBase thisBody]),[zenodoBase thisBody]);
        
        dlNow = 1;
        disp("FTdata files successfully downloaded.")
    else
        dlNow = 0;
        disp("Files will not be downloaded now. Magnetic induction functions will have errors.")
    end
    
    return
end