function checkSeaFreeze(compatVer)
    try
        seaVer = SeaFreeze_version;
    catch
        if exist('SeaFreeze.m','file')
            % We get here if SeaFreeze is on the path, but
            % SeaFreeze_version was not found. This function was added in
            % SeaFreeze version 0.9.2, the earliest PlanetProfile
            % version to adopt SeaFreeze.
            seaVer = '0.9.1';
        else
            if exist(fullfile('Thermodynamics/SeaFreeze/Matlab'),'dir') % If user has SeaFreeze but simply not added to the path
                addpath(fullfile('Thermodynamics/SeaFreeze/Matlab'));
                disp('SeaFreeze was not found on the path; SeaFreeze/Matlab was found and added to the path.')
            elseif exist(fullfile('Thermodynamics/SeaFreeze/SeaFreeze/Matlab'),'dir') % In case user has cloned the SeaFreeze git repo into the SeaFreeze folder
                addpath(fullfile('Thermodynamics/SeaFreeze/SeaFreeze/Matlab'));
                disp('SeaFreeze was not found on the path; Thermodynamics/SeaFreeze/SeaFreeze/Matlab was found and added to the path.')
            end
            try
                seaVer = SeaFreeze_version;
            catch
                error(['ERROR: SeaFreeze was not found on the path. Is it installed?' ...
                    ' Download from: https://github.com/Bjournaux/SeaFreeze and place the ' ...
                    'contents of the Matlab version into PlanetProfile/Thermodynamics/SeaFreeze/.' ...
                    ' SeaFreeze is a strict dependency; aborting.'])
            end
        end
    end
    
    % Get major and minor versions as numbers
    seaDigits = split(seaVer,'.');
    seaMaj = str2double(char(seaDigits(1)));
    seaMin = str2double(char(seaDigits(2)));
    seaSub = str2double(char(seaDigits(3)));
    compatDigits = split(compatVer,'.');
    compatMaj = str2double(char(compatDigits(1)));
    compatMin = str2double(char(compatDigits(2)));
    compatSub = str2double(char(compatDigits(3)));
    
    % Compare the version numbers
    if compatMaj == seaMaj && compatMin == seaMin && compatSub == seaSub
        return
    elseif compatMaj > seaMaj || (compatMaj == seaMaj && compatMin > seaMin) || (compatMaj == seaMaj && compatMin == seaMin && compatSub > seaSub)
    warning(['WARNING: This version of PlanetProfile is compatible with a newer version of SeaFreeze, ' ...
        'version ' compatVer ', but the installed version of SeaFreeze is ' seaVer '. Upgrading to version ' ...
        compatVer ' is strongly recommended, as missing features are likely to cause errors.'])
        % Wait a few seconds to be very sure user has a chance to see this message.
        pause(3);
    else
        warning(['WARNING: This version of PlanetProfile is compatible with SeaFreeze version ' ...
            compatVer ', but SeaFreeze version ' seaVer ' is installed. If you encounter errors, ' ...
            'try installing the compatible version of SeaFreeze.'])
    end
end