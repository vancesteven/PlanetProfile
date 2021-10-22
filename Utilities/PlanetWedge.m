function PlanetWedge(input)%filename or PPstruct)


%filename
if isstr(input)
%     loadFile and get Radii from density changes
elseif isstruct(input)
    r = input.r_m*1e-3;
    phase = input.phase;
    radii = getRadii(r,phase);
else
    radii = Test;
end

colors = getColors;

if  ~isempty(radii.Ionosphere)
    nIonLayers = size(radii.Ionosphere,2);
    for iL = 1:nIonLayers
        thisIonColor = -(iL-1)*(colors.IonosphereTop - colors.IonosphereBot)/nIonLayers + colors.IonosphereTop;
        Ionosphere = wedgeR(radii.Ionosphere);Ionosphere.FaceColor = colors.Ionosphere;
    end
end
if  ~isempty(radii.R)
    IceI = wedgeR(radii.R); IceI.FaceColor=colors.IceI;
else
    error('R not found')
end
if  ~isempty(radii.Ocean)
    nOceanLayers = size(radii.Ocean,2);
    for iL=1:nOceanLayers
        thisOceanColor = -(iL-1)*(colors.OceanTop - colors.OceanBot)/nOceanLayers + colors.OceanTop;
        Ocean = wedgeR(radii.Ocean(iL));Ocean.FaceColor = thisOceanColor;
        Ocean.EdgeColor = 'none';
    end
    Ocean = wedgeR(radii.Ocean(1));Ocean.FaceColor = 'none';
end
if ~isempty(radii.IceIII)
    IceIII = wedgeR(radii.IceIII);IceIII.FaceColor = colors.IceIII;
end
if ~isempty(radii.IceV)
    IceV = wedgeR(radii.IceV);IceV.FaceColor = colors.IceV;
end
if ~isempty(radii.IceVI)
    IceVI = wedgeR(radii.IceVI);IceVI.FaceColor = colors.IceVI;
end
if ~isempty(radii.IceVII)
    IceVII = wedgeR(radii.IceVII);IceVII.FaceColor = colors.IceVII;
end
if ~isempty(radii.Rock)
    Rock = wedgeR(radii.Rock);Rock.FaceColor = colors.Rock;
end
if ~isempty(radii.Core)
    Core =  wedgeR(radii.Core);Core.FaceColor = colors.Core;
end
set(gca,'XLim',radii.R/2*[-1 1],'YLim',[0 radii.R*1.02],'XTick',[],'XColor','none','Ycolor','none')
axis tight

end % main function

% here's the test case 
% eventually match radii and densities from planetprofile outputs
% either read in the PP output file, or send the data from PlanetProfile

function radii = Test
figure(1);clf; hold on
R = 2556; % Ganymede
R = 1354; % Triton 
R = 1561; % Europa 
radii.R = R;
radii.Ocean = 1541:-1:1441;
% radii.IceIII = R-200;
% radii.IceV = R-300;
% radii.IceVI = R-600;
radii.Ionosphere = [R+100];
radii.IceIII = [];
radii.IceV = [];
radii.IceVI = [];
radii.Rock = 1440;
radii.Core = [500];
radii.IceVII = [];
end


function P = wedgeR(R)
width = pi/7;
P = wedge(pi/2-width,pi/2+width,0,0,R);
end

function colors = getColors
colors.IonosphereTop = [1 0 1]; % matlab's magenta
colors.Ionosphere = [1 0 1]; % matlab's magenta
colors.IonosphereBot = [1 0 1]; % matlab's magenta
colors.IceI = [150, 226, 241]/255;
colors.IceII = [3, 169, 252]/255;
colors.IceIII = [150, 226, 241]/255;
colors.IceV = [150, 226, 241]/255;
colors.IceVI =  [150, 226, 241]/255;
colors.IceVII =  [44, 115, 150]/255;
colors.OceanTop = [134, 149, 201]/255; % Darker and richer
%colors.OceanTop = [158, 169, 208]/255; % Lighter but more contrast
colors.OceanBot = [45, 55, 100]/255;
colors.Rock = [101, 46, 11]/255;
colors.Core = [141, 122, 121]/255;
end

function radii = getRadii(r,phase)
radii.R = r(1);
radii.Ionosphere = [];
radii.Core = r(min(find(phase>=100))); 
radii.Rock = r(min(find(phase>=50 & phase<100)));
radii.Ocean = r(find(phase==0));
% radii.Ih = r(min(find(phase==1)));
radii.IceII = r(min(find(phase==2)));
radii.IceIII = r(min(find(phase==3)));
radii.IceV = r(min(find(phase==5)));
radii.IceVI = r(min(find(phase==6)));
radii.IceVII = r(min(find(phase==7)));
end
