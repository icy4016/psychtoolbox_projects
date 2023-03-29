% Clear the workspace
close all;
clearvars;
sca;
Screen('Preference', 'SkipSyncTests', 1); 
% Here we call some default settings for setting up Psychtoolbox
PsychDefaultSetup(2);
ListenChar(2);
% Seed the random number generator. Here we use the an older way to be
% compatible with older systems. Newer syntax would be rng('shuffle'). Look
% at the help function of rand "help rand" for more information
rand('seed', sum(100 * clock));

% Screen Number
screenNumber = max(Screen('Screens'));

% Define black, white and grey
white = WhiteIndex(screenNumber);
grey = white / 2;
black = BlackIndex(screenNumber);

% Open the screen
[window, windowRect] = PsychImaging('OpenWindow', screenNumber, grey, [], 32, 2,...
    [], [],  kPsychNeed32BPCFloat);

% Flip to clear
Screen('Flip', window);

% Query the frame duration
ifi = Screen('GetFlipInterval', window);

% Maximum priority level
topPriorityLevel = MaxPriority(window);

% Get the centre coordinate of the window
[xCenter, yCenter] = RectCenter(windowRect);

%--------------------
% Gabor information
%--------------------

% Dimensions
gaborDimPix = 55;

% Sigma of Gaussian
sigma = gaborDimPix / 6;

% Obvious Parameters
orientation = 90;
contrast = 0.5;
aspectRatio = 1.0;

% Spatial Frequency (Cycles Per Pixel)
% One Cycle = Grey-Black-Grey-White-Grey i.e. One Black and One White Lobe
numCycles = 3;
freq = numCycles / gaborDimPix;

% Build a procedural gabor texture
gabortex = CreateProceduralGabor(window, gaborDimPix, gaborDimPix,...
    [], [0.5 0.5 0.5 0.0], 1, 0.5);

% Positions of the Gabors
dim = 3;
%  [x, y] = meshgrid(-dim:dim, [-3 3]); % 2 horizontal rows
% [x, y] = meshgrid([-3 3], -dim:dim); % 2 vertical columns
[x] = [-3 -2 -1 0 1 2 3; -3 0 0 0 0 0 3; ...
    -3 0 0 0 0 0 3; -3 0 0 0 0 0 3; -3 0 0 0 0 0 3; ...
    -3 0 0 0 0 0 3; -3 -2 -1 0 1 2 3];
[y] = [-3 -3 -3 -3 -3 -3 -3; -2 0 0 0 0 0 -2; ...
    -1 0 0 0 0 0 -1; 0 0 0 0 0 0 0; ...
    1 0 0 0 0 0 1; 2 0 0 0 0 0 2; 3 3 3 3 3 3 3];
% [x , y] = 
% Calculate the distance in "Gabor numbers" of each gabor from the center
% of the array
% dist = sqrt(x.^2 + y.^2);

% Select only the finite values
% x = x(isfinite(x));
% y = y(isfinite(y));
    
% xPos = xPos(xPos~=960);
% yPos = yPos(yPos~=600);
xPos = x .* gaborDimPix + xCenter;% + further_extend_x;
yPos = y .* gaborDimPix + yCenter;% + further_extend_y;

% further_extend_x = [-30*ones(1,7) zeros(1,5*7) 30*ones(1,7)]';

new_xPos = [xPos(1:7)-gaborDimPix xPos(43:49)+gaborDimPix]';
new_yPos = [yPos(1:7) yPos(43:49)]';

% further_extend_y = [-30 0 0 0 0 0 30 -20 0 0 0 0 0 20 -10 0 0 0 0 0 10 0 0 0 0 0 0 0 ...
%     -10 0 0 0 0 0 10 -20 0 0 0 0 0 20 -30 0 0 0 0 0 30]';

new_xPos = [new_xPos; xPos(1); xPos(8); xPos(15); xPos(22); xPos(29); xPos(36); xPos(43); ...
    xPos(1); xPos(8); xPos(15); xPos(22); xPos(29); xPos(36); xPos(43)];
new_yPos = [new_yPos; [yPos(1); yPos(8); yPos(15); yPos(22); yPos(29); yPos(36); yPos(43)]-gaborDimPix; ...    
    [yPos(7); yPos(14); yPos(21); yPos(28); yPos(35); yPos(42); yPos(49)]+gaborDimPix];
% % vertical_rows = [xPos(1:7) yPos(1:7)
% xPos = xPos(1:7) +further_extend_x;
% Count how many Gabors there are
nGabors = numel(new_xPos);

% Make the destination rectangles for all the Gabors in the array
baseRect = [0 0 gaborDimPix gaborDimPix];
allRects = nan(4, nGabors);
for i = 1:nGabors
    allRects(:, i) = CenterRectOnPointd(baseRect, new_xPos(i), new_yPos(i));
end

% Drift speed for the 2D global motion
degPerSec = 360 * 4;
degPerFrame =  degPerSec * ifi;

% Randomise the Gabor orientations and determine the drift speeds of each gabor.
% This is given by multiplying the global motion speed by the cosine
% difference between the global motion direction and the global motion.
% Here the global motion direction is 0. So it is just the cosine of the
% angle we use. We re-orientate the array when drawing
% gaborAngles = 90 .* ones(1, nGabors); %rand(1, nGabors) .* 180 - 90;
gaborAngles = [90 .* ones(1,nGabors/2) zeros(1, nGabors/2)];

change_ratio = degPerFrame * [ ...
    0.6 0.1 -0.4 -0.9 -1.4 -1.9 -2.4 ...
    0.6 0.1 -0.4 -0.9 -1.4 -1.9 -2.4 ...
    -2.4 -1.9 -1.4 -0.9 -0.4 0.1 1.6 ...
    -2.4 -1.9 -1.4 -0.9 -0.4 0.1 1.6];

degPerFrameGabors = [cosd(gaborAngles(1:nGabors/2)-90) cosd(gaborAngles(nGabors/2+1:nGabors))] ...
    .* degPerFrame + change_ratio;

% Randomise the phase of the Gabors and make a properties matrix. We could
% if we want have each Gabor with different properties in all dimensions.
% Not just orientation and drift rate as we are doing here.
% This is the power of using procedural textures
phaseLine = rand(1, nGabors) .* 360;
propertiesMat = repmat([NaN, freq, sigma, contrast, aspectRatio, 0, 0, 0],...
    nGabors, 1);
propertiesMat(:, 1) = phaseLine';

% Perform initial flip to gray background and sync us to the retrace:
vbl = Screen('Flip', window);

% Numer of frames to wait before re-drawing
waitframes = 1;

% Animation loop
while ~KbCheck

    % Set the right blend function for drawing the gabors
    Screen('BlendFunction', window, 'GL_ONE', 'GL_ZERO');

    % Batch draw all of the Gabors to screen
    Screen('DrawTextures', window, gabortex, [], allRects, gaborAngles - 90,...
        [], [], [], [], kPsychDontDoRotation, propertiesMat');

    % Change the blend function to draw an antialiased fixation point
    % in the centre of the array
    Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

    % Draw the fixation point
    Screen('DrawDots', window, [xCenter; yCenter], 10, grey, [], 2);


    % Flip our drawing to the screen
    vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);

    % Increment the phase of our Gabors
    phaseLine = phaseLine + degPerFrameGabors;
    propertiesMat(:, 1) = phaseLine';

end
ListenChar(0);
% Clean up
sca;
close all;
% clear all;