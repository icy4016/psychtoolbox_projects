% Clear the workspace
% cck 03/28/17 042317 (maing 9 vertical lines)
% 042317 creating 17 columns and turning 90 deg and 45 deg
% purpose: create back and forth square 7x4 Gabor matrices! 
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
dim = 8;
%  [x, y] = meshgrid(-dim:dim, [-3 3]); % 2 horizontal rows
% [x, y] = meshgrid([-3 3], -dim:dim); % 2 vertical columns
[x] = meshgrid(-dim:dim, [-4:1:4]);
[y] = [meshgrid([-4:1:4], -dim:dim)]';

% -4 -3 -2 -1 0 1 2 3 4; -4 -3 -2 -1 0 1 2 3 4; ...
%     -4 -3 -2 -1 0 1 2 3 4;-4 -3 -2 -1 0 1 2 3 4; -4 -3 -2 -1 0 1 2 3 4; ...
%     -4 -3 -2 -1 0 1 2 3 4;-4 -3 -2 -1 0 1 2 3 4];
% [y] = [-3 -3 -3 -3 -3 -3 -3 -3 -3;-2  -2 -2 -2 -2 -2 -2 -2 -2; ...
%        -1 -1 -1 -1 -1 -1 -1 -1 -1;0 0 0 0 0 0 0 0 0; ...
%        1 1 1 1 1 1 1 1 1; 2 2 2 2 2 2 2 2 2; 3 3 3 3 3 3 3 3 3];

% [x , y] = 
% Calculate the distance in "Gabor numbers" of each gabor from the center
% of the array
% dist = sqrt(x.^2 + y.^2);

% Select only the finite values
% x = x(isfinite(x));
% y = y(isfinite(y));

% Center the annulus coordinates in the centre of the screen
xPos = x .* gaborDimPix + xCenter;
yPos = y .* gaborDimPix + yCenter;

% Count how many Gabors there are
nGabors = numel(xPos);

% Make the destination rectangles for all the Gabors in the array
baseRect = [0 0 gaborDimPix gaborDimPix];
allRects = nan(4, nGabors);
for i = 1:nGabors
    allRects(:, i) = CenterRectOnPointd(baseRect, xPos(i), yPos(i));
end

% Drift speed for the 2D global motion
degPerSec = 360 * 4;
degPerFrame =  degPerSec * ifi;

% Randomise the Gabor orientations and determine the drift speeds of each gabor.
% This is given by multiplying the global motion speed by the cosine
% difference between the global motion direction and the global motion.
% Here the global motion direction is 0. So it is just the cosine of the
% angle we use. We re-orientate the array when drawing
% change_ratio = [2.5 0 0 0 0 0 2.5 2 0 0 0 0 0 2 ...
%     1.5 0 0 0 0 0 1.5 1 0 0 0 0 0 1 1.5 0 0 0 0 0 1.5 ...
%     2 0 0 0 0 0 2 2.5 0 0 0 0 0 2.5];
change_ratio(1,:) = degPerFrame * repmat([1.1 0.6 0.1 -0.4 -0.9 -1.4 -1.9 -2.4 -2.9],1,17);
change_ratio(2,:) = degPerFrame * repmat([-2.9 -2.4 -1.9 -1.4 -.9 -0.4 0.1 0.6 1.1],1,17);
%    [0.6 0.1 -0.4 -0.9 -1.4 -1.9 -2.4 -2.4 -1.9 -1.4 -.9 -.4 0.1 0.6 ...
%     0.6 0.1 -0.4 -0.9 -1.4 -1.9 -2.4 -2.4 -1.9 -1.4 -.9 -.4 0.1 0.6 ...
%     0.6 0.1 -0.4 -0.9 -1.4 -1.9 -2.4 -2.4 -1.9 -1.4 -.9 -.4 0.1 0.6 ...
%     0.6 0.1 -0.4 -0.9 -1.4 -1.9 -2.4 -2.4 -1.9 -1.4 -.9 -.4 0.1 0.6 ...
%     0.6 0.1 -0.4 -0.9 -1.4 -1.9 -2.4 -2.4 -1.9 -1.4 -.9 -.4 0.1 0.6 ...
%     0.6 0.1 -0.4 -0.9 -1.4 -1.9 -2.4 -2.4 -1.9 -1.4 -.9 -.4 0.1 0.6 ...
%     0.6 0.1 -0.4 -0.9 -1.4 -1.9 -2.4 -2.4 -1.9 -1.4 -.9 -.4 0.1 0.6 ...
%     0.6 0.1 -0.4 -0.9 -1.4 -1.9 -2.4 -2.4 -1.9 -1.4 -.9 -.4 0.1 0.6 ...
%     0.6 0.1 -0.4 -0.9 -1.4 -1.9 -2.4];
%     1 .4 -.2 -.8 -1.4 -2 -2.6 0 0 0 0 0 0 0 -2.6 -2 -1.4 -.8 -.2 .4 1 ...
%     -2.4 -1.9 -1.4 -.9 -.4 .1 0.6 0.6 0.1 -0.4 -0.9 -1.4 -1.9 -2.4];

% degPerFrameGabors = cosd(gaborAngles-90) .* degPerFrame + change_ratio;

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

HideCursor;

% Numer of frames to wait before re-drawing
waitframes = 1;
start_time = GetSecs;
i = 1;
alternating_count = 1;

gaborAngles = 90 .* ones(1, nGabors); %rand(1, nGabors) .* 180 - 90;
change45 = 1;

% Animation loop
% while ~KbCheck
tdown = 0;    

while (~tdown)
    % Set the right blend function for drawing the gabors
    Screen('BlendFunction', window, 'GL_ONE', 'GL_ZERO');
       
    Screen('FillRect',window, 127, [0 0 100 100]); % has to change the background because of the low contrast cond

    [S1,ERRMSG] = sprintf('%d degree', gaborAngles(1,1));
    [S2,ERRMSG] = sprintf('press ''r'' to change 45 deg in Gabor angle');
    [S3,ERRMSG] = sprintf('Esc=quit');  
    
    Screen('DrawText', window, S1, 10, 25, [255 255 255], [0 0 0]);
    Screen('DrawText', window, S2, 10, 55, [255 255 255], [0 0 0]);
    Screen('DrawText', window, S3, 10, 85, [255 255 255], [0 0 0]);

    % Batch draw all of the Gabors to screen
    Screen('DrawTextures', window, gabortex, [], allRects, gaborAngles,...
        [], [], [], [], kPsychDontDoRotation, propertiesMat');

    % Change the blend function to draw an antialiased fixation point
    % in the centre of the array
    Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

    % Draw the fixation point
%     Screen('DrawDots', window, [xCenter; yCenter], 60, grey, [], 2);
    [keyisdown,secs,keyCode] = KbCheck;
    if keyisdown %& any(keyCode(keycodes))
        if find(keyCode) == KbName('r') % increasing contrast by 5%
            gaborAngles = (90-45*(change45)) .* ones(1, nGabors);% = con_value + 5;
            change45 = change45+1;
        elseif find(keyCode) == KbName('Escape');
            tdown = 1;                       % ESCAPE key to quit               
        end            
    end
    
    while KbCheck; end %clearing key buffer
    
    % Flip our drawing to the screen
    vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);

    thetime = GetSecs;
    
    if thetime - start_time > 2 * alternating_count %use 2 sec as duration for alternation
        i = 3-i;
        alternating_count = alternating_count + 1;
    end
    if mod(nGabors,2)
        degPerFrameGabors = [cosd(gaborAngles(1:floor(nGabors/2))-90) cosd(gaborAngles(floor(nGabors/2)+1:nGabors))] ...
            .* degPerFrame + change_ratio(i,:);
    else
        degPerFrameGabors = [cosd(gaborAngles(1:nGabors/2)-90) cosd(gaborAngles(nGabors/2+1:nGabors))] ...
            .* degPerFrame + change_ratio(i,:);
    end

    phaseLine = phaseLine + degPerFrameGabors;
    propertiesMat(:, 1) = phaseLine';

end
ListenChar(0);
% Clean up
sca;
ShowCursor;
close all;
% clear all;