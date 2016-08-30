function [respMat]=MattExp1a(BlockType,SzeMaj,SzeMin,masks,feedback)

% modifications:
% 1/4/16: Implimented Training (1) and Experiment (0) options on feedback.
% 1/4/16: Implimented fixed increment for large vs small, based on small size
% 1/4/16: Fixed audio deactivate bug that would crash program when not in
% feedback mode.

if nargin < 4
   feedback=1;
end

if nargin < 3
    SzeMaj=10;
    SzeMin=30;
end

if nargin < 1
    BlockType = 3;
end

if exist('masks')==0
   load masks; 
end

% Number of trials per condition. 
if feedback==1
    TrialsPerCondition = 4; % only run 24 trials per condition demo mode
else
    TrialsPerCondition = 12; % 12 trials (i.e., cells) per stimulus condition =96 trials
end

dim = 2; % number of elements -dim:0:+dim e.g., 2 = 5*5 array
IncSze=10; % increase of size between small and large elements

% interval time in seconds 
FixTimeSecs = 1; % time fixation is presented
StimTimeSecs = .2; % time stimulus is presented
MaskTimeSecs = StimTimeSecs.*3; % time mask is presented - three times longer than SOA for stim.


%% Setup PTB with some default values
PsychDefaultSetup(2);
Screen('Preference', 'SkipSyncTests', 1);

% Seed the random number generator. Here we use the an older way to be
% compatible with older systems. Newer syntax would be rng('shuffle'). Look
% at the help function of rand "help rand" for more information
rand ('seed', round(sum(100 * clock)));

% Set the screen number to the external secondary monitor if there is one
% connected
screenNumber = max(Screen('Screens'));

% Define black, white and grey
white = WhiteIndex(screenNumber);
grey = white / 2;
black = BlackIndex(screenNumber);

% Open the screen
[window, windowRect] = PsychImaging('OpenWindow', screenNumber, black, [], 32, 2);

% Flip to clear
Screen('Flip', window);

% Query the frame duration
ifi = Screen('GetFlipInterval', window);

% Set the text size
Screen('TextSize', window,20);

% Query the maximum priority level
topPriorityLevel = MaxPriority(window);

% Get the centre coordinate of the window
[xCenter, yCenter] = RectCenter(windowRect);

% Get the size of the on screen window in pixels
% For help see: Screen WindowSize?
[screenXpixels, screenYpixels] = Screen('WindowSize', window);

% Set the blend funciton for the screen
Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

DrawFormattedText(window, 'Generating Experiment','center', 'center', white);
vbl = Screen('Flip', window);

%----------------------------------------------------------------------
%  Timing Information for the computer screen being tested
%----------------------------------------------------------------------

% interval time in frames
FixTimeFrames = round(FixTimeSecs / ifi);
StimTimeFrames = round(StimTimeSecs / ifi);
MaskTimeFrames = round(MaskTimeSecs / ifi);

% Numer of frames to wait before re-drawing
waitframes = 1;

%----------------------------------------------------------------------
%                       Keyboard information
%----------------------------------------------------------------------
% Define the keyboard keys that are listened for. We will be using the left
% and right arrow keys as response keys for the task and the escape key as
% a exit/reset key
escapeKey = KbName('ESCAPE');
leftKey = KbName('LeftArrow');
rightKey = KbName('RightArrow');

%----------------------------------------------------------------------
%                     Calculate number of trials (counterbalanced)
%----------------------------------------------------------------------

% Experiment Conditions
    %randColour=1; % 1 = blue 75%, 2 = red 75%
    %randOri=1; % 1 = vertical 75%, 2=horizontal 75%
    %randSize=2;% 1 = big 75%, 2=small 75%
switch BlockType
    case 0, % demo training block
        
    case 1, % colour only
        [randColour,randOri,randSize,randBlock] = BalanceFactors(TrialsPerCondition, 1, [1 2],[1 2],[1 2],1);
    case 2, % orientation only
        [randColour,randOri,randSize,randBlock] = BalanceFactors(TrialsPerCondition, 1, [1 2],[1 2],[1 2],2);
    case 3, % size only
        [randColour,randOri,randSize,randBlock] = BalanceFactors(TrialsPerCondition, 1, [1 2],[1 2],[1 2],3);
    case 4, % col, sze, and ori
        % note:  need to take round(TrialsPerCondition/3), else ntrials*3
        TrialsPerCondition=round(TrialsPerCondition/3);
        [randColour1,randOri1,randSize1,randBlock1] = BalanceFactors(TrialsPerCondition, 1, [1 2],[1 2],[1 2],[1]);
        [randColour2,randOri2,randSize2,randBlock2] = BalanceFactors(TrialsPerCondition, 1, [1 2],[1 2],[1 2],[2]);
        [randColour3,randOri3,randSize3,randBlock3] = BalanceFactors(TrialsPerCondition, 1, [1 2],[1 2],[1 2],[3]);
        randColour=[randColour1;randColour2;randColour3];
        randOri=[randOri1;randOri2;randOri3];
        randSize=[randSize1;randSize2;randSize3];
        randBlock=[randBlock1;randBlock2;randBlock3];
    case 5, % col and ori
        [randColour,randOri,randSize,randBlock] = BalanceFactors(TrialsPerCondition, 1, [1 2],[1 2],[1 2],[1 2]);
    case 6, % col and sze
        [randColour,randOri,randSize,randBlock] = BalanceFactors(TrialsPerCondition, 1, [1 2],[1 2],[1 2],[1 3]);
    case 7, % ori and sze
        [randColour,randOri,randSize,randBlock] = BalanceFactors(TrialsPerCondition, 1, [1 2],[1 2],[1 2],[2 3]);
    case 8, % [col ori] [col sze] & [ori sze]
        TrialsPerCondition=round(TrialsPerCondition/3);
        [randColour1,randOri1,randSize1,randBlock1] = BalanceFactors(TrialsPerCondition, 1, [1 2],[1 2],[1 2],[1 2]);
        [randColour2,randOri2,randSize2,randBlock2] = BalanceFactors(TrialsPerCondition, 1, [1 2],[1 2],[1 2],[1 3]);
        [randColour3,randOri3,randSize3,randBlock3] = BalanceFactors(TrialsPerCondition, 1, [1 2],[1 2],[1 2],[2 3]);
        randColour=[randColour1;randColour2;randColour3];
        randOri=[randOri1;randOri2;randOri3];
        randSize=[randSize1;randSize2;randSize3];
        randBlock=[randBlock1;randBlock2;randBlock3];
    case 9, % [col ori sze]
        [randColour,randOri,randSize,randBlock] = BalanceFactors(TrialsPerCondition, 1, [1 2],[1 2],[1 2],[1 2 3]);
end
    
% Duplicate the condition matrix to get the full number of trials
condMatrix = [randColour';randOri';randSize';randBlock'];

% Get the size of the matrix
[~, numTrials] = size(condMatrix);

% Randomise the conditions
shuffler = Shuffle(1:numTrials);
condMatrixShuffled = condMatrix(:, shuffler);

%----------------------------------------------------------------------
%                     Make a response matrix
%----------------------------------------------------------------------
% This is a five row matrix the first row will record the colour bias,
% the second row the orientation bias, the third row the size bias, fourth row the key
% they respond with, fifth row is trial result (1 for correct, 0 for incorrect)
% and the final row the time they took to make there response (i.e., reaction time).
respMat = nan(8, numTrials);

%----------------------------------------------------------------------
%                     Load Pregenerate Mondarian Masks 
%----------------------------------------------------------------------
%masks = mkMondrian(screenXpixels, screenYpixels,numTrials,1,3);

%----------------------------------------------------------------------
%                     Pregenerate Response Sounds (only if feedback)
%----------------------------------------------------------------------
if feedback==1
    % Initialize Sounddriver
    InitializePsychSound(1);

    % sound characteristics
    nrchannels = 2;
    freq = 48000;
    beepLengthSecs = .2;
    beepPauseTime = .2;
    repetitions=1;
    startCue = 0;
    waitForDeviceStart = 1;
    pahandle = PsychPortAudio('Open', [], 1, 1, freq, nrchannels);
    PsychPortAudio('Volume', pahandle, 1);
    correctBeep = MakeBeep(800, beepLengthSecs, freq);
    incorrectBeep = MakeBeep(200, beepLengthSecs, freq);
end


%----------------------------------------------------------------------
%                       Experimental loop
%----------------------------------------------------------------------

% Animation loop: we loop for the total number of trials
for trial = 1:numTrials

    % Cue to determine whether a response has been made
    respToBeMade = true;

    % If this is the first trial we present a start screen and wait for a key-press
    if trial == 1
        switch BlockType
            case 1, % colour only
                if feedback==1
                    DrawFormattedText(window, 'Training: Report Colour Block: Press Any Key To Begin','center', 'center', white);
                else
                    DrawFormattedText(window, 'Experiment: Report Colour Block: Press Any Key To Begin','center', 'center', white);
                end
            case 2, % orientation only
                if feedback==1
                    DrawFormattedText(window, 'Training: Report Orientation Block: Press Any Key To Begin','center', 'center', white);
                else
                    DrawFormattedText(window, 'Experiment: Report Orientation Block: Press Any Key To Begin','center', 'center', white);
                end
            case 3, % size only
                if feedback==1
                   DrawFormattedText(window, 'Training: Report Size Block: Press Any Key To Begin','center', 'center', white);
                else
                   DrawFormattedText(window, 'Experiment: Report Size Block: Press Any Key To Begin','center', 'center', white); 
                end
            case 4, % col, sze, and ori
                % note:  need to take round(TrialsPerCondition/3), else ntrials*3
                if feedback==1
                   DrawFormattedText(window, 'Training: Report either Colour, Orientation, or Size Block: Press Any Key To Begin','center', 'center', white);
                else
                   DrawFormattedText(window, 'Experiment:Report either Colour, Orientation, or Size Block: Press Any Key To Begin','center', 'center', white); 
                end
            case 5, % col and ori
                if feedback==1
                   DrawFormattedText(window, 'Training: Report Colour or Orientation Block: Press Any Key To Begin','center', 'center', white);
                else
                   DrawFormattedText(window, 'Experiment:Report Colour or Orientation Block: Press Any Key To Begin','center', 'center', white); 
                end
            case 6, % col and sze
                if feedback==1
                   DrawFormattedText(window, 'Training: Report Colour or Size Block: Press Any Key To Begin','center', 'center', white);
                else
                   DrawFormattedText(window, 'Experiment:Report Colour or Size Block: Press Any Key To Begin','center', 'center', white); 
                end
            case 7, % ori and sze
                if feedback==1
                   DrawFormattedText(window, 'Training: Report Size or Orientation Block: Press Any Key To Begin','center', 'center', white);
                else
                   DrawFormattedText(window, 'Experiment:Report Size or Orientation Block: Press Any Key To Begin','center', 'center', white); 
                end
            case 8, % [col ori] [col sze] & [ori sze]
                if feedback==1
                    DrawFormattedText(window, 'Training: Report Colour and Orientation OR Colour and Size OR Size and Orientation Block: Press Any Key To Begin','center', 'center', white);
                else
                   DrawFormattedText(window, 'Experiment:Report Colour and Orientation OR Colour and Size OR Size and Orientation Block: Press Any Key To Begin','center', 'center', white); 
                end
            case 9, % [col ori sze]
                if feedback==1
                    DrawFormattedText(window, 'Training: Report Colour or Size or Orientation Block: Press Any Key To Begin','center', 'center', white);
                else
                    DrawFormattedText(window, 'Experiment:Report Colour or Size or Orientation Block: Press Any Key To Begin','center', 'center', white);
                end
        end
        vbl = Screen('Flip', window);
        KbStrokeWait;
    end
    
    %% -------------  pre-generate up the stimulus  ------------- %
    % variables (either 1 or 2)
    randColour=condMatrixShuffled(1,trial); % 1 = blue 75%, 2 = red 75%
    randOri=condMatrixShuffled(2,trial); % 1 = vertical 75%, 2=horizontal 75%
    randSize=condMatrixShuffled(3,trial);% 1 = big 75%, 2=small 75%
    
    % Use the meshgrid command to create our base dot coordinates. This will
    % simply be a grid of equally spaced coordinates in the X and Y dimensions,
    % centered on 0,0
    % For help see: help meshgrid
    [x, y] = meshgrid(-dim:1:dim, -dim:1:dim);

    % Here we scale the grid so that it is in pixel coordinates. We just scale
    % it by the screen size so that it will fit. This is simply a
    % multiplication. 
    pixelScale = screenYpixels / (dim * 2 + 2);
    x = x .* pixelScale;
    y = y .* pixelScale;

    % Calculate the number of squares
    numSquares = numel(x);
    
    % Set the color of our squares
    SquareColors = zeros(3, numSquares);
    for a=1:numSquares
        if randColour==1
            SquareColors(:,a)=[0;0;1];
        else
            SquareColors(:,a)=[1;0;0];
        end
    end

    r=1:numSquares;
    r=Shuffle(r);
    for a=1:round((numSquares/100)*25)
        if randColour==1
            SquareColors(:,r(a))=[1;0;0];
        else
            SquareColors(:,r(a))=[0;0;1];    
        end
    end
    
    % Set the size and orientation of the Rectangle randomly between 20 and SzeMin pixels
    for a=1:numSquares
        if randSize==1
            if randOri==1
                SquareSizesX = ones(1,numSquares) .* SzeMaj+IncSze;
                SquareSizesY  = ones(1,numSquares) .* SzeMin+IncSze;
            else
                SquareSizesX = ones(1,numSquares) .* SzeMin+IncSze;
                SquareSizesY  = ones(1,numSquares) .* SzeMaj+IncSze;
            end
        else
            if randOri==1
                SquareSizesX = ones(1,numSquares) .* SzeMaj;
                SquareSizesY  = ones(1,numSquares) .* SzeMin;
            else
                SquareSizesX = ones(1,numSquares) .* SzeMin;
                SquareSizesY  = ones(1,numSquares) .* SzeMaj;
            end
        end
    end

    r=1:numSquares;
    r=Shuffle(r);
    for a=1:round((numSquares/100)*25)
        if randSize==1
            if randOri==1
                SquareSizesX(:,r(a))= SzeMin;
                SquareSizesY(:,r(a))= SzeMaj;
            else
                SquareSizesX(:,r(a))= SzeMaj;
                SquareSizesY(:,r(a))= SzeMin;
            end
        else
            if randOri==1
                SquareSizesX(:,r(a))= SzeMin+IncSze;
                SquareSizesY(:,r(a))= SzeMaj+IncSze;
            else
                SquareSizesX(:,r(a))= SzeMaj+IncSze;
                SquareSizesY(:,r(a))= SzeMin+IncSze;
            end
        end
    end

    % Make the matrix of positions for the Square. This need to be a two row
    % vector. The top row will be the X coordinate of the dot and the bottom
    % row the Y coordinate of the dot. Each column represents a single dot. For
    % help see: help reshape
    SquarePositionMatrix = [reshape(x, 1, numSquares); reshape(y, 1, numSquares)];

    % We can define a center for the dot coordinates to be relaitive to. Here
    % we set the centre to be the centre of the screen
    SquareCenter = [xCenter yCenter];

    % Make our rectangle coordinates
    allRects = nan(4, numSquares);
    for i = 1:numSquares
        baseRect=[0 0 SquareSizesY(i) SquareSizesX(i)];
        allRects(:, i) = CenterRectOnPointd(baseRect, SquarePositionMatrix(1,i)+xCenter, SquarePositionMatrix(2,i)+yCenter);
    end
    
    % make texture from mask
    %MaskTexture = Screen('MakeTexture', window,rand(768,1024,3)); % random noise mask
    MaskTexture = Screen('MakeTexture', window,masks{randi(100,1),1});

    Priority(topPriorityLevel);
    vbl = Screen('Flip', window);

    %% Flip again to sync us to the vertical retrace at the same time as
    % drawing our fixation point
    for frame = 1:FixTimeFrames
        Screen('DrawDots', window, [xCenter; yCenter], 10, white, [], 2);
        vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
    end

    % Now present the Stimulus
    for frame = 1:StimTimeFrames
        Screen('FillRect', window, SquareColors, allRects);
        vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
    end
    
    % Now present the Mask
    for frame = 1:MaskTimeFrames
        Screen('DrawTexture', window, MaskTexture, [], [], 0);
        vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi); 
    end
    
    Screen('FillRect', window, black);
    vbl = Screen('Flip', window);
    
    %%
    Test=condMatrixShuffled(4,trial);
    
    %%
    tStart = GetSecs;
    while respToBeMade == true       
        % Make a baseRect for test images
        if Test==1
            TestRectL = [0 0 SzeMaj SzeMaj];
            TestRectR = TestRectL;
            colorL=[0 0 1];
            colorR=[1 0 0];
        elseif Test==2 % orientation
            TestRectL = [0 0 SzeMin SzeMaj]; % vertical
            TestRectR = [0 0 SzeMaj SzeMin]; % horizontal
            colorL=[1 1 1];
            colorR=colorL;
        elseif Test==3
            TestRectL = [0 0 SzeMaj+IncSze SzeMin+IncSze];
            TestRectR = [0 0 SzeMaj SzeMin];
            colorL=[1 1 1];
            colorR=colorL;
        end

        % Screen X positions of our three rectangles
        TestPos = [screenXpixels * 0.45 screenXpixels * 0.55];

        % Draw the rect to the screen
        Screen('FillRect', window, colorL, CenterRectOnPointd(TestRectL, TestPos(1), yCenter));
        Screen('FillRect', window, colorR, CenterRectOnPointd(TestRectR, TestPos(2), yCenter));
        
        %% UNCOMMENT BELOW TO SEE STIM PARAMETERS ON RESPONSE SCREEN
%         Screen('DrawText', window, ['Col ',num2str(randColour)], 100, 50, [0, 0, 255, 255]);
%         Screen('DrawText', window, ['Ori ',num2str(randOri)], 100, 150, [0, 0, 255, 255]);
%         Screen('DrawText', window, ['Size ',num2str(randSize)], 100, 250, [0, 0, 255, 255]);
       
        
        
        % Check the keyboard. The person should press the
        [keyIsDown,secs, keyCode] = KbCheck;
        if keyCode(escapeKey)
            ShowCursor;
            sca;
            return
        elseif keyCode(leftKey)
            response = 1; % left
             if Test==1
                if response==randColour
                    trialResult=1;
                else
                    trialResult=0;
                end
            elseif Test==2
                if response==randOri
                    trialResult=1;
                else
                    trialResult=0;
                end  
            elseif Test==3
                if response==randSize
                    trialResult=1;
                else
                    trialResult=0;
                end
             end
             
            if feedback==1
                if trialResult==1
                    PsychPortAudio('FillBuffer', pahandle, [correctBeep; correctBeep]);
                else
                    PsychPortAudio('FillBuffer', pahandle, [incorrectBeep; incorrectBeep]);
                end
                PsychPortAudio('Start', pahandle, repetitions, startCue, waitForDeviceStart);
                WaitSecs(beepLengthSecs)
                PsychPortAudio('Stop', pahandle);
            end 
            respToBeMade = false;
        elseif keyCode(rightKey)
            response = 2; % right
            if Test==1
                if response==randColour
                    trialResult=1;
                else
                    trialResult=0;
                end
            elseif Test==2
                if response==randOri
                    trialResult=1;
                else
                    trialResult=0;
                end  
            elseif Test==3
                if response==randSize
                    trialResult=1;
                else
                    trialResult=0;
                end
            end
            if feedback==1
                if trialResult==1
                    PsychPortAudio('FillBuffer', pahandle, [correctBeep; correctBeep]);
                else
                    PsychPortAudio('FillBuffer', pahandle, [incorrectBeep; incorrectBeep]);
                end
                PsychPortAudio('Start', pahandle, repetitions, startCue, waitForDeviceStart);
                WaitSecs(beepLengthSecs)
                PsychPortAudio('Stop', pahandle);
            end 
            respToBeMade = false;
        end 

        % Flip to the screen
        vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
    end
    tEnd = GetSecs;
    rt = tEnd - tStart;

    % Record the trial data into out data matrix    
    respMat(1, trial) = BlockType;
    respMat(2, trial) = randColour;
    respMat(3, trial) = randOri;
    respMat(4, trial) = randSize;
    respMat(5, trial) = Test; % test being probed....
    respMat(6, trial) = response;
    respMat(7, trial) = trialResult;
    respMat(8, trial) = rt*1000; % in ms
end

% Close the audio device
if feedback==1
   PsychPortAudio('Close', pahandle);
end

% End of experiment screen. We clear the screen once they have made their response
if feedback==1
    DrawFormattedText(window, 'Training Block Completed \n\n Press Any Key To Exit','center', 'center', white);
else
   DrawFormattedText(window, 'Experiment Block Completed \n\n Press Any Key To Exit','center', 'center', white); 
end
Screen('Flip', window);
KbStrokeWait;

% Shut down realtime-mode:
Priority(0);


% We're done: Close all windows and textures:
Screen('CloseAll');
