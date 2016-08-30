% MasterMattExp1

% Version History: Ver. 1a
% 1/4/16: Modified mask load code to check to see if masks exist (then
% load), or if not (generate). Thus, masks only generated on first run of
% experiment on a computer. 
% 1/4/16: Modified csv write code to put all results from each subject into
% its own csv file. Thus, all conditions for each subject in different
% file. 


partname=input('What are your initials?  ','s');

%% size of stimuli - from Pilot data
SzeMaj=30;
SzeMin=10;

% if program has been run before, masks exist, else generate bank of masks.
if exist('masks','file')==7
    disp('Loading Masks, please wait....');
    load masks;
    disp('Masks Loaded.');
else
    screensize = get( 0, 'Screensize' );
    disp('Generating Masks, please wait....');
    masks = mkMondrian(screensize(3),screensize(4),100,1,3);
    save masks -v7.3;
    disp('Masks Generated.');
end

% setup data file
% format for data output
formatString	= '%d, %d, %d, %d, %d, %d, %d, %4.0f, \n';
dataFileName = [partname,'_MattExp1.csv'];
if exist(dataFileName, 'file') == 0
    dataFile = fopen(dataFileName, 'a');
    fprintf(dataFile, '%s \n', 'BlockType, randColour, randOri, randSize, TestCondition, Response, trialResult, RT(ms)');
    fclose(dataFile);
end

%% setup blocks
BLOCKS=[1 2 3 5 6 7 9];
BLOCKS=Shuffle(BLOCKS);

for n=1:size(BLOCKS,2)

    % first lets begin with a training block
    [respMat_training]=MattExp1a(BLOCKS(n),SzeMaj,SzeMin,masks,1);

    % now the experiment block
    [respMat]=MattExp1a(BLOCKS(n),SzeMaj,SzeMin,masks,0);
    
    % now save data
    dataFile = fopen(dataFileName, 'a');
    for d=1:size(respMat,2)
        fprintf(dataFile, formatString, respMat(1,d),respMat(2,d),respMat(3,d),respMat(4,d),respMat(5,d),respMat(6,d),respMat(7,d),respMat(8,d));
    end
    fclose(dataFile);
end