%% use for analyzing all folders in the rootdir 
% go through all folders under the root and find 
% fishnum, areanum, and amplitude and analyze all the data
% somehow only works with more than one fish folder in the dir
% or it'll loop infinitely


root_dir = '/Volumes/LabData/2P/_nMLF to ana/DONE';
ana_which = 'p';
review_fig = 0;

%%

fishFolders = {dir(root_dir).name};
fishFlags = strfind(fishFolders,'fish');  % find folders with fish
fishIdx = find(~cellfun(@isempty,fishFlags));

for i = fishIdx
    currentFish = fishFolders{i};
    fish_dir = strcat(root_dir,'/',currentFish);
    fishnum = extractAfter(fish_dir,'fish');    % get the fish number
    fishnum = str2double(fishnum);
    
    areaList = dir(fish_dir);
    isSub = [areaList.isdir];
    areaFolders = {areaList(isSub).name};
    areaFlags = strfind(areaFolders,'area');
    areaIdx = find(~cellfun(@isempty,areaFlags));
    
    for ii = areaIdx
        currentArea = areaFolders{ii};
        area_dir = strcat(fish_dir,'/',currentArea);
        areanum = extractAfter(area_dir,'area');    % get the area number
        areanum = str2double(areanum(1));  % only support one digit area number
        
%         trialList = dir(area_dir);   % each area only have one amp, extrat amp
%         isSub = [trialList.isdir];
%         trialFolders = {trialList(isSub).name};
%         trialFlags = strfind(trialFolders,'_p');
%         trialIdx = find(~cellfun(@isempty,trialFlags));
%         expName = string(trialFolders(trialIdx(1)));
%         amplitude = extractAfter(expName,'_p');    % get the area number
%         amplitude = str2double(amplitude);

        % now, run the analyze
        % header = area_dir
        % analyzeAndCombineTrials_vol(header,fishnum,areanum,ana_which,review_fig)
    end
    
end