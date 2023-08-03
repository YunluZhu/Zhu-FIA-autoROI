


% root_dir = '/Volumes/LabDataPro/CB';

function metadata_table = generate_metadata(root_dir)
%%
table_fish_num = [];
table_area_num = [];
table_area_name = {};
table_rep_num = [];
table_Zframes = [];
table_X = [];
table_Y = [];
table_FR = [];
table_date = [];

fishFolders = {dir(root_dir).name};
fishFlags = strfind(fishFolders,'fish');  % find folders with fish
fishIdx = find(~cellfun(@isempty,fishFlags));
total_fish_num = length(fishIdx);
% fprintf('Found %d fish folders\n', total_fish_num);

for i = fishIdx
    currentFish = fishFolders{i};
    % disp(['In ', currentFish]);
    fish_dir = strcat(root_dir,'/',currentFish);
    this_fishNum = split(extractAfter(fish_dir,'fish')," ");    % get the fish number
    this_fishNum = str2num(this_fishNum{1});
    exp_date = split(extractBefore(fish_dir,'_fish'), "/");
    exp_date = str2num(exp_date{end});
    areaList = dir(fish_dir);
    isSub = [areaList.isdir];
    areaFolders = {areaList(isSub).name};
    areaFlags = strfind(areaFolders,'area');
    areaIdx = find(~cellfun(@isempty,areaFlags));
    total_area_num = length(areaIdx);
    % fprintf('- Found %d area folders\n', total_area_num);



    for ii = areaIdx  % loop through each area folder
        currentArea = areaFolders{ii};
        area_dir = strcat(fish_dir,'/',currentArea);
        thisArea = split(extractAfter(area_dir,'area'),"_");    
        thisAreaNum = str2double(thisArea{1}); 
        thisAreaName = thisArea{2};

        table_fish_num(end+1) = this_fishNum;
        table_area_num(end+1) = thisAreaNum;
        table_area_name{end+1} = thisAreaName;

        repList = dir(area_dir);
        isSub = [repList.isdir];
        repFolders = {repList(isSub).name};
        repFlags = strfind(repFolders,'f');
        repIdx = find(~cellfun(@isempty,repFlags));
        thisArea_rep_num = length(repIdx);
        table_rep_num(end+1) = thisArea_rep_num;

        enter_repeat = repIdx(1);
        currentRep = repFolders{enter_repeat};
        rep_dir = strcat(area_dir,'/',currentRep);

        cd(rep_dir)
        [Xpix,Ypix,FrameRate,Zframes] = findkeyparams_vol();
        table_Zframes(end+1) = Zframes;
        table_X(end+1) = Xpix;
        table_Y(end+1) = Ypix;
        table_FR(end+1) = str2double(FrameRate);
        table_date(end+1) = exp_date;
        

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

metadata_table = table(table_date', table_fish_num', table_area_num', convertCharsToStrings(table_area_name)', table_rep_num', ...
table_X', table_Y', table_Zframes', table_FR');
metadata_table.Properties.VariableNames = ["exp_date", "fish_num", "area_num", "area_name", "repeat_number", "x", "y", "z", "frame_rate"];
cd(root_dir)
date = datetime("today");
formatOut = 'yymmdd';
writetable(metadata_table, strcat(datestr(date,formatOut), " metadata.csv"));

end



