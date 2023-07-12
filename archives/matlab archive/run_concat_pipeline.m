%% nMLF 2P anlaysis 
% modified to accomodate nMLF pipeline where ROIs are generated using particle 
% analysis based on motion-corrected concatenated dataset
% 
% what this script does is to:
% 
% 1. get experiment parameters, including time frames and slices per repeat
% 
% per area/condition
% 
% 2. get area/condition number, and repeat number per area by looking for
% 
% area and rep folders
% 
% 3. slice master hyperstack into individual repeats
% 
% 4. calculate dFF for individual repeats using theoriotical baseline
% 
% 5. run quality control and average repeats per area
% 
% 

root_dir = '/Volumes/LabDataPro/2P nMLF speed/Calcium imaging/light_2analyze';
cd(root_dir);
%%
metadata_table = generate_metadata(root_dir);
%%
sti_per_rep = 4;

fishFolders = {dir(root_dir).name};
fishFlags = strfind(fishFolders,'fish');  % find folders with fish
fishIdx = find(~cellfun(@isempty,fishFlags));
total_fish_num = length(fishIdx);
% fprintf('Found %d fish folders\n', total_fish_num);

for i = fishIdx
    currentFish = fishFolders{i};
    % disp(['In ', currentFish]);
    fish_dir = strcat(root_dir,'/',currentFish);
    cd(fish_dir);
    dd = dir;

    if contains([dd(:).name],'RoiSet.zip')

        this_fishNum = split(extractAfter(fish_dir,'fish')," ");    % get the fish number
        this_fishNum = str2double(this_fishNum{1});
        
        rows = metadata_table.fish_num == this_fishNum;
        thisFish_rep_number = sum(metadata_table.repeat_number(rows));

        FrameRate = metadata_table.frame_rate(rows);
        FrameRate = FrameRate(1)
        Zframes = metadata_table.z(rows);
        Zframes = Zframes(1)
        ypix = metadata_table.y(rows);
        ypix = ypix(1)       
        xpix = metadata_table.x(rows);
        xpix = xpix(1)

        Tall_hyp1 = floor(FrameRate/(Zframes+1)*20)*sti_per_rep*thisFish_rep_number;  
        Tall_hyp2 = round(FrameRate/(Zframes+1)*20)*sti_per_rep*thisFish_rep_number;
        
        
        ROI_file='RoiSet.zip'; %path to ROI file (.zip)
    
        % % Need ReadImageJROI.m script
        [a,ROI,zPos] = ReadImageJROI_vol(ROI_file,[ypix,xpix]);%need to reformat a using poly2mask
        [a,ROI] = fixROIs(a,ROI); % ROI are the binary masks for which pixels to analyze vs ignore; is a d1xd2xnumberofROIs matrix
        % 
        numROIs=size(ROI,3);
        
        file_names = {dd.name};
        imgFlag = startsWith(file_names,'masterHyperstack');  
        imgName = file_names{imgFlag};
        % read image
        img = TIFFStack(imgName);  % TIFFStack() is used here to read hyperstack
        [~,~,Tall,Zall] = size(img);  
        % % sanity check 
        assert(ismember(Tall,[Tall_hyp1,Tall_hyp2]),sprintf('WARNING: time slices does not match frame rate, check this folder: %s', currentFish))
        % 
        rawF = zeros(Tall,numROIs);
        %% open stack, save the deltaF/F profile of the entire trial
        
        
        % % slow
        for j=1:Tall  % loop through all the t
     
            for p=1:numROIs  % loop through all the rois
                thisROI=int16(ROI(:,:,p));
                sizeofROI=nnz(thisROI);
                rawF(j,p)= (sum(sum(img(:,:,j,zPos(p)).* uint16(thisROI)))/sizeofROI);
            end
        end
      
        % % use area number and repeat number
        % % for the concatenated raw matrix where each row is a frame and each
        % column is an ROI, generate area column and repeat number column
    
        listOfRepeatNum = metadata_table.repeat_number(rows);
        listOfAreaNum = metadata_table.area_num(rows);
    
        framesPerRepeat = Tall / thisFish_rep_number;
        repeat_col = [];
        area_col = [];
        for q = 1:length(listOfRepeatNum)
            
            listOfAreaNum(q)
            this_area_col = repelem(listOfAreaNum(q), framesPerRepeat*listOfRepeatNum(q));
            area_col = [area_col, this_area_col];
            for repeat=1:listOfRepeatNum(q)
                repeat_col = [repeat_col, repelem(repeat, framesPerRepeat)];
            end
        end
        rawF_tabl = array2table(rawF);
        rawF_tabl.area = area_col';
        rawF_tabl.repeat = repeat_col';

        
        % get baseline
        
        G = findgroups(rawF_tabl.area,rawF_tabl.repeat);
        baseline_cell = splitapply(@extract_baseline, table2array(rawF_tabl(:,1:end-2)),G);
        baseline_array = vertcat(baseline_cell{:});
        
        
        % get dFF
        dFF = (rawF - baseline_array)./baseline_array;
        
        % save rawF, rawF_tabl, baseline_array
        save("rawF.mat", "rawF")
        save("baseline.mat", "baseline_array")
        save("dFF.mat", "dFF")
        writetable(rawF_tabl, "rawF_df.csv")

        % % -------- may need to add ROI QC -------
    end
end


%%



%%
function baseline_heightMatch = extract_baseline(response)
    sti_repeat_number = 4;
    vol_rate = 1.2879;%FrameRate/(Zall+1);
    baseline_sec = 3;
    sti_dur = height(response)/sti_repeat_number;
    baseline_dur_idx = round(vol_rate*baseline_sec);
    basedex1 = (sti_dur-baseline_dur_idx):sti_dur;
    basedex2 = (sti_dur*2-baseline_dur_idx):sti_dur*2;
    basedex3 = (sti_dur*3-baseline_dur_idx):sti_dur*3;
    basedex4 = (sti_dur*4-baseline_dur_idx):sti_dur*4;
    baseline = response(cat(2,basedex1,basedex2,basedex3,basedex4),:);
    meanBaseline = mean(baseline,1);
    baseline_heightMatch = {repelem(meanBaseline, height(response),1)};
end