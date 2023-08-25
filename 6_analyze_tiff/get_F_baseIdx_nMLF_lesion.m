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

root_dir = '/Volumes/LabDataPro/2P nMLF speed/Calcium imaging/2analyze_lesion';
root_dir = '/Volumes/LabDataPro/2P nMLF speed/Calcium imaging/2analyze_lesion_manualROI';
if_bkg_adj = false;


%%
cd(root_dir);
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
        currentFish
        this_fishNum = split(extractAfter(fish_dir,'fish')," ");    % get the fish number
        this_fishNum = str2num(this_fishNum{1});
        this_expDate = split(extractBefore(fish_dir,'_fish'), "/");
        this_expDate = str2num(this_expDate{end});
         
        rows = (metadata_table.fish_num == this_fishNum) & (metadata_table.exp_date == this_expDate);
        thisFish_rep_number = sum(metadata_table.repeat_number(rows));

        FrameRate = metadata_table.frame_rate(rows);
        FrameRate = FrameRate(1);
        Zframes = metadata_table.z(rows);
        Zframes = Zframes(1);
        ypix = metadata_table.y(rows);
        ypix = ypix(1);       
        xpix = metadata_table.x(rows);
        xpix = xpix(1);

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
        
        hw = waitbar(0, ['Reading ROIs for ' currentFish]);
        hw.Children.Title.Interpreter = 'none';
        % slow
        % for j=1:Tall  % loop through all the t
        %     for p=1:numROIs  % loop through all the rois
        %         thisROI=int16(ROI(:,:,p));
        %         sizeofROI=nnz(thisROI);
        %         rawF(j,p)= (sum(sum(img(:,:,j,zPos(p)).* uint16(thisROI)))/sizeofROI);
        %     end
        % 
        %     if mod(j,floor(Tall/20))<1e-2
        %         waitbar(j/Tall,hw);
        %     end
        % end
        for p=1:numROIs  % loop through all the rois
            thisROI=int16(ROI(:,:,p));
            sizeofROI=nnz(thisROI);
            for j=1:Tall  % loop through all the t
                rawF(j,p)= (sum(sum(img(:,:,j,zPos(p)).* uint16(thisROI)))/sizeofROI);
            end

            % if mod(j,floor(p/20))<1e-2
                waitbar(p/numROIs,hw);
            % end
        end
        close(hw)
     

        %% Create ROI metadata table
        
        ROI_metadata = [];
        fields_to_remove = ["vnRectBounds", "nVersion", "vnPosition", "mnCoordinates", "strType", "nStrokeWidth", "nStrokeColor", "nFillColor", "bSplineFit", "nPosition"];
        for sel_roi = 1:length(a)
            thisROImetadata = a{sel_roi};
            thisROImetadata.id = sel_roi;
            thisROImetadata.zPos = thisROImetadata.nPosition;
            thisROImetadata.nTop = thisROImetadata.vnRectBounds(1);
            thisROImetadata.nLeft = thisROImetadata.vnRectBounds(2);
            thisROImetadata.nBottom = thisROImetadata.vnRectBounds(3);
            thisROImetadata.nRight = thisROImetadata.vnRectBounds(4);
            thisROImetadata.xCenter = mean([thisROImetadata.nLeft, thisROImetadata.nRight]);
            thisROImetadata.yCenter = mean([thisROImetadata.nTop, thisROImetadata.nBottom]);
            
            thisROI = int16(ROI(:,:,sel_roi));
            thisROImetadata.area = nnz(thisROI);
        
            for sel_field = 1:length(fields_to_remove)
                thisROImetadata = rmfield( thisROImetadata , fields_to_remove(sel_field));
            end
            
            if isempty(ROI_metadata)
                ROI_metadata = struct2table(thisROImetadata);
            else
                ROI_metadata = [ROI_metadata; struct2table(thisROImetadata)];
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
            
            listOfAreaNum(q);
            this_area_col = repelem(listOfAreaNum(q), framesPerRepeat*listOfRepeatNum(q));
            area_col = [area_col, this_area_col];
            for repeat=1:listOfRepeatNum(q)
                repeat_col = [repeat_col, repelem(repeat, framesPerRepeat)];
            end
        end
        rawF_tabl = array2table(rawF);
        rawF_tabl.area = area_col';
        rawF_tabl.repeat = repeat_col';

        %% get ksDensity dFF
        % get groups

        G = findgroups(rawF_tabl.area);

        [dFF_ksDensity_cell, dF_ksDensity_cell, baseline_KsDensity] = splitapply(@extract_dFF_ksDensity, rawF,G);
        dFF_ksDensity = vertcat(dFF_ksDensity_cell{:});
        dF_ksDensity = vertcat(dF_ksDensity_cell{:}); 


        save("dFF_ksDensity.mat", "dFF_ksDensity")
        save("dF_ksDensity.mat", "dF_ksDensity")
        save("baseline_KsDensity.mat", "baseline_KsDensity")
        
        %% get last few second baseline
        baseline_cell = splitapply(@extract_nMLF_baseline, table2array(rawF_tabl(:,1:end-2)),G);
        baseline_idxBase_cat = vertcat(baseline_cell{:});
        baseline_idxBase = unique(baseline_idxBase_cat, 'rows');

        % get dFF
        dF_idxBase = (rawF - baseline_idxBase_cat);
        dFF_idxBase = dF_idxBase./baseline_idxBase_cat;

        save("rawF.mat", "rawF")

        save("baseline_idxBase.mat", "baseline_idxBase")
        save("dFF_idxBase.mat", "dFF_idxBase")
        save("dF_idxBase.mat", "dF_idxBase")

        writetable(rawF_tabl, "rawF_df.csv")
        writetable(ROI_metadata, "ROI_metadata.csv")

        if if_bkg_adj % obsolete
            % adjust for background fluo change in area with light
            % calculate median baseline fluo and adjust accordingly 
            % turn off for stimulus that does not change background fluo level


            baseline_tabl = array2table(baseline_idxBase_cat);
            baseline_tabl.area = rawF_tabl.area;
            B = findgroups(baseline_tabl.area);
            baseline_median = splitapply(@median, baseline_idxBase_cat, B);
            baseline_adj_forArea = baseline_median - repelem(baseline_median(1,:), length(listOfAreaNum), 1);

            baseline_adj = [];

            for q = 1:length(listOfRepeatNum)
                this_area_adj = repelem(baseline_adj_forArea(q, :), framesPerRepeat*listOfRepeatNum(q), 1);
                baseline_adj = vertcat(baseline_adj, this_area_adj);
            end    
            rawF_adj = rawF - baseline_adj;

            baseline_cell_postAdj = splitapply(@extract_baseline, rawF_adj,G);
            baseline_array_postAdj = vertcat(baseline_cell_postAdj{:});
            baseline_idxBase_adj = unique(baseline_array_postAdj, 'rows');

            dFF_idxBase_adj = (rawF_adj - baseline_array_postAdj)./baseline_array_postAdj;

            % get ksDensity dFF
            [dFF_ksDensity_adj_cell, baseline_KsDensity_adj] = splitapply(@extract_dFF_ksDensity, rawF_adj,G);
            dFF_ksDensity_adj = vertcat(dFF_ksDensity_adj_cell{:});

               % save rawF, rawF_tabl, baseline_array
            save("rawF_adj.mat", "rawF_adj")
            save("baseline_idxBase_adj.mat", "baseline_idxBase_adj") 
            save("dFF_idxBase_adj.mat", "dFF_idxBase_adj")
            save("dFF_ksDensity_adj.mat", "dFF_ksDensity_adj")
            save("baseline_KsDensity_adj.mat", "baseline_KsDensity_adj")


        end

        % % -------- may need to add ROI QC -------
    end
end