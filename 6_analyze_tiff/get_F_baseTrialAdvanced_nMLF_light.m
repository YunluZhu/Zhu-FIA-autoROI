%% nMLF 2P anlaysis Base Trials advanced KSD
% For extracting traces of experiments with baseline condition
% 
% hardcoded for recognizing dark/light in the area folder name
% 
% returns 2 sets of traces:
%%  
% # dFF estimated by KS density
% # dFF calculated using the mean of baseline area activity
%% 
% Note that KS density estimation of baseline is applied to the condition (light 
% or dark) with baseline included, therefore it uses the ORIGINAL method (extract_dFF_ksDensity_ori) 
% 
% Note that results are quite similar between 1 & 2
% 
% 


root_dir = '/Volumes/LabDataPro/2P nMLF speed/Calcium imaging/new_2analyze';

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



        % get baseline

        A = findgroups(rawF_tabl.area);
        meanF_area = splitapply(@median, rawF,A);
        % find baseline areas
        listOfAreaNames = metadata_table.area_name(rows);
        base_rows = contains( listOfAreaNames , 'base' ) + contains( listOfAreaNames , 'Base' );
        rawF_tabl.LDcond = repelem("limbo", height(rawF_tabl))';
        dFF_res = {};
        dF_res = {};
        idff = 1;
        for condition = ["dark", "light"]
            condition_char = char(condition);
            condition_cap = replace(condition, condition_char(1), upper(condition_char(1)));
            cond_rows = contains(listOfAreaNames, condition ) + contains(listOfAreaNames, condition_cap);
            this_cond_base_rows = base_rows .* cond_rows;
            this_cond_F_rows = ~base_rows .* cond_rows;
            this_cond_F_areaNum = listOfAreaNum(logical(this_cond_F_rows));
            this_cond_rawF = rawF_tabl(rawF_tabl.area == this_cond_F_areaNum, 1:end-3);
            this_LDcond_areaNum = listOfAreaNum(logical(cond_rows));
            rows_satisfied = rawF_tabl.area == this_LDcond_areaNum(1) | rawF_tabl.area == this_LDcond_areaNum(2);
            rawF_tabl(rows_satisfied, "LDcond") = repelem({condition},height(rawF_tabl(rows_satisfied, :)))';
            
            dF = bsxfun(@minus, table2array(this_cond_rawF), meanF_area(logical(this_cond_base_rows),:));
            dFF = bsxfun(@rdivide, dF, meanF_area(logical(this_cond_base_rows),:));

            this_dF_tabl = horzcat(array2table(dF), rawF_tabl(rawF_tabl.area == this_cond_F_areaNum, end-2: end));
            this_dFF_tabl = horzcat(array2table(dFF), rawF_tabl(rawF_tabl.area == this_cond_F_areaNum, end-2: end));
            dF_res{idff} = this_dF_tabl;
            dFF_res{idff} = this_dFF_tabl;
            idff = idff+1;
        end

        dFF_baseTrial = vertcat(dFF_res{:});
        dF_baseTrial = vertcat(dF_res{:});

        % get ksDensity dFF
        A2 = findgroups(rawF_tabl.LDcond);
        [dFF_ksDensity_cell, dF_ksDensity_cell, baseline_KsDensity] = splitapply(@extract_dFF_ksDensity_advanced, rawF,A2);
        dFF_ksDensity = vertcat(dFF_ksDensity_cell{:}); % NOTE the order of conditions here is different Therefore to adjust
        dF_ksDensity = vertcat(dF_ksDensity_cell{:}); % NOTE the order of conditions here is different Therefore to adjust

        A2s = sort(A2);
        
        if mean(A2s == A2) ~= 1
            index_for_sorting = 1:length(A2);
            index_2apply_cell = splitapply(@doNothing, index_for_sorting',A2);
            index_2apply = vertcat(index_2apply_cell{:});
            dFF_ksDensity = dFF_ksDensity(index_2apply,:);

            dF_ksDensity = dF_ksDensity(index_2apply,:);
        end

        save("rawF.mat", "rawF")
        save("dFF_ksDensity.mat", "dFF_ksDensity")
        save("dF_ksDensity.mat", "dF_ksDensity")
        save("baseline_KsDensity.mat", "baseline_KsDensity")

        writetable(rawF_tabl, "rawF_df.csv")
        writetable(dFF_baseTrial, "dFF_baseTrial.csv")
        writetable(dF_baseTrial, "dF_baseTrial.csv")

        writetable(ROI_metadata, "ROI_metadata.csv")


        % % -------- may need to add ROI QC -------
    end
end