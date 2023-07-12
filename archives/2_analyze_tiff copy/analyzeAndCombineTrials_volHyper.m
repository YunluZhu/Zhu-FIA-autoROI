
function [] = analyzeAndCombineTrials_volHyper (header,fishnum,areanum,review_fig)
	% header = "/Volumes/LabDataPro 1/2P nMLF speed/analyzed/230609 fish1/area1_pre";
    % fishnum = 1;
    % areanum = 1;

    cd(header);
	d = dir;
	isub = [d(:).isdir]; %# returns logical vector
	folders = {d(isub).name}'; % contains the name of all subfolders in that directory
	
	%% --cycle through each folder and save the mean fluorescence and dff for each
	% ROI for each trial
	numtrials=0;
	for i=1:length(folders)
        if contains(folders{i},'f')
			currentfolder = fullfile(header, folders{i});
            cd(currentfolder)
            dd = dir;

            if contains([dd(:).name],'RoiSet.zip') %if there is no ROISet.zip created in a folder, skip that folder
                % check how many stimulus are there
                numstims = 0;
                for K = 1 : length(dd)
                    file = dd(K).name;
                    if startsWith(file, 'max_Image_scan_')
                        numstims = numstims + 1;  % get the number of stimulus, each is saved as a separate tif
                    end
                end
                if numstims == 0
                    warning("No hyperstack found in: " + folders{i})
                else
                    disp(numstims + " stimulus found in: " + folders{i})
                end

                numROIs = analyzeROIs_volHyper(currentfolder,folders{i},header);  % the analyzeROIs_vol function reads the ori_.tif hyperstack file & extracts fluo of all the ROIs & calculates the dFF using the mean of two baselines
    			numtrials=numtrials+1;
            end    


        end
	end

end   
    
    
	
