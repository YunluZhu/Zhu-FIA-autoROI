
function [] = analyzeAndCombineTrials_vol (header,fishnum,areanum,ana_which,review_fig)
	
	cd(header);
	d = dir;
	isub = [d(:).isdir]; %# returns logical vector
	folders = {d(isub).name}'; % contains the name of all subfolders in that directory
	
	%% --cycle through each folder and save the mean fluorescence and dff for each
	% ROI for each trial
	for i=1:length(folders)
		if contains(folders{i},'f')
			currentfolder = fullfile(header, folders{i});
            cd(currentfolder)
            dd = dir;
		    if contains([dd(:).name],'RoiSet.zip') %if there is no ROISet.zip created in a folder, skip that folder
                numROIs=analyzeROIs_vol(currentfolder,folders{i},header);  % the analyzeROIs_vol function reads the ori_.tif hyperstack file & extracts fluo of all the ROIs & calculates the dFF using the mean of two baselines
            end    
        end
	end
	
	%% --find the trials of the same stimulus amplitude from name and average together over
	 % trials of the same stimulus type. Save the data in the header path.
	
	% --pitch stimulus (w/ dynamic+static component)
	if contains(ana_which,'p')
		numtrials=0;
		for i=1:length(folders)
			if ~isempty(strfind(folders{i},sprintf('p')))
				numtrials=numtrials+1;
			end
		end
	
		dFF_allTrials_pitch = {};
		rawF_allTrials_pitch = {};
		dFF_mean_pitch = {};
		rawF_mean_pitch = {};
		dFF_sem_pitch = {};
		
		countedtrials=0;
	
		for j = 1:numtrials  
		  if j==1
			  dirname=(sprintf('f%d_a%d_p',fishnum,areanum));
          else
              dirname=(sprintf('f%d_a%d_p_00%d',fishnum,areanum,(j-1)));
		  end
		  
		  try cd(strcat(header,'/',sprintf(dirname)));
		  try 
		      load 'dFF.mat';
		      load 'rawF.mat';
		      
		      % concatenate all the trials from the same ROI
		      countedtrials=countedtrials+1;
		      for k=1:numROIs
		          if j==1
				      dFF_allTrials_pitch{k}=[dFF(:,k)];
				      rawF_allTrials_pitch{k}=[rawF(:,k)];
		          else
				      dFF_allTrials_pitch{k}(:,countedtrials)=[dFF(:,k)];
				      rawF_allTrials_pitch{k}(:,countedtrials)=[rawF(:,k)];
		          end
		      end
		  catch
		  end 
		  catch 
		  end
		end % done with concatenating rawF
		
		% quick dirty mean
		for k = 1:numROIs
	        dFF_mean_pitch{k} = mean(dFF_allTrials_pitch{k},2);  % mean across all trials
			rawF_mean_pitch{k} = mean(rawF_allTrials_pitch{k},2);  % mean across all trials
			dFF_sem_pitch{k} = sem(dFF_allTrials_pitch{k}')';

        end
        
        metadata = load('metadata.mat');  % load from the last folder, should be same across trials

		% save data
		cd(header)
		save('dFF_allTrials_pitch.mat','dFF_allTrials_pitch');
		save('rawF_allTrials_pitch.mat','rawF_allTrials_pitch');
		save('dFF_mean_pitch.mat','dFF_mean_pitch');
		save('rawF_mean_pitch.mat','rawF_mean_pitch');
		save('dFF_sem_pitch.mat','dFF_sem_pitch');
        save('metadata.mat','metadata');
        
	end  % end of pitch

    %% let's plot
	if contains(ana_which,'p')
        for k = 1:numROIs
            fprintf('ROI: %d\n',k)
            [xout_pitch,yout_pitch]=makeRibbon([1:size(dFF_mean_pitch{k},1)],dFF_mean_pitch{k}-dFF_sem_pitch{k},dFF_mean_pitch{k}+dFF_sem_pitch{k});
            
            p = figure;

            hold all; 
            plot(dFF_mean_pitch{k},'r') % blue = response to the step "pitch" stimulus
            fill(xout_pitch,yout_pitch,'b')
            alpha(0.5)
            xlim([0 size(dFF_mean_pitch{k},1)]);

            set(gca,'TickDir','out','FontSize',30,'FontName','Planscribe NF','LineWidth',5,'Box','off')
            legend({'pitch','sem'},'FontSize',10,'Box','off')
            title(sprintf('pitch dFF %d.fig',k))     

            hold off
            savefig(sprintf('pitch_dFF_%d.fig',k))     
   
            if review_fig  % do you want to review figures? if so code will pause until figure closed
                drawnow;
                try
                    waitforbuttonpress
                    % Close figure or leave it open
                    close(p)
                    disp('mouse or key pressed')
                catch
                    disp('figure closed')
                end
            else
                close(p)
            end
            
        end
    end   
    
    
	%% Remember to write control
end
	
