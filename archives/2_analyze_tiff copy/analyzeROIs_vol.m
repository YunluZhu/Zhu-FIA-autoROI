%% this function is called by analyzeAndCombineAllTrials, so give the experiment header to that function and then this one will work for each individual trial

function [numROIs] = analyzeROIs_vol(folderpath,foldername,parentDir)

    %% Extract key para
    [Xpix,Ypix,FrameRate,Zframes] = findkeyparams_vol();
    FrameRate = str2double(FrameRate);
    % calculate hypothetical T all
    Tall_hyp1 = floor(FrameRate/(Zframes+1)*15)*4;  % assuming using 1 flyback and 95 stimulus protocol
    Tall_hyp2 = round(FrameRate/(Zframes+1)*15)*4;
    
    %% Import ROI
    height = Ypix;
    width = Xpix;
    
    ROI_file='RoiSet.zip'; %path to ROI file (.zip)

    % Need ReadImageJROI.m script
    [a,ROI,zPos] = ReadImageJROI_vol(ROI_file,[height,width]);%need to reformat a using poly2mask
    [a,ROI] = fixROIs(a,ROI); % ROI are the binary masks for which pixels to analyze vs ignore; is a d1xd2xnumberofROIs matrix

    numROIs=size(ROI,3);
    cd(parentDir)  % go back to the parent folder (area folder) to find the ori Tiff file

    tmp = TIFFStack(['ori_' foldername '.tif']);  % TIFFStack() is used here to read hyperstack
    [~,~,Tall,Zall] = size(tmp);  
    
    % sanity check 
    assert(ismember(Tall,[Tall_hyp1,Tall_hyp2]),sprintf('WARNING: time slices does not match frame rate, check this folder: %s',fullfile(parentDir,foldername)))
    
    rawF = zeros(Tall,numROIs);
    
    %% open stack, save the deltaF/F profile of the entire trial
    for i=1:Tall  % loop through all the t
        for j=1:numROIs  % loop through all the rois
            thisROI=int16(ROI(:,:,j));
            sizeofROI=nnz(thisROI);
            rawF(i,j)= (sum(sum(tmp(:,:,i,zPos(j)).* uint16(thisROI)))/sizeofROI);
        end
    end

    %% Extract baseline and responses
    dur_dex = Tall/4;

    
    
    
    
    basedex1 = 1:dur_dex;
    basedex2 = (3*dur_dex+1):Tall;
    downdex = (dur_dex+1):(2*dur_dex);
    updex = (2*dur_dex+1):(3*dur_dex);
    
    baseline1 = rawF(basedex1,:);
    baseline2 = rawF(basedex2,:);
    downres = rawF(downdex,:);
    upres = rawF(updex,:);

	%% a quick calculation of dFF
    meanBaseline = mean((baseline1+baseline2)./2,1);
    for j=1:numROIs
        dFF(:,j)=(rawF(:,j)-meanBaseline(j))./meanBaseline(j);
    end
    
    down_dFF = dFF(downdex,:);
    up_dFF = dFF(updex,:);

    %% save variables in folder
    cd(folderpath)
    save('dFF_down.mat','down_dFF')
    save('dFF_up.mat','up_dFF')
    
    save('base1.mat','baseline2')
    save('base2.mat','baseline2');
    
    save('down_res.mat','downres');
    save('up_res.mat','upres');
    
    save('ROImasks.mat','ROI');
    save('downIDX.mat','downdex');
    save('upIDX.mat','updex');

    save('dFF.mat','dFF');
    save('rawF.mat','rawF');
    save('metadata.mat','Xpix','Ypix','FrameRate','Zframes');
end