%% this function is called by analyzeAndCombine, so give the experiment header to that function and then this one will work for each individual trial

function [numROIs] = analyzeROIs_volHyper(folderpath,foldername,parentDir)

    % folderpath = "/Volumes/LabDataPro 1/2P nMLF speed/analyzed/230609 fish1/area1_pre/f1-3_from5-10-20-30_005";
    % foldername = "f1-3_from5-10-20-30_005";
    % parentDir = "/Volumes/LabDataPro 1/2P nMLF speed/analyzed/230609 fish1/area1_pre";
    % cd(folderpath)
    %% Extract key para
    [Xpix,Ypix,FrameRate,Zframes] = findkeyparams_vol();
    FrameRate = str2double(FrameRate);
    % calculate hypothetical T all
    Tall_hyp1 = floor(FrameRate/(Zframes+1)*20)*4;  
    Tall_hyp2 = round(FrameRate/(Zframes+1)*20)*4;
    
    %% Import ROI
    height = Ypix;
    width = Xpix;
    
    ROI_file='RoiSet.zip'; %path to ROI file (.zip)

    % Need ReadImageJROI.m script
    [a,ROI,zPos] = ReadImageJROI_vol(ROI_file,[height,width]);%need to reformat a using poly2mask
    [a,ROI] = fixROIs(a,ROI); % ROI are the binary masks for which pixels to analyze vs ignore; is a d1xd2xnumberofROIs matrix

    numROIs=size(ROI,3);
    % cd(folderpath)  % go back to the parent folder (area folder) to find the ori Tiff file

    tmp = TIFFStack('hyperConcat_allRegions.tif');  % TIFFStack() is used here to read hyperstack
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
    dur_dex = round(Tall/26.6);
    sti_dur = Tall/4;
    basedex1 = (sti_dur-dur_dex):sti_dur;
    basedex2 = (sti_dur*2-dur_dex):sti_dur*2;
    basedex3 = (sti_dur*3-dur_dex):sti_dur*3;
    basedex4 = (sti_dur*4-dur_dex):sti_dur*4;
    res1dex = 1:(sti_dur-dur_dex);
    res2dex = sti_dur+1:(2*sti_dur-dur_dex);
    res3dex = sti_dur*2+1:(3*sti_dur-dur_dex);
    res4dex = sti_dur*3+1:(4*sti_dur-dur_dex);

    
    baseline = rawF(cat(2,basedex1,basedex2,basedex3,basedex4),:);
    res1 = rawF(res1dex,:);
    res2 = rawF(res2dex,:);
    res3 = rawF(res3dex,:);
    res4 = rawF(res4dex,:);

	%% a quick calculation of dFF
    meanBaseline = mean(baseline,1);
    for j=1:numROIs
        dFF(:,j)=(rawF(:,j)-meanBaseline(j))./meanBaseline(j);
    end
    
    dff1 = dFF(res1dex,:);
    dff2 = dFF(res2dex,:);
    dff3 = dFF(res3dex,:);
    dff4 = dFF(res4dex,:);

    %% save variables in folder
    cd(folderpath)
    save('dFF1.mat','dff1')
    save('dFF2.mat','dff2')
    save('dFF3.mat','dff3')
    save('dFF4.mat','dff4') 

    save('raw1.mat','res1')
    save('raw2.mat','res2')
    save('raw3.mat','res3')
    save('raw4.mat','res4') 

    save('base.mat','baseline')

    save('ROImasks.mat','ROI');

    save('dFF.mat','dFF');
    save('rawF.mat','rawF');
    save('metadata.mat','Xpix','Ypix','FrameRate','Zframes');
end