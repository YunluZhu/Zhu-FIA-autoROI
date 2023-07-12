function[a,ROIs]=fixROIs(a,ROIs)
[x1,y1,n1]=size(ROIs); %get the size of the ROIs X*Y*number of ROIs
for k=1:n1
    test=a{k};%get each individual structure 
    try
    newROI=poly2mask(test.mnCoordinates(:,1),test.mnCoordinates(:,2),x1,y1);%use the real coordinates saved in 'a' to draw the ROI properly using it as a mask
    catch
        disp('well one of your ROIs is not a polygon!!!');
        %add a blank ROI? 
        ROIs(:,:,k)=[];
    end
%mncoordinates gives the outline of the shape, this converts it to a binary
%images
    ROIs(:,:,k)=newROI;%replace each ROI with the new one.
end
    