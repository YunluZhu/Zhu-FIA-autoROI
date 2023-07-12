function [Xpix,Ypix,FrameRate,Zframes]=findkeyparams_vol()

%use this function to get the FOV and the framerate of aquisition 
%you have to parse the XML file from thor. 
%this sucks and is definitely not elegent 
matexppath = pwd;
matexpfile = 'Experiment.xml';
cd(matexppath);
xmlexp=parseXML5D('Experiment.xml');
kids=xmlexp.Children;
% for x= 1:length(kids)-1
%     disp(kids(x).Name);
% end
FrameRate=kids(26).Attributes(23).Value; %LSM variable, attribute 11 is framerate
Xpix=kids(26).Attributes(36).Value;%hardcoded because I'm not good at this
Ypix=kids(26).Attributes(37).Value;
Xpix=str2double(Xpix);
Ypix=str2double(Ypix);
Zframes=kids(18).Attributes(6).Value;
Zframes=str2double(Zframes);

end