// ROI recognization
// input directory should be fish folder


input = getDirectory("Input directory containing motion corrected z stacks");

print("----");

process_fishFolder(input);

function process_fishFolder(input) {
	close_all_windows();
    list = getFileList(input);
    flag = 0;
    for (i = 0; i < list.length; i++) {
        if (endsWith(list[i], "_mcor.tif")) {
			flag = 1;
        }
    }
    if (flag == 1) {
    	concatZstack(input);
        getROIs(input);
        print("done");
        close_all_windows();
        print("wd closed");
    }
	for (i = 0; i < list.length; i++) {
        if (endsWith(list[i], "/")) {
        	print(list[i]);
        	process_fishFolder(input+list[i]);	
	     }
	}
}

function concatZstack(fish_folder_dir) {
	z = 0;
	file_list = getFileList(fish_folder_dir);
	for (j = 0; j < file_list.length; j++) {
		if (endsWith(file_list[j], "_mcor.tif")) {
			this_file_dir = fish_folder_dir + file_list[j];
			z = z + 1;
			run("Bio-Formats", "open=["+this_file_dir+"] color_mode=Default rois_import=[ROI manager]  view=Hyperstack stack_order=XYCZT");
			}
	}
	z_num = z;
	Stack.getDimensions(width, height, channels, raw_slices, raw_frames);
	run("Concatenate...", "all_open open");
	run("Stack to Hyperstack...", "order=xytzc channels=1 slices="+z_num+" frames="+raw_frames+" display=Color");
	fish_folder_name_list = split(fish_folder_dir, "/");
	fish_folder_name = fish_folder_name_list[fish_folder_name_list.length-1];
	fish_name_list = split(fish_folder_name, " ");
	fileName = "masterHyperstack" + "_" + fish_name_list[0];
	outFile = fish_folder_dir + fileName;
	saveAs("Tiff", outFile);
	run("Re-order Hyperstack ...", "channels=[Channels (c)] slices=[Frames (t)] frames=[Slices (z)]");
	run("Z Project...", "projection=[Average Intensity] all");
	fileName = getTitle();
	outFile = fish_folder_dir + fileName;
	saveAs("Tiff", outFile);
	close_all_windows();
}	   

function getROIs(fish_folder_dir) {
	close_all_windows();
	file_list = getFileList(fish_folder_dir);
	for (ii = 0; ii < file_list.length; ii++) {
		if (startsWith(file_list[ii], "AVG_master")) {
			this_file_dir = fish_folder_dir + file_list[ii];
			run("Bio-Formats", "open=["+this_file_dir+"]");
			run("Gamma...", "value=0.5 stack");
			setAutoThreshold("Triangle dark no-reset stack");
			run("Convert to Mask", "method=Triangle background=Dark black");
			run("Watershed", "stack");
			run("Analyze Particles...", "size=20-Infinity circularity=0.30-1.00 show=Nothing add stack");
			roiNumber = roiManager("count");
			roiManager("Select", newArray(roiNumber));
			fileName = "RoiSet.zip";
			outFile = fish_folder_dir + fileName;
			roiManager("Save", outFile);
			roiManager("Select", newArray(roiNumber));
			roiManager("Deselect");
			roiManager("Delete");
		}
	}
}

function close_all_windows() { 
      while (nImages>0) { 
          selectImage(nImages); 
          close(); 
      } 
} 
  