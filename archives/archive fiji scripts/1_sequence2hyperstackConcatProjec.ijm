// Go through all the subfolders
	// concatenate a list of stacks over Z&T (OME Tiff from 2P) into one hyperstack
	// save the hyperstack
	// max projection over T
	// save tmax
	// zNumber is the number of intended slices in the flattened OME TIFF


input = getDirectory("Input directory");
z = getString("How many slices?: ", "6");

print("----");
 
processFolder(input);

function processFolder(input) {

	list_of_oriStack = newArray(1);
    list = getFileList(input);
    for (i = 0; i < list.length; i++) {
        if (startsWith(list[i], "ChanA")) {
        	print(input);
            seq2hyperZprojTproj(input, z);
            close_all_windows();
            break;
        }
    }
	for (i = 0; i < list.length; i++) {
        if (endsWith(list[i], "/")) {
        	processFolder(input+list[i]);
        }
    }
    	
}
 
function close_all_windows() { 
      while (nImages>0) { 
          selectImage(nImages); 
          close(); 
      } 
} 
  
function seq2hyperZprojTproj(folder, zNumber) {
	// 
	// 
	File.openSequence(folder, "virtual filter=ChanA_0");
	Stack.getDimensions(width, height, channels, raw_slices, raw_frames);
	frameNumber = raw_slices / zNumber;
	frameNumber = round(frameNumber);
	run("Stack to Hyperstack...", "order=xytzc channels=1 slices="+zNumber+" frames="+frameNumber+" display=Color");
	outFile = folder + "hyperConcat_allRegions";
	saveAs("Tiff", outFile);
	// get z max for each response
//	run("Duplicate...", "duplicate");
	run("Z Project...", "projection=[Max Intensity] all");
	run("Stack Splitter", "number=4");
	outFile = folder + "max_Image_scan_4_region_0_0";
	saveAs("Tiff", outFile);
	close();
	outFile = folder + "max_Image_scan_3_region_0_0";
	saveAs("Tiff", outFile);
	close();
	outFile = folder + "max_Image_scan_2_region_0_0";
	saveAs("Tiff", outFile);
	close();
	outFile = folder + "max_Image_scan_1_region_0_0";
	saveAs("Tiff", outFile);
	close();
	// close max proj
	close();
	// get t max for each z
	run("Re-order Hyperstack ...", "channels=[Channels (c)] slices=[Frames (t)] frames=[Slices (z)]");
	run("Z Project...", "projection=[Max Intensity] all");
	run("Re-order Hyperstack ...", "channels=[Channels (c)] slices=[Frames (t)] frames=[Slices (z)]");
	outFile = folder + "tmax_" + "hyperConcat_allRegions";
	saveAs("Tiff", outFile);
}
