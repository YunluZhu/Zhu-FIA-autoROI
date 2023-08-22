// Go through all the subfolders
	// concatenate a list of stacks over Z&T (OME Tiff from 2P) into one hyperstack
	// concatenate individual tiff 2p data into one hyperstack
	// save the hyperstack
	// zNumber is the number of intended slices in the flattened OME TIFF
// FISH LEVEL DIR recommended
input = getDirectory("Input directory");
z = getString("How many slices?: ", "6");

print("----");
 
processFolder(input);

function processFolder(input) {


	flag=0;
	list_of_oriStack = newArray(1);
    list = getFileList(input);
    for (i = 0; i < list.length; i++) {
        if (startsWith(list[i], "Image_scan_")) { //if OME tiffs
        	print(input);
			print(list[i]);
//            stack2hyperZprojTproj(input, list[i], z);
//            close_all_windows();
            if (flag==0) {
            	list_of_oriStack[flag] = input+list[i];
            }
            else {
            	newFileDir = newArray(1);
            	newFileDir[0] = input+list[i];
            	list_of_oriStack = Array.concat(list_of_oriStack, newFileDir);
            }
            flag = 1;
        }
        ///////
        else if (startsWith(list[i], "ChanA_0")) { //if individual tiffs
        	if (indexOf(input, "anatomy") == -1) {
	        	print(input);
	            seq2hyperZprojTproj(input, z);
	            close_all_windows();
	            break;
        	}
        }
        
        ///////
    }
    if (flag > 0) {
    	stackList2hyperConcat(input, list_of_oriStack, z);
    	close_all_windows();
    	flag = 0;
    	list_of_oriStack = newArray(1);
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
  
function stack2hyperZprojTproj(folder, imgName, zNumber) {
	// convert an OME TIFF to hyperstack then perform z projection. save file
	// zNumber is the number of intended slices in the flattened OME TIFF
	run("Bio-Formats", "open=" + "[" + folder + imgName + "]");
	zNumber = parseInt(zNumber);
	Stack.getDimensions(width, height, channels, raw_slices, raw_frames);
	frameNumber = raw_frames / zNumber;
	frameNumber = round(frameNumber);
//	open();
	run("Stack to Hyperstack...", "order=xyczt(default) channels=1 slices="+zNumber+" frames="+frameNumber+" display=Color");
//	outFile = folder + "hyper_" + imgName;
//	outFile = substring(outFile, 0, lengthOf(outFile)-1);
//	saveAs("Tiff", outFile);
	run("Z Project...", "projection=[Max Intensity] all");
	outFile = folder + "max_" + imgName;
	outFile = substring(outFile, 0, lengthOf(outFile)-1);
	saveAs("Tiff", outFile);
    close();
//    run("Re-order Hyperstack ...", "channels=[Channels (c)] slices=[Frames (t)] frames=[Slices (z)]");
//    run("Z Project...", "projection=[Max Intensity] all");
//	run("Re-order Hyperstack ...", "channels=[Channels (c)] slices=[Frames (t)] frames=[Slices (z)]");
//	outFile = folder + "tmax_" + imgName;
//	outFile = substring(outFile, 0, lengthOf(outFile)-1);
//	saveAs("Tiff", outFile);
}

function stackList2hyperConcat(folder, list_of_imgDir, zNumber) {
	// concatenate a list of stacks over Z&T (OME Tiff from 2P) into one hyperstack
	// save the hyperstack
	// max projection over T
	// save tmax
	// zNumber is the number of intended slices in the flattened OME TIFF
	list_of_imgDir = Array.sort(list_of_imgDir);
	for (j=0; j < list_of_imgDir.length; j++) {
		run("Bio-Formats", "open=" + "[" + list_of_imgDir[j] + "]");
		zNumber = parseInt(zNumber);
		Stack.getDimensions(width, height, channels, raw_slices, raw_frames);
		frameNumber = raw_frames / zNumber;
		frameNumber = round(frameNumber);
		run("Stack to Hyperstack...", "order=xyczt(default) channels=1 slices="+zNumber+" frames="+frameNumber+" display=Color");
	}
	run("Concatenate...", "all_open open");
	outFile = folder + "hyperConcat_allRegions";
	saveAs("Tiff", outFile);
//	run("Re-order Hyperstack ...", "channels=[Channels (c)] slices=[Frames (t)] frames=[Slices (z)]");
//    run("Z Project...", "projection=[Max Intensity] all");
//	run("Re-order Hyperstack ...", "channels=[Channels (c)] slices=[Frames (t)] frames=[Slices (z)]");
//	outFile = folder + "tmax_" + "hyperConcat_allRegions";
//	saveAs("Tiff", outFile);

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
//	run("Z Project...", "projection=[Max Intensity] all");
//	run("Stack Splitter", "number=4");
//	outFile = folder + "max_Image_scan_4_region_0_0";
//	saveAs("Tiff", outFile);
//	close();
//	outFile = folder + "max_Image_scan_3_region_0_0";
//	saveAs("Tiff", outFile);
//	close();
//	outFile = folder + "max_Image_scan_2_region_0_0";
//	saveAs("Tiff", outFile);
//	close();
//	outFile = folder + "max_Image_scan_1_region_0_0";
//	saveAs("Tiff", outFile);
//	close();
	// close max proj
//	close();
	// get t max for each z
//	run("Re-order Hyperstack ...", "channels=[Channels (c)] slices=[Frames (t)] frames=[Slices (z)]");
//	run("Z Project...", "projection=[Max Intensity] all");
//	run("Re-order Hyperstack ...", "channels=[Channels (c)] slices=[Frames (t)] frames=[Slices (z)]");
//	outFile = folder + "tmax_" + "hyperConcat_allRegions";
//	saveAs("Tiff", outFile);
}

