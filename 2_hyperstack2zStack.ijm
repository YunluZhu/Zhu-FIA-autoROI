// Go through all the subfolders
	// find hyperstacks, split z slices, save individual z slices


input = getDirectory("Input directory");

print("----");
z = getString("How many slices?: ", "6");

splitHyperstack(input);


function splitHyperstack(input) {
	// this function splits hyperstack into multiple stacks by slice number
    list = getFileList(input);
    for (i = 0; i < list.length; i++) {
        if (startsWith(list[i], "hyperConcat")) {
        	print(input);
            hyperZsplit(input);
            close_all_windows();
            break;
        }
    }
    
    
	for (i = 0; i < list.length; i++) {
        if (endsWith(list[i], "/")) {
//        	children = getFileList(input+list[i]);
//        	if (children.length > 2) {
        		splitHyperstack(input+list[i]);	
//	        }
         }
    }
}
 
function close_all_windows() { 
      while (nImages>0) { 
          selectImage(nImages); 
          close(); 
      } 
} 
  
function hyperZsplit(folder) {
	// 
	// 
	run("Bio-Formats", "open=["+folder+"/hyperConcat_allRegions.tif] color_mode=Default rois_import=[ROI manager] split_focal view=Hyperstack stack_order=XYCZT");
//	Stack.getDimensions(width, height, channels, raw_slices, raw_frames);
	while (nImages>0) {
		fileName = getTitle();
		outFile = folder + "/" + fileName;
		saveAs("Tiff", outFile);
		close();
	}
}



