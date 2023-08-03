// Go through all the subfolders
	// find hyperstacks, split z slices, save individual z slices
	// input directory needs to be higher or equal to Fish level folder

input = getDirectory("Input directory");

print("----");
z = getString("How many slices?: ", "6");


process_combRep(input)
process_combArea(input)

function process_combRep(input) {
    list = getFileList(input);
    for (i = 0; i < list.length; i++) {
        if (startsWith(list[i], "f")) {
//        	print(input);
            concatRep_byZ(input);
            close_all_windows();
            break;
        }
    }
    
	for (i = 0; i < list.length; i++) {
        if (endsWith(list[i], "/")) { 	
        	process_combRep(input+list[i]);	
       		close_all_windows();
         }
	}
}

function concatRep_byZ(area_folder) {
	rep_list = getFileList(area_folder);
	for (zz=0; zz<z; zz++) {
		for (j = 0; j < rep_list.length; j++) {
			if (startsWith(rep_list[j], "f")) {

				this_file_dir = area_folder + rep_list[j] + "hyperConcat_allRegions.tif - Z=" + zz + ".tif";
				run("Bio-Formats", "open=["+this_file_dir+"] color_mode=Default rois_import=[ROI manager] split_focal view=Hyperstack stack_order=XYCZT");
			}
		}
		if (nImages>0) {
			print(area_folder);
			print(rep_list[0]);
			run("Concatenate...", "all_open open");
			fileName = "z"+zz+"_allRep";
			outFile = input + fileName;
			saveAs("Tiff", outFile);
			close();
		}
		
	}	   
}


function process_combArea(input) {
    list = getFileList(input);
    for (i = 0; i < list.length; i++) {
        if (indexOf(list[i], "area") >= 0) {
            concatArea_byZ(input);
            close_all_windows();
            break;
        }
    }
    
	for (i = 0; i < list.length; i++) {
        if (endsWith(list[i], "/")) {
        	process_combArea(input+list[i]);	
        	close_all_windows();
         }
	}
}
	
	


function concatArea_byZ(fish_folder) {
	area_list = getFileList(fish_folder);
	for (zz=0; zz<z; zz++) {
		for (a = 0; a < area_list.length; a++) {
			if (startsWith(area_list[a], "area")) {
				this_file_dir = fish_folder + area_list[a] + "z" + zz + "_allRep.tif";
				run("Bio-Formats", "open=["+this_file_dir+"] color_mode=Default rois_import=[ROI manager] split_focal view=Hyperstack stack_order=XYCZT");
			}
		}
		if (nImages>1) {
			print(fish_folder);
			print(area_list[0]);
			run("Concatenate...", "all_open open");
			fileName = "z"+zz+"_allArea";
			outFile = input + fileName;
			saveAs("Tiff", outFile);
			close();
		}
		else if (nImages==1) {
			print(fish_folder);
			print(area_list[0]);
			fileName = "z"+zz+"_allArea";
			outFile = input + fileName;
			saveAs("Tiff", outFile);
			close();
		}
	}	   
}

function close_all_windows() { 
      while (nImages>0) { 
          selectImage(nImages); 
          close(); 
      } 
} 
  