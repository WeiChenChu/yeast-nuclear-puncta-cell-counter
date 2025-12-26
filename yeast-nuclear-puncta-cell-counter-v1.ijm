/*
 FIJI/ ImageJ Macro written by Wei-Chen CHU, ICOB imaging Core, Academia Sinica, Taiwan
 Personal Website: weichenchu.com
 
 Dependency: CLIJ, CLIJ2, PTBIOP, BioVoxxel 3D box, Cellpose3(Python enviornment, important: change to your cellpose enviornment path, in line 68), 
 Last update: 2025-12-26
 
 Count the cell that show green punta within a nuclear membrane in the Yeast
 
 */

#@ File(label = "Input directory", style = "directory") input
#@ File(label = "Output directory", style = "directory") output
#@ String(label = "File suffix", value = ".czi") suffix

run("Colors...", "foreground=black background=white selection=yellow");

//Cell size filter
cell_minimum_size = 1000;
cell_maximum_size = 20000;

//Denoise parameter
median_radius = 3.0;
TopHat_radius = 20.0;

//Signal Detection parameters
red_segment_method = "Otsu";
number_of_dilations_and_erotions = 6.0;

green_signal_prominence = 400;

// Init GPU
run("CLIJ2 Macro Extensions", "cl_device=");
Ext.CLIJ2_clear();

processFolder(input);

// function to scan folders/subfolders/files to find files with correct suffix
function processFolder(input) {
    list = getFileList(input);
    list = Array.sort(list);
    for (i = 0; i < list.length; i++) {
        if (File.isDirectory(input + File.separator + list[i]))
            processFolder(input + File.separator + list[i]);
        if (endsWith(list[i], suffix))
            processFile(input, output, list[i]);
    }
}

function processFile(input, output, file) {
    // Do the processing here by adding your own code.
    // Leave the print statements until things work, then remove them.
    run("Bio-Formats Importer", "open=" + input + File.separator + file + " color_mode=Composite rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
    
    img = getTitle();
    main_file_name = File.nameWithoutExtension;
    
    //Duplicate a cropped image
    
    makeRectangle(315, 0, 1300, 1216);
	
	run("Duplicate...", "duplicate channels=1-3 title=input_img");
		run("Split Channels");

	//Cellpose Cyto3 Segmentation
	//Important: change to your cellpose3 enviornment path!!!!!!!
	selectImage("C3-input_img");
	run("Cellpose ...", "env_path=C:\\Users\\Owner\\miniforge3\\envs\\cellpose3 env_type=conda model=cyto3 model_path=path\\to\\own_cellpose_model diameter=55 ch1=0 ch2=-1 additional_flags=--use_gpu");
	
	
	cell_label = getTitle();
    
    selectImage("C1-input_img");
    img_C1 = getTitle();
    Ext.CLIJ2_push(img_C1);

    selectImage("C2-input_img");
    img_C2 = getTitle();
    Ext.CLIJ2_push(img_C2);

    //Workflow for C2 (Green signal) =========================================

    // Median2D Box
    Ext.CLIJ2_median2DBox(img_C2, img_median, median_radius, median_radius);

    //Ext.CLIJ2_pull(image_3);

    // Subtract Images
    Ext.CLIJ2_subtractImages(img_C2, img_median, img_sub_from_median);

    //Ext.CLIJ2_release(img_C2);
    Ext.CLIJ2_release(img_median);

    Ext.CLIJ2_pull(img_sub_from_median);
    //saveAs("tif", output + File.separator + main_file_name + "_sub_from_median.tif");

    run("Find Maxima...", "prominence=" + green_signal_prominence + " output=[Single Points]");
    run("Subtract...", "value=254");
    green_spot = getTitle();
    setMinAndMax(0, 1);
    //saveAs("tif", output + File.separator + main_file_name + "_spot.tif");

    Ext.CLIJ2_release(img_sub_from_median);
    Ext.CLIJ2_push(green_spot);


    //Workflow for C1 (Red, nuclear membrane)======================

    // Median2D Box
    Ext.CLIJ2_median2DBox(img_C1, img_C1_median, median_radius, median_radius);
    Ext.CLIJ2_release(img_C1);
    //Ext.CLIJ2_pull(img_C1_median);

    // Top Hat Box
    Ext.CLIJ2_topHatBox(img_C1_median, img_C1_median_TopHat, TopHat_radius, TopHat_radius, 0);
    Ext.CLIJ2_release(img_C1_median);
    //Ext.CLIJ2_pull(img_C1_median_TopHat);

    // Automatic Threshold
    Ext.CLIJ2_automaticThreshold(img_C1_median_TopHat, red_mask, red_segment_method);
    //Ext.CLIJ2_release(img_C1_median_TopHat);


    // Closing Diamond
    Ext.CLIJ2_closingDiamond(red_mask, red_mask_closing, number_of_dilations_and_erotions);
    Ext.CLIJ2_addImages(green_spot, red_mask, mask_merge);
    //Ext.CLIJ2_pull(mask_merge);
    //`get`MinAndMax(min, max);
    //setMinAndMax(min, max);
    //run("Fire");

    //Processing cell_labe
    Ext.CLIJ2_push(cell_label);
    Ext.CLIJ2_excludeLabelsOnEdges(cell_label, cell_label_exclude_edge);
    Ext.CLIJ2_excludeLabelsOutsideSizeRange(cell_label_exclude_edge, cell_label_filtered, cell_minimum_size, cell_maximum_size);

    Ext.CLIJ2_maximumIntensityMap(img_C2, cell_label_filtered, green_intensity_map);
    Ext.CLIJ2_pull(green_intensity_map);
    Ext.CLIJ2_excludeLabelsWithValuesWithinRange(green_intensity_map, cell_label_exclude_edge, cell_label_exclude_edge_exclude_overexpose, 16383, 16383);

    Ext.CLIJ2_maximumIntensityMap(mask_merge, cell_label_exclude_edge_exclude_overexpose, positive_map);
    Ext.CLIJ2_pull(positive_map);
    Ext.CLIJ2_excludeLabelsWithValuesOutOfRange(positive_map, cell_label_exclude_edge_exclude_overexpose, positive_label, 2, 2);
    
    //Ext.CLIJ2_pull(positive_map);
    //getMinAndMax(min, max);
    //setMinAndMax(min, max);
    //run("Fire");
    Ext.CLIJ2_pull(cell_label_filtered);
    run("glasbey_on_dark");
    saveAs("tif", output + File.separator + main_file_name + "_label.tif");
    
    Ext.CLIJ2_pull(cell_label_exclude_edge_exclude_overexpose);
    run("glasbey_on_dark");
    getStatistics(area, mean, min, max, std);
    cell_count = max;

    Ext.CLIJ2_pull(positive_label);
    run("glasbey_on_dark");
    getStatistics(area, mean, min, max, std);
    cell_positive_count = max;
	
	if (cell_positive_count>=1) {
		run("Labels to 2D Roi Manager", "process_all_slices=false");
    	roiManager("save", output + File.separator + main_file_name + "_ROI.zip");	
   	    selectImage(img);
    	run("Duplicate...", "duplicate channels=1-3");  
    	roiManager("show all without labels");
	} else {
    	selectImage(img);
    	run("Duplicate...", "duplicate channels=1-3");  
    	//roiManager("show all without labels");
	}
    
    
    Stack.setChannel(1);
	run("Magenta");
    //run("Enhance Contrast...", "saturated=0.35");
    setMinAndMax(500, 4500);
    Stack.setChannel(2);
    //run("Enhance Contrast...", "saturated=0.35");
    setMinAndMax(1000, 7000);
    saveAs("tif", output + File.separator + main_file_name + "_positive_overlay.tif");
    
    //run("Flatten");
    //saveAs("tif", output + File.separator + main_file_name + "_positive_label.tif");
    
    print("----------------------------------------------");
    print(img + "  / cell_count_total     : " + cell_count);
    print(img + "  / cell_count_positive : " + cell_positive_count);
    print(img + "  / Positive_Ratio : " + cell_positive_count / cell_count);

    Ext.CLIJ2_clear();
	run("Close All");

}