/*This script analyzes vinculin foci number normalized by DAPI cell count from IF data
 * of sponges treated with various chemicals.
 * Number of cells is identified via counting DAPI objects.
 * Vinculin signal will be thresholded.
 * Obtain number of vinculin foci normalized by number of cells.
 */

// 1. Get the directory containing your TIFF files
inputDir = getDirectory("Choose Source Directory ");

// Get a list of all files in that directory
list = getFileList(inputDir);

// Turn on Batch Mode (runs much faster, hides images)
setBatchMode(true);

// Clear the Results table before starting so you get fresh data
run("Clear Results");
run("Set Measurements...", "area mean redirect=None decimal=3");

// 2. Start the Loop
for (i = 0; i < list.length; i++) {
    
    // Check if the file is a TIFF (ignores system files like .DS_Store)
    if (endsWith(list[i], ".tif") || endsWith(list[i], ".tiff")) {
        
    // Open the image
    open(inputDir + list[i]);

	file_name = getTitle();
	run("Split Channels");
	
	// Counting number of cells
	selectImage("C4-" + file_name);
	resetMinAndMax;
	run("8-bit");
	run("Gaussian Blur...", "sigma=5");
	run("Subtract Background...", "rolling=100");
	
	
	setAutoThreshold("Default dark");
	run("Convert to Mask");
	run("Watershed");
	run("Analyze Particles...", "size=10-Infinity circularity=0.10-1.00 display exclude clear summarize overlay");

	// 6. Cleanup
	//setBatchMode(false);
	  
	// Optional: Print count to log
	DAPI_count = nResults;
	
	selectImage("C1-" + file_name);
	close("\\Others");
	
	resetMinAndMax;

	// Manual correction for two images with aberrant intensity range
	if (file_name == "MAX_20251102_B22p18_P1022_TryptamineIF_Dish4_DAPI_ATub488_Pha568_EmVinc647_40xW_sponge3_zstack_tent-02-Airyscan Processing-52_Z80-537.tif" ||
	    file_name == "MAX_20251102_B22p18_P1022_TryptamineIF_Dish4_DAPI_ATub488_Pha568_EmVinc647_40xW_sponge3_zstack_tent-01-Airyscan Processing-51_Z280-663.tif") {
	        setMinAndMax(50, 250);
	}

	run("8-bit");
	run("Gaussian Blur...", "sigma=2");
	run("Subtract Background...", "rolling=10");
	
	setAutoThreshold("Default dark");
	run("Convert to Mask");
	run("Watershed");
	run("Analyze Particles...", "size=0.5-Infinity display exclude clear summarize overlay");
	
	Vinc_count = nResults;
	close("*");

    }
}

// Turn off Batch Mode to see the final tables
setBatchMode(false);

// Confirmation message
showMessage("Batch Processing Complete!");

