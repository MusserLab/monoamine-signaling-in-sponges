/*This script analyze IF data of sponges treated with various chemicals and stained for Vinculin.
 * Number of cells is identified via counting DAPI objects.
 * Vinculin signal will be thresholded 
 * Obtain number of vinculin foci and size of foci
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

