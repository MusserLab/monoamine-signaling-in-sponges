/*This script measures phospho-myosin light chain (pMYL) fluorescence intensity
 * from IF data of sponges treated with various chemicals.
 * Images are segmented using percentile thresholding and Analyze Particles.
 * Mean fluorescence intensity is measured per ROI and saved as CSV.
 */

// 1. Get the directory containing your TIFF files
inputDir = getDirectory("Choose Source Directory ");

outDir = getDirectory("Output Directory");

// Get a list of all files in that directory
list = getFileList(inputDir);

// Turn on Batch Mode (runs much faster, hides images)
//setBatchMode(true);

// Clear the Results table before starting so you get fresh data
run("Clear Results");
run("Set Measurements...", "area mean display redirect=None decimal=3");

// 2. Start the Loop
for (i = 0; i < list.length; i++) {
    
    // Check if the file is a TIFF (ignores system files like .DS_Store)
    if (endsWith(list[i], ".tif") || endsWith(list[i], ".tiff")) {
        
    // Open the image
    open(inputDir + list[i]);
    
	file_name = File.nameWithoutExtension;
	
	run("Split Channels");
	close;
	close;
	
	run("Gaussian Blur...", "sigma=2");
	setAutoThreshold("Percentile dark no-reset");
	
	run("Convert to Mask");
	run("Analyze Particles...", "size=1-Infinity add");
	wait(5000);
	
	close;
	
	roiManager("Deselect");
	roiManager("Measure");

	savePath = outDir + file_name + "_pMYL_MFI.csv";
	saveAs("Results", savePath);
	
	close("*");
	run("Clear Results");
	roiManager("Delete");
	
	}
}

// Turn off Batch Mode to see the final tables
//setBatchMode(false);

// Confirmation message
showMessage("Batch Processing Complete!");
	