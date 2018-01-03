# Examples of Running the Pipeline

Here we demonstrate three walkthroughs of how to run the MOSFIRE pipeline. We include a longslit reduction, a slitmask reduction, and a long2pos reduction in the H and K bands. Example datasets can be downloaded from this link, and includes all three data types.

## Getting Started

After downloading and unzipping all of the test data, make a directory in your preferred location to perform your reduction and run handle.

    mkdir reduced
    cd reduced
    mospy handle /home/[yourhomedir]/Data/test_dataset/*.fits
	
You should see five new directories after handle is done.
	
	LONGSLIT-3x0.7 <-- Longslit observations
	LONGSLIT-46x0.7 <-- Longslit calibrations
	MOSFIRE_DRP_MASK <-- Slitmask calibrations and observations
	long2pos <-- long2pos calibrations
	long2pos_specphot <-- long2pos_specphot observations

## Longslit Reduction

Move to the Longslit observation directory and copy the calibrations to the observation directory:

	cd LONGSLIT-3x0.7/2012nov28/K
	cp ../../../LONGSLIT-46x0.7/2013oct15/K/*.txt ./
	
Run the autodriver to create the driver file Longslit_HIP17971.py

