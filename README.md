
The following Matlab scripts were developed to analyze light microscopy data, calculate fibril persistence length and quantify ER membrane curvature measurements.

Bäuerlein FJB et al. (2017) Cell



# ER_mobility_analysis_part1.m
# ER_mobility_analysis_part2.m
# ER_mobility_analysis_part3.m

Matlab script in three parts to calculate the accumulation and variance of the endoplasmic reticulum around inclusion bodies (IB).

	•	part 1: loads all ER-movies, identifies individual IBs, defines IB boarders and creates FIJI script for bleach correction
	•	run FIJI script for bleach correction of the ER-movies before continuing with part2
	•	part 2: calculates ER-variance, defines extracellular and nuclear regions to exclude from averaging, calculates rotational averages and saves analyses to output files
	•	part 3: cohort analysis for all IBs found and plotting of radial averages of ER and IB signal plus ER variance 

Necessary data: 

	1)	single exposure in the GFP-channel to identify Htt-GFP inclusion bodies 
	2)	10Hz time series of the ER-channel over 20-30s



# RSDA.m

Ring Shape Detection Algorithm (RSDA).

This Matlab script analyses light microscopy data acquired with the FEI CorrSight microscope using the MAPS-Software:
the script analyses the EM-grid overview taken with x times y tiles, identifies inclusion bodies by detecting ringlike structures that are formed by mCherryUb around the IB and marks the detected IBs in the bright field micrographs by a surrounding rectangle. Finally an Excel list is written, with location, morphological properties, the distance of IBs to the grid bars and a quality score to optimize the selection for FIB lamella preparation.



# Persistencelength_full.m

This Matlab script was developed to calculate the persistence length of fibrils traced in cryo-EM data with Amira 6.0.1, using the XTracing module.

Needed data is the „Attribute“ dataset derived from the XTracing module saved as .xml file. After manually transforming to .xls format, the data is first transformed by the Persistencelength_full.m script with the Amira_lineset_conversion.m script (see below) and by scripts previously developed by M. Jasnin & E. Villa (Jasnin et al. (2013) PNAS) to space the datapoint equidistantly and following the persistence length is calculated with the core script Persistencelength.m.



# Persistencelength.m

This Matlab script is the core of the Persistencelength_full.m script and performes the actual calculation of the persistence length of fibrils.



# Amira_lineset_conversion.m

Matlab script, that aggregates data from the Amira 6.0.1 file.

