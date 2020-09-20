Here we provide raw images and codes to support the article.
The complete dataset of raw images is more than 2 Tb. Here we are limited to 50Gb. For the full dataset contact the first author Nicola Pellicciotta.
We choose to provide a full dataset of two culture at DIV 16, one trated with shear flow and a control without flow.
The video with propelled particles are in the directory FL,
 The bright field images without particles are stored in BF, unfortunately only few example because too large in size.
The results of the analysis of this dataset is reported in the direcoty analysis (avaialble for each culture).

Moreover we provide the code to analyse these data.
The routine is to:

Step 1: for each field of view (fov) getting the cilia beating direction from the FL images. This is done with PIV and the code is Step1_PIVanalysis.mat

Step 2: for each fov getting ciliated cell position and CBF from the BF movies. Gather the cilia beating direction and cilia posion and frequency in a unique figure and matlab class (Res.mat). this is done in Step2_gatherResults.mat

The results of these analysis are stored in the analysis folder for each culture.

These routines are repeated for each experiment and results are then plotted to get trends. In the folder code4figures we report the code that we used to make the figures in the papers starting from a matlab file "all_results*.mat", where are gathered all the analysis. 

The code may improve in the future, with more comments. please check my Nicola's github page for the latest update.Please contact us for any problem.

