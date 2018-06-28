When using this code, please reference:


Kraskov, A., Stögbauer, H., & Grassberger, P. (2004). Estimating
mutual information. Phys Rev E, 69(6), 066138. http://doi.org/10.1103/PhysRevE.69.066138

Harald Stogbauer, Alexander Kraskov, Sergey A. Astakhov and Peter Grassberger, Least Dependent Component Analysis Based on Mutual Information. Physical Review E, 70, 066123 (2004).

(And our paper).

This code includes files originally copyrighted in 2009, by Sergey Astakhov, Alexander Kraskov, Harald Stoegbauer, and Peter Grassberger.




The contents of this folder are as follows:
1. A folder titled ‘KraskovStoegbauerGrassberger’. This contains the files that perform the nearest-neighbor mutual information estimation that we have termed the KSG estimator. We are here redistributing code originally distributed by the authors of the estimator. You will need to do a few things in this folder before you can use the estimator:
	a. Open the file ‘MIxnyn.m’, and edit the file and change '[INSERT PATH HERE]' with the path to the directory where the downloaded code is
	b. In terminal, navigate to the kraskovGrassberger folder and compile the code, using something like:
	   gcc -c miutils.C -o miutils.o
	   gcc MIxnyn.C -o MIxnyn miutils.o
	
2. The .m files described in the text of Holmes & Nemenman, 2018. For examples on how to run these on data, see the ExampleScripts, which demonstrate how to use these.
	findMI_KSG_subsampling.m
	findMI_KSG_stddev.m
	findMI_KSG_bias_kN.m
	reparamaterize_data.m

3. A folder titled ‘ExampleScripts’, which contains scripts to perform similar analyses to the figures in Holmes & Nemenman 2018.  
	a. gaussianExample.m - this script performs mutual information estimation on correlated gaussian data. It generates figures similar to figure 2 and figure 3 in Holmes & Nemenman, which show how error bars are calculated, how mutual information estimates depends on k, and how they depend on N.
	b. logNormalExample.m - much like the previous script, this shows how error bars are calculated and how the estimates depend on k. It also demonstrates the effects of renormalizing the data. This generates the equivalent to figure 4 in Holmes & Nemenman.
	c. higherDimensionalExample.m - This script performs a very similar analysis to gaussianExample.m, but with higher  dimensional inputs. Equivalent to Fig. 6 in Holmes & Nemenman.
	d. finchDataAnalysis.m -- not currently included. This one will probably just exactly generate the corresponding figure.
	e. finchData.mat - this contains the data used by finchDataAnalysis.m. 
	f. NfkBDataAnalysis.m -- not currently included. This one will probably just exactly generate the corresponding figure.
	g. NfkBData.mat

4. (Do we want to put this in?) A folder of scripts that will exactly generate each figure in Holmes and Nemenman 2018.  I would include here probably the script that generates the data for a figure, the data itself, and then the script that actually generates the figure. 


**** Question: do we want to include the various code by KSG that they distributed with this that is entirely irrelevant?

	