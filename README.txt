======BASIC DESCRIPTION:=====

This software implements the modified Kraskov-Stoegbauer-Grassberger (KSG)
estimator of mutual information between two continuous variables, with
modifications introduced by Holmes and Nemenman in

1. Holmes, CM. & Nemenman, I.  Estimation of mutual information for
real-valued data with error bars and controlled bias. Submitted, 2018.

Please cite the reference [1] above when using this software. The software
is COPYRIGHTED by Holmes and Nemenman, 2018, and is distributed under
GPL 3. Please read the license information in the attached file.

The original KSG estimator was introduced in

2. Kraskov, A., Stoegbauer, H., & Grassberger, P. (2004). Estimating
mutual information. Phys Rev E, 69(6), 066138.

3. Stoegbauer, H., Kraskov, A., Astakhov, S.A., & Grassberger,
P. (2004), Least Dependent Component Analysis Based on Mutual
Information. Physical Review E, 70, 066123.

The code for the KSG estimator was copyrighted and distributed under
GPL 3 in 2009, by Astakhov, Kraskov, Stoegbauer, and Grassberger. The
original KSG code is distributed with no changes together  with this
software for ease of use and completeness. Please cite references [2]
and [3] when using this software as well.

All comments, suggestions, and proposed modifications should be
communicated to the authors:
	     Caroline M. Holmes <cholmes@princeton.edu>
	     Ilya Nemenman <ilya.nemenman@emory.edu>

The main GitHub repository for this software is:
    	     https://github.com/EmoryUniversityTheoreticalBiophysics/ContinuousEntropyEstimation
The software can be downloaded and modified there.

The software requires a C compiler with standard libraries, as well as
a newer version of Matlab (we have run it using Matlab R2016b, but
newer and older versions may also work).


=====CONTENT OF THE PACKAGE=====

The contents of this folder are as follows:

1. A folder titled ‘kraskovStoegbauerGrassberger’. This contains the
files that perform the nearest-neighbor mutual information estimation
using the KSG estimator (see Refs. [2,3]). The files are distributed
with no change and copyrighted by their original authors. You will need
to do a few things in this folder before you can use the estimator:

	a. Open the file ‘MIxnyn.m’ in this directory, edit the file
	    and change '[INSERT PATH HERE]' in the file with the path
	    to the directory where the current software package is
	    installed.
	    
	b. In terminal, navigate to the kraskovGrassberger folder and
	    compile the code, using something like:
              gcc -c miutils.C -o miutils.o
	      gcc MIxnyn.C -o MIxnyn miutils.o
	
2. The .m files described in the text of Ref. [1]. For examples of how
to run these files on data, see the ExampleScripts below.

	findMI_KSG_subsampling.m - [CMH: BRIEF DESCRIPTION]
	findMI_KSG_stddev.m - [CMH: BRIEF DESCRIPTION]
	findMI_KSG_bias_kN.m - [CMH: BRIEF DESCRIPTION]
	reparamaterize_data.m - [CMH: BRIEF DESCRIPTION]

3. A folder titled ‘ExampleScripts’, which contains scripts to perform
analyses similar to those in the figures of [1].

	a. gaussianExample.m - this function performs mutual information
	    estimation on correlated gaussian data. It generates
	    figures similar to Figs. 2, 3 in Ref. [1], which show how
	    error bars are calculated, and how mutual information
	    estimates depends on k and on N.
	
	b. logNormalExample.m - this shows how error bars are
	    calculated and how mutual information estimates depend on
	    k and N for log-normal bivariate date. It also
	    demonstrates the effects of reparameterizing the data and
	    generates the equivalent of Fig. 4 in Ref. [1].

	c. higherDimensionalExample.m - This script performs analyses
	    similar to (a) and (b) but with higher dimensional
	    Gaussian inputs. It is equivalent to Fig. 6 in Ref. [1].

	d. NfkBDataAnalysis.m -- [CMH: INCLUDE AND DESCRIBE]
	
	e. NfkBData.mat -- this is the data file for data used in
            Fig. 7 in Ref. [1]. The data lists single cell NF-kappaB
            activation (P65 nuclear localization) for in response to
            dozes of the TNF stimulus. See Ref. [1] and

	    	  Cheong, R., Rhee, A., Wang, C. J., Nemenman, I., &
	    	  Levchenko, A. (2011). Information transduction
	    	  capacity of noisy biochemical signaling
	    	  networks. Science 334(6054), 354–358.

	    for additional details about the data.
		  
	f.  finchDataAnalysis.m -- [CMH: INCLUDE AND DESCRIBE]

	g. finchData.mat - this file contains the data used by
	    finchDataAnalysis.m.  [CMH: DESCRIBE THE DATA, REFERENCE
	    EXPERIMENTAL PAPER]
	
