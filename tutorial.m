% All code here is based on work in the following two papers.
%
% Holmes, C.M. & Nemenman, I.  Estimation of mutual information for
% real-valued data with error bars and controlled bias. 
% Phys. Rev. E, 2019.
%
% A. Kraskov, H. Stogbauer, and P. Grassberger,
% Estimating mutual information.
% Phys. Rev. E 69 (6) 066138, 2004


%%

%%%%%% SECTION ONE %%%%%%

%We want to learn how to estimate information, so let's start wiht a simple
%example, fully written out for you.

%You do not have to write any code in this section; instead, read the
%comments carefully, and make sure that you understand each line. 

%First, let's generate the data for our simple example so that we know
%ahead of time what answer we should get. The easiest thing here is two
%correlated gaussian variables.


%So: generating correlated gaussian variables:
numPoints = 4000; 
correlation = 0.8; %strong correlation

muX = 0;         %X mean of the bivariate gaussian
muY = 0;         %Y mean of the bivariate gaussian
sigmaX = 1;      %X standard deviation   
sigmaY = 1;      %Y standard deviation

%generate a correlated X and Y
x1 = normrnd(0, 1, numPoints, 1);
x2 = normrnd(0, 1, numPoints, 1);
x3 = correlation .* x1 + (1 - correlation^2)^.5 .* x2;
X = muX + x1 * sigmaX;
Y = muY + x3 * sigmaY;

%check what the data looks like:
figure;
plot(X, Y,'.');


%for gaussians, we can explicitly calculate the information. This equation
%holds for all of the cases below that are built on gaussians with a
%correlation; it is a useful point of comparison consistently.
mutualInformationAnalytical = - 1/2* log2(1-correlation^2);

%%

%okay, we are now ready to actually start working with the mutual
%information code.

%First, a line that we only need to run once per install. This compiles the
%C code.
mex MIxnyn.C

%Next, we will set where the temporary saved file will go:
pathToSave = pwd;

%%

%choose a list of k values to try
listOfKs = [1 2 3 4 5 7 10 15 20]; 

%choose a list of data split sizes to use for checking dependence on the size of the data set.
listSplitSizes = [1 2 3 4 5 6 7 8 9 10]; 

%We can run everything without very much coding at all here; this is a one
%line command that will take your data and run it through the pipeline.
%Note the flags 1,1,1, which mean that many plots will be generated while
%this runs.
[values, errorBars] = findMI_KSG_bias_kN(X,Y,listOfKs, listSplitSizes, 1, 1, 1,pathToSave); %runs the actual estimates. the last three entries determine whether or not various plots are generated during the estimation process.
%the output 'values' here corresponds to the estimated mutual information
%for each choice of k, and the ouput 'errorbars' is the corresponding error
%bar.


%note: unfortunately, our C code does everything in nats, not bits (a
%different unit of information). I do not convert back to bits until the
%end:



%before you move on, a few questions to think about:

%1. How much information is there between the variables X and Y?
%2. How does that compare to the analytical result?
%3. Is it a lot? A little? Compare it to the notion of one bit.
%4. Do you understand what each of the plots that came out at line 75 are
%saying, and why we are looking at them?
%5. What value of k would you use for this data?


%%

%%%%%% SECTION TWO %%%%%%

%Now that you understand what you need to do to estimate information, try
%it for yourself! Let's work with a similarly easy data set:


numPoints = 4000; 
correlation = 0.5; %medium strength correlation

muX = 0;         %X mean of the bivariate gaussian
muY = 0;         %Y mean of the bivariate gaussian
sigmaX = 1;      %X standard deviation   
sigmaY = 1;      %Y standard deviation

%generate a correlated X and Y
x1 = normrnd(0, 1, numPoints, 1);
x2 = normrnd(0, 1, numPoints, 1);
x3 = correlation .* x1 + (1 - correlation^2)^.5 .* x2;
X = muX + x1 * sigmaX;
Y = muY + x3 * sigmaY;

%check what the data looks like:
figure;
plot(X, Y,'.');


%for gaussians, we can explicitly calculate the information:
mutualInformationAnalytical = - 1/2* log2(1-correlation^2);
%%

%here, write out the code to estimate information! Ask yourself the same
%questions as you did in the first section.



%%

%%%%%% SECTION THREE %%%%%%

%Once you're sure you understand section two, on to section three! Here, we
%will explore the function reparameterize_data.m

%this function takes data and maps it to a gaussian version of that data,
%while preserving the ordering of the data. 

%to look at this, we need to think about non-gaussian variables. Let's take
%log-linear variables for now, although our method is more general than
%that.

correlation = 0.5;
N = 10000;

x1 = normrnd(0, 1, N, 1);
x2 = normrnd(0, 1, N, 1);
x3 = correlation .* x1 + (1 - correlation^2)^.5 .* x2;
X = [exp(x1)];
Y = [exp(exp(x3))]; 
%%

%We can reparameterize X and Y like this:
X_reparam = reparameterize_data(X);
Y_reparam = reparameterize_data(Y);

%Here, your job is to compare the information estimates for these two
%cases. Note again that the true information remains unchanged!

%Ask ask yourself the same questions that were written out at the end of
%section 1.


%%

%%%%%% SECTION FOUR %%%%%%

%Let's look at that same data again, with far fewer samples. What happens
%now? Are you still able to estimate the information well?

correlation = 0.5;
N = 100;

x1 = normrnd(0, 1, N, 1);
x2 = normrnd(0, 1, N, 1);
x3 = correlation .* x1 + (1 - correlation^2)^.5 .* x2;
X = [exp(x1)];
Y = [exp(exp(x3))]; 

%Estimate the information. What do you need to check in order to know that
%you've done this well?


%%

%%%%%% SECTION FIVE %%%%%%


%Let's look at some real data!

%ISIData.mat comes from Sam Sober's group, and is a cut of the data that
%was featured in the Srivastava, Holmes, et al., PNAS (2017).

%The experiments recorded from a motor unit in breathing anesthetized
%finches. This ISIdata is selected to be phase matched across different
%breathing cycles. We are then looking at two consecutive ISIs (X) and the
%following two consecutive ISIs (Y). Our goal here is to see if there is
%fine timescale structure in the code.

%If you're interested, you could try to replicate figure 1C from the PNAS
%paper with this data. 

%Otherwise, just try to estimate the information! Note that this data is higher
%dimensional. Our code is built to work the same with higher dimensional
%inputs, but how do you expect the higher dimensionality to affect
%information estimation problems?

load('ISIData.mat')

%%

%%%%%% SECTION SIX %%%%%%

%this one is intentionally open ended! 

%bird1.mat is some of the data from bird 1 in the 2017 paper (not entire
%data length for size/convenience reasons)

%here, there are two variables, and you can't just take the mutual
%information between them - you need to create two variables to compare.

%what are you interested in looking at? Can you estimate it well? 

%You can think about relationships between the pressure and spikes, but
%also about structure within each!





