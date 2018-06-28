% Copyright 2018 Caroline Holmes, Ilya Nemenman
%-----------------------------------------------------------------------------------------
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should receive a copy of the GNU General Public License
% along with this program.  See also <http://www.gnu.org/licenses/>.
%-----------------------------------------------------------------------------------------
% Please reference
% 
%(OUR PAPER)
%
% A. Kraskov, H. Stogbauer, and P. Grassberger,
% Estimating mutual information.
% Phys. Rev. E 69 (6) 066138, 2004
%
% in your published research.


function [means, errorBars] = findMI_KSG_bias_kN(X,Y,listOfKs, listSplitSizes, do_plot, do_plot_stddev, do_plot_subsampling)
%{
uses kraskov-grassberger to estimate I(X;Y)
X,Y - continuous 1-d variables - should each be given as vectors. These
*must* be the same length as each other

listOfKs - k is the number of nearest neighbors to count when using the kg
estimator. this is a vector of the different k values that we would like to
test.

listSplitSizes - divisions in the data when looking at the dependence of
the mutual information on size. For example, 2 means that we should be
looking at two non-overlapping subsets of the data, each of size N/2. The
first entry should always be 1.

means - values of estimated I(X;Y) at each k value in listOfKs

errorBars - errors for the estimated values of I(X;Y) at each k

do_plot is a binary variable. If it is 1, plots will automatically be
generated.

do_plot_stddev and do_plot_subsampling are equivalents for
findMI_KSG_stddev and findMI_KSG_subsampling
%}

means = zeros(length(listOfKs),1); %the estimated value of the mutual information at each k
errorBars = zeros(length(listOfKs),1); %the errorbars at each k
N = length(X);

for i = 1:length(listOfKs) %go through each k value, find the values of the information and the error estiates
    [MIs] = findMI_KSG_subsampling(X,Y, listOfKs(i), listSplitSizes, do_plot_subsampling);
    means(i) = mean(MIs{1,2});
    errorBars(i) = findMI_KSG_stddev(MIs, N, do_plot_stddev); 
end


%plot this to see the relationship between k, N, and the estimates of the
%information. 

if do_plot == 1
figure
errorbar(listOfKs, means, errorBars)
xlabel('k')
ylabel('mutual information, in bits')
title('mutual information across k values')
end

end