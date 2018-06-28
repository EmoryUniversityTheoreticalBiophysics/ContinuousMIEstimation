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


function [MIs] = findMI_KSG_subsampling(X,Y, kvalue, listSplitSizes, do_plot)
%{

uses kraskov-grassberger to estimate I(X;Y)
X,Y - continuous 1-d variables - should each be given as vectors. These
*must* be the same length as each other

kvalue - number of nearest neighbors to count when using the kg estimator.

listSplitSizes - divisions in the data when looking at the dependence of
the mutual information on size. For example, 2 means that we should be
looking at two non-overlapping subsets of the data, each of size N/2. The
first entry should always be 1.

do_plot is a binary variable. If it is 1, plots will automatically be
generated.

MIs - a cell array, with n rows, where each row contains data from a different
amount of the data used for each point. There are two columns - the left entries contain the
relevant splitSizes for that row, and the right entries contain the
estimated mutual informations for each of the subdivisions of the data.

%}

n = length(listSplitSizes); 
MIs = cell(n,2);

N = length(X);

means = zeros(n,1);
errorBars = zeros(n,1);

for i = 1:n
    MIs{i,1} = listSplitSizes(i);
end

for i = 1:n
    
    %randomly choose which subset each data point belongs to:
    a = randperm(length(X)); 
    k = listSplitSizes(i); %number of divisions we need the data to be split into
    l = round(linspace(0,length(X),k + 1));
    
    MI_T = zeros(k,1); %MI_T is just a temporary place to hold the different mutual information values for each combination subset
    
    for j = 1: k %move through each of the subsets of equal size
        xT = X(a(l(j) + 1:l(j+1))); % the 'X' values that belong to this subset
        yT = Y(a(l(j) + 1:l(j+1))); % the 'Y' values that belong to this subset
        MI_T(j) = MIxnyn(xT,yT,kvalue); %call the estimator function.
    end
    MIs{i,2} = MI_T./log(2); %we need this log(2) to convert from nats to bits.
    means(i) = mean(MIs{i,2});
    errorBars(i) = std(MIs{i,2});
end

if do_plot == 1
figure
errorbar(listSplitSizes./N, means, errorBars)
xlabel('1/N')
ylabel('mutual information, in bits')
title(horzcat('mutual information dependence on N, for k = ', num2str(kvalue)))
end


end