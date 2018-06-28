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


function [errorEstimate] = findMI_KSG_stddev(MIs, N, do_plot)
%{
MIs - the output from findMI_twoContDistrubutions
N - total size of the data: how many data points were in the original X
and Y variables? (where we're estimating I(X;Y))

do_plot is a binary variable. If it is 1, plots will automatically be
generated.

errorEstimate - this is a single value, and the estimate of the error at
for the full data size.
%}

Divisions = MIs(:,1);  
MIs = MIs(:,2); %Yes, frustrating labeling, but now MIs only contains the mutual informations, while Divisions contains what we've called 'listSplitSizes' in the past.


n = length(MIs); %we don't need B. n is the same as it was in findMI_twoContDistributions -it's the number of split sizes that we care about. 


Nsamples = zeros(n, 1); 
for i = 1:n
    Nsamples(i) = floor(N./(Divisions{i,1})); %just giving us how many data points were actually in each subset. needed for labeling purposes and not much else.
end

%everything above here should be changed to just be passed from
%findMI_twoContDistributions

Xvals = 1./Nsamples; %because we want to look at dependence on 1/N. 
xActual = Xvals(1);  %'actual' meaning with using the entire available data set
Xvals = Xvals(2:end); %We only want to be doing our fittings for where subsets of the data is used (where there is a nonzero variance)



listVariances = zeros(n - 1,1); %n-1 because there is no variance in the first row, which we assume is where the entire data set was used
for i = 2:n
    listVariances(i - 1) = var(MIs{i,1});  
end


%now need to fit line with slope = 1 to the above: we want the line to be
%of the form that log(variance) = log(Xvals) + b, and we only want to find
%b

%we are only using this line to guide intuition, so we will do this without
%worrying about how to properly weight the points.
b = mean(log10(listVariances) -log10(Xvals));

%create the line we'll plot
Xs = linspace(xActual./2,max(Xvals), 100);
Ys = Xs .* (10^b); %and then we'll plot everything in log space
 

%okay, now let's do the estimation of the error for the actual mutual
%information. We won't do this by projecting along the line, although our
%answer shouldn't be too different from that. 

%instead, we will use the chi^2 calculation to estimate the variance
k = cell2mat(Divisions(2:end)); %because we can't use the first entry here (Divisions = 1, variance = 1) to help us predict the actual variance. 
variancePredicted = sum((k-1)./k.*listVariances)./sum((k-1));
Sml=variancePredicted*N;
varS = 2*Sml^2/sum((k-1));
stdvar = sqrt(varS/N^2); %the error bars on our estimate of the variance


if do_plot == 1
    
figure
hold on
plot(Xvals, listVariances, 'x')
plot(Xs,Ys)

plot(xActual, variancePredicted, 'ks','MarkerFaceColor','k') %draws the projected point
plot([xActual xActual], [variancePredicted - stdvar, variancePredicted + stdvar],'k-')

ax = gca;
ax.XScale = 'log';
ax.YScale = 'log';
xlim([min(Xs) max(Xs)]);

end

errorEstimate = variancePredicted^.5;





end