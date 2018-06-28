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


%{
We're going to start by estimating the mutual information between two
correlated gaussian-distributed variables - we can do this analytically,
which allows us to compare the estimated and actual values. 

This function will move across values of 'k', and will plot an equivalant
to figure 2 at each k. It will also output the dependence of the estimated
mutual information on the value of k.
%}

close all;

numPoints = 200;
muX = 0;
muY = 0;
sigmaX = 1;
sigmaY = 1;
correlation  = 0.5;



%generate a correlated X and Y
x1 = normrnd(0, 1, numPoints, 1);
x2 = normrnd(0, 1, numPoints, 1);
x3 = correlation .* x1 + (1 - correlation^2)^.5 .* x2;
X = muX + x1 * sigmaX;
Y = muY + x3 * sigmaY;


disp(horzcat('correlation = ', num2str(correlation)))


%calculate actual expected mutual information between X and Y
mutualInformation_analytical = log2((1/(1-correlation^2)).^.5);
disp(horzcat('analytical value of the mutual information = ', num2str(mutualInformation_analytical),' bits'))


%now we will do the estimate - we call this fucntion but it will call the
%others. This will generate quite a lot of plots.
listOfKs = [1 2 3 4 5 7 10 15 20];
listSplitSizes = [1 2 3 4 5 6 7 8 9 10];
[values, errorBars] = findMI_KSG_bias_kN(X,Y,listOfKs, listSplitSizes, 1, 1, 1);
%[values, errorBars] = mutualInformation_twoContDistributions_check_k_dependence(X,Y,listOfKs, listSplitSizes);
disp(horzcat('estimated values of mutual information at different k"s are ', num2str(values'), ' and standard deviations are ', num2str(errorBars')))
disp(horzcat('the values of k used are ', num2str(listOfKs)))
 