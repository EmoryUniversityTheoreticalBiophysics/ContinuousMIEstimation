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
% Holmes, C.M. & Nemenman, I.  Estimation of mutual information for
% real-valued data with error bars and controlled bias. 
% Submitted, 2018.
%
% A. Kraskov, H. Stogbauer, and P. Grassberger,
% Estimating mutual information.
% Phys. Rev. E 69 (6) 066138, 2004
%
% in your published research.


function [transformedValues] = reparameterize_data(listOfValues)

%function [transformedValues] = reparameterize_data(listOfValues)
%
%This function reparamaterizes the data to a standard normal distribution. 
%
% listOfValues -- one dimensional real valued array of data. For example, 
% it could be one of the components of the multidimensional X variable 
% in our information estimation. This reparamaterization, if being
% used, should be performed on each dimension independently.

%This reparamaterization does not change the ordering of the points -
%instead, it just changes the distances along the axis.


[A,I] = sort(listOfValues); % data solrted by their values and their index
            % in the original data set

listPositions = zeros(length(listOfValues),1);

for i = 1:length(listOfValues)
    listPositions(I(i)) = i;
end

%reassign each i'th value to where the i'th value would have been in a 
%standard normal dataset using the erfinv function.

listPositions = listPositions-mean(listPositions);
listPositions = 2.*listPositions./(length(listPositions));
transformedValues = sqrt(2)*erfinv(listPositions);
