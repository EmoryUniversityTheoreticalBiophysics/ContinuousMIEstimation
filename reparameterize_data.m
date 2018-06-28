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


function [transformedValues] = reparameterize_data(listOfValues)

%This reparamaterizes our data to a gaussian. listOfValues should be
%one-dimensional. for example, it could be the first dimension of the X
%variable in our information estimation. This reparamaterization, if being
%used, should be performed on each dimension independently.

%This reparamaterization does not change the ordering of the points -
%instead, it just changes the distances along the axis so that you have a
%gaussian distribution. 

[A,I] = sort(listOfValues); 

listPositions = zeros(length(listOfValues),1);

for i = 1:length(listOfValues)
    listPositions(I(i)) = i;
end

%reassign to be in range where can use erfinv: need to be between -1 and 1.


listPositions = listPositions-mean(listPositions);
listPositions = 2.*listPositions./(length(listPositions));

transformedValues = erfinv(listPositions);
end
