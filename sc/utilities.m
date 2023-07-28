% This program allows for the computation of water and heat flow through
% a mine network
%     Copyright (C) 2022  Durham University
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <https://www.gnu.org/licenses/>.

classdef utilities
    % An static class which contains utility functions

    %   Methods:
    %   1. perlin_noise(im) - Generates Perlin noise and adds it to the input image 'im'.
    %       Input: 'im' - Input image of size (n x m) to which Perlin noise will be added.
    %       Output: 'im' - Output image with added Perlin noise.
    %
    %   2. normalize(array, amax, amin) - Normalizes the input array between 'amin' and 'amax'.
    %       Input: 'array' - Input array to be normalized.
    %              'amax'  - Maximum value for normalization.
    %              'amin'  - Minimum value for normalization.
    %       Output: 'normx' - Normalized array with values scaled between 0 and 1.
    %
    %   3. distance(n1, n2) - Calculates the Euclidean distance between two points n1 and n2.
    %       Input: 'n1' - Coordinates of the first point (2-element vector: [x1, y1]).
    %              'n2' - Coordinates of the second point (2-element vector: [x2, y2]).
    %       Output: 'x' - Euclidean distance between points n1 and n2.

    methods(Static)
        function im = perlin_noise(im)
        
            [n, m] = size(im);
            i = 0;
            w = sqrt(n*m);
            
            while w > 3
                i = i + 1;
                d = interp2(randn(n, m), i-1, 'spline');
                im = im + i * d(1:n, 1:m);
                w = w - ceil(w/2 - 1);
            end
        end

        function normx = normalize(array, amax, amin)
            normx = (array - amin) / (amax - amin);
        end

        function x = distance(n1, n2)
            x = sqrt((n2(1)-n1(1))^2 + (n2(2)-n1(2))^2);
        end
    end
end

