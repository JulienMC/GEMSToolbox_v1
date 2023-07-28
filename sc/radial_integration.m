%% This program allows for the computation of water and heat flow through a mine network
%%     Copyright (C) 2022  Durham University
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
%%
% A scaling function to approximate the critical radius of a pipe in 3D
% r_cs = (2,1); contains the min and max critical radii for that pipe
% r_ca = (2,1); contains the matching angles at which these critical radii
% are observed.
% Density = int; is the density of the integration a density of 100 offers
% a stable solution.
% r_as = array; all the radii at 360 angle around the pipe
function r_as = radial_integration(r_cs, a_cs, r_inf, density)
%density = 1000; for accurate results
verbose = 0;
r_cmin = r_cs(1);
r_cmax = r_cs(2);
a_cmin = a_cs(1);
a_cmax = a_cs(2);
a_mean = (a_cmin+a_cmax)/2;
a_mean_pi = a_mean+pi;
angles = a_cmin:(2*pi)/(density-1):(2*pi)+a_cmin;%linspace(0,2*pi,density)+ a_cmin;
r_as = nan(density,1); % to store the representative radii at various integraiton angles

    for i = 1:density
        a = angles(i);
        if (a >= a_cmin && a <= a_mean) || (a > a_mean_pi && a <= a_cmin + 2*pi)
            % region where to use scaled r_cmin
            r_as(i) = abs(r_cmin / cos(a - a_cmin));
        elseif (a > a_mean && a <= a_cmax) || (a > a_cmax && a <= a_mean_pi) 
            % region where to use scaled r_cmax
            r_as(i) = abs(r_cmax / cos(abs(a_cmax - a)));
        end
        if (a > a_cmax + pi/2 && a < a_cmin + 3/2 * pi && a_cmax - a_cmin < pi) || (a > a_cmin + pi/2 && a < a_cmax - pi/2 && a_cmax - a_cmin > pi) 
            % define the infinite heat supply region when the interference
            % angles are not at 180 degrees.
            r_as(i) = r_inf;
        end
        r_as(i) = min(r_as(i),r_inf);
    end
end