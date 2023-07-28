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

function [rii_c, rii_mode_frac, rii_minmax, r0] = calcRadii(v,t,r, L, ip,pipe_nodes,xtotal,d, int_d, k_r, rho_r, Cp_r, k_f,nu_f,rho_f,Cp_f,eps, pipe_centres, xs1, ys1, zs1, np)
    % compute the r0 used for the R&D method
    % this represents how far the thermal interference
    % would reach into the rock.
    r0 = Calcr0(2, r(ip), L(ip), v(ip), t, k_r, rho_r, Cp_r, k_f,nu_f,rho_f,Cp_f,eps);

    % we compute the critical radii around the pipe
    rii_c = critical_radii(ip,pipe_nodes,xtotal,d,r0,int_d,pipe_centres, xs1, ys1, zs1, np);
    rii_c(rii_c(:) <= r(ip)*1.2) = r(ip)*1.2; % we make sure that the smallest radii is no less than 20% greater than the pipe radius, to allow the FD scheme to work.
    % Computes the radii integral and takes it ratio to r0
    % value, to be used as the weighting factor for the
    % averaging.
    rii_min = min(rii_c);
    if r0 <= rii_min
        rii_mode_frac = 1;
    else
        rii_mode_frac = sum(log(rii_c)-log(r(ip)))/((log(r0)-log(r(ip)))*int_d); % evaluated against OGS models
    end
    rii_minmax = [rii_min, max(rii_c)];
end

