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

% Calculates heat transfer through a planar opening experiencing a water
% flow. See Loredo et al. 2017 for details.

function Tout = heatplate_PB(Tin, Tr, x, P, q, c_w, lambda_r, rho_r, c_r, t, A, rho_w)
    % Calculate T change relative to Tr (negative for cooling, positive for
    % heating)

    % Inputs:
    %   Tin: scalar - Inflow water temperature [C] 
    %   Tr: scalar - Background rock temperature [C]
    %   x: scalar - Distance from injection [m]
    %   P: scalar - Effective perimeter of the flow area [m] 
    %   q: scalar - Volumetric flow rate [m^3/s]
    %   c_w: scalar - Specific heat of water [J/kg.C]
    %   lambda_r: scalar - Thermal conductivity of rock [W/m.C]
    %   rho_r: scalar - Density of rock [kg/m^3]
    %   c_r: scalar - Specific heat of rock [J/kg.C]  
    %   t: scalar - Time [s]
    %   A: scalar - Cross-sectional area of flow [m^2]
    %   rho_w: scalar - Density of water [kg/m^3]
    %
    % Outputs:
    %   Tout: scalar - Outflow water temperature [C]

    dT = Pruess_Bodvarsson(Tin, Tr, x, P, q, c_w, lambda_r, rho_r, c_r, t, A, rho_w);
    Tout = dT + Tr;
end

function dT = Pruess_Bodvarsson(Tin, Tr, x, P, q, c_w, lambda_r, rho_r, c_r, t, A, rho_w)
    T1 = x*P/(2*q*c_w);
    T2 = ((lambda_r*rho_r*c_r)/(t-x*A*rho_w/q))^0.5;
    U = isreal(T1 * T2);
    dT = (Tin - Tr) * erfc(T1*T2)*U;
end