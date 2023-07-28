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

function [Tout, r0] = pipeheat (r, l, Tin, v, t, Tr, k_r, rho_r, Cp_r, k_f, nu_f,rho_f, Cp_f, eps)
% Description:
% This function calculates the temperature change (Tout) in a pipe segment due to heat exchange with the pipe wall
% and computes the radius of the cooling 'halo' (r0). The method used for the calculations is based on
% the work by Rodriguez and Diaz (2009).
%
% Parameters:
%   r: Radius of the pipe (meters).
%   l: Length of the pipe segment (meters).
%   Tin: Initial temperature of the fluid inside the pipe (Celsius).
%   v: Fluid velocity inside the pipe (m/s).
%   t: Time duration of heat exchange (seconds).
%   Tr: Temperature of the rock surrounding the pipe (Celsius).
%   k_r: Thermal conductivity of the rock (W/m·K).
%   rho_r: Density of the rock (kg/m³).
%   Cp_r: Specific heat capacity of the rock (J/kg·K).
%   k_f: Thermal conductivity of the fluid (W/m·K).
%   nu_f: Kinematic viscosity of the fluid (m²/s).
%   rho_f: Density of the fluid (kg/m³).
%   Cp_f: Specific heat capacity of the fluid (J/kg·K).
%   eps: Darcy-Weisbach roughness (unitless).
%
% Returns:
%   Tout: The final temperature of the fluid inside the pipe segment after heat exchange (Kelvin).
%   r0: The radius of the rock surrounding the pipe after heat exchange (meters).
%
% Note:
% The function uses the 'heat_transfer_coeff' function to calculate the heat transfer coefficient for the fluid.
% The value of 'verbose' controls whether additional information is displayed during the calculations.

verbose = 0;

VF = pi*r^2*v;    % volumetric flow rate = cross_pipe surface * velocity (m^3/s)

% fluid heat_transfer_coeff
h_f = heat_transfer_coeff (k_f, r, l, v, nu_f, rho_f, Cp_f, eps);

% Update r0 = radius of rock surrounding pipe that cooled using one of two methods:
% 1) Original method from (Rodriguez and Diaz, 2009), Equation A17.
% 2) Improved version of this:
%     a) Revised calculation of total ground heat release Equation A12
%     b) Avoid assumption in r0-derivation that T_pi = (T0+T_e)/2
imethod = 2;
if imethod == 1
    r0 = r*sqrt(1+4*h_f/(rho_r*Cp_r*r)*t);
elseif imethod == 2
    r0_in = 2*r;
    niter=1; nitermax=100;
    while 1
        niter = niter + 1;
        if niter>=nitermax
           fprintf('WARNING: pipeheat.m, r0 calc., imethod=2: no convergence in %d iterations\n',niter);
           break;
       end
       L = log(r0_in/r); % Ask JvH for doc showing derivation
       r0_out =  r * sqrt(1 + 4*k_r/(rho_r*Cp_r*r^2)*t + L);
       if abs(r0_in - r0_out)/r0_in<1e-5
           break;
       end
       r0_in = r0_out;
    end
    r0 = r0_out;
end
if verbose
   fprintf(' pipeheat: t=%f, r0=%f\n',t,r0)
end

% If the volumetric flow rate in the pipe is negligeable, then the following computation
% exceeds the rock temperature. So we set the Tout to rock temp and exit
% function here. (we don't do it before as we want r0 to compute.)
VF_th = 1e-13;
if VF < VF_th
    Tout = Tr;
    return;
end

% effective heat transfer coeff for fluid + wall:
U = (1/h_f + r/k_r*log(r0/r))^-1;
% Outflow T using eqn 4 in Rodriguez & Diaz (2009):
coef1 = 2*pi*r*l*U*Tr;
coef2 = rho_f*Cp_f*VF;
coef3 = pi*r*l*U;
Tout = (coef1 + (coef2-coef3)*Tin) / (coef2+coef3);

if verbose
    L = log (r0/r);
    Tp = (k_r/L*Tr+r*h_f*Tout)/(k_r/L + r*h_f);
    fprintf('   Tw=%f, Tp=%f, T0=%f\n',Tout, Tp, Tr)
end

