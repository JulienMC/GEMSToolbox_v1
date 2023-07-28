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
function r = res_coef_Brown03(Q,L,d, eps, rho_f, nu_f)
    % This function calculates the resistance coefficient (r) of a pipe 
    % using the Darcy-Weisbach formulation as described in:
    %
    % Brown, S.F., 2002. The history of the Darcy-Weisbach equation for pipe flow resistance
    % https://ascelibrary.org/doi/10.1061/40650\%282003\%294 https://ascelibrary.org/doi/abs/10.1061/40650\%282003\%294
    % The resistance coefficient relates pressure drop and flow rate:
    %    dP = r * Q
    % Where:
    %   r = f * 8/(pi^2*g) * L/d^5 
    %   f = Darcy-Weisbach friction factor 
    %   L = Pipe length [m]
    %   d = Pipe diameter [m] 
    %   g = Gravity [m/s^2]
    %
    % Inputs:
    %   Q: Scalar - Volumetric flow rate [m^3/s] 
    %   L: Scalar - Pipe length [m]
    %   d: Scalar - Pipe diameter [m]
    %   eps: Scalar - Pipe roughness [m] 
    %   rho_f: Scalar - Fluid density [kg/m^3]
    %   nu_f: Scalar - Fluid kinematic viscosity [m^2/s]
    %
    % Outputs:
    %   r: Scalar - Resistance coefficient
    %   
    % The friction factor f is calculated using the pipe_friction_factor function.
    %

   % Calculating resistance coeff of a pipe, using the Darcy-Weisbach
   % formulation as described in Brown, 2003:
   %    r=f*8/(pi^2*g)*L/d^5
   %            f    = the friction factor (p.189 of EPAnet manual)
   %            d    = diameter of the pipe (m)
   %            L    = length of the pipe (m)
   %            g    = gravitational acceleration 
   % 
   % Note that this calculation is based on values/formulae for pipes
  
   mu  = nu_f*rho_f;  % dynamic viscosity of water (Pa s)
   % Reynolds number = rho*v*d/mu, with v=Q/(pi*r^2) = 4*Q/(pi*d^2)
   %    (for flow in pipe, see https://en.wikipedia.org/wiki/Reynolds_number)
   Re  = rho_f*abs(4*Q/(pi*d^2))*d/mu;
   f = pipe_friction_factor(Re, d, eps);
   g = 9.8;
   r = f*8/(pi^2*g)*L/d^5;
   %fprintf("Debug: Q, = %e, Re = %e, f = %e, r = %e\n",Q,Re,f,r);
end