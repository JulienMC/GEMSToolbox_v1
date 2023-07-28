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
function r0 = Calcr0(imethod, r, l, v, t, k_r, rho_r, Cp_r, k_f,nu_f,rho_f,Cp_f,eps)
% % Loading physical properties
% k_r = PhysicalProperties.k_r;       % rock heat conductivity (W/m,K)
% rho_r = PhysicalProperties.rho_r;   % rock density (kg/m^3)
% Cp_r = PhysicalProperties.Cp_r;     % rock specific heat (J/kg,K)
% k_f = PhysicalProperties.k_f;       % water heat conductivity (W/m,K)
% nu_f = PhysicalProperties.nu_f;     % water kinematic viscocity (m^2/s)
% rho_f = PhysicalProperties.rho_f;   % water density (kg/m^3)
% Cp_f = PhysicalProperties.Cp_f;     % water specific heat (J/kg,K)
% eps = PhysicalProperties.eps;       % Darcy-Weisbach roughness

% fluid heat_transfer_coeff
h_f = heat_transfer_coeff (k_f, r, l, v, nu_f, rho_f, Cp_f, eps);

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
end