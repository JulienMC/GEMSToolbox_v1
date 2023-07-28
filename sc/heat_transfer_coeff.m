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
function h = heat_transfer_coeff (k, r, l, v, nu, rho, Cp, eps)
% This function calculates the heat transfer coefficient for fluid flow in a 
% pipe using correlations for laminar and turbulent flow.
%
% Inputs:
%   k: scalar - Thermal conductivity [W/m.K]
%   r: scalar - Pipe radius [m] 
%   l: scalar - Pipe length [m]
%   v: scalar - Fluid velocity [m/s]
%   nu: scalar - Kinematic viscosity [m^2/s]
%   rho: scalar - Fluid density [kg/m^3]  
%   Cp: scalar - Specific heat capacity [J/kg.K]
%   eps: scalar - Pipe roughness [m]
%
% Outputs:
%   h: scalar - Heat transfer coefficient [W/m^2.K]

% Note that here we use the kinematic viscosity = dynamic viscosity / fluid
% density. Hence, why Re and Pr have the density term inverted.
Re = (2*v*r/nu); % Reynolds number
Pr = (rho*Cp*nu/k); % Prandtl number
f = pipe_friction_factor(Re, r*2, eps); % Darcy weisbach friction factor

%Nu = abs(Dittus_boelter(Re,Pr)); 

if Re > 3000 && Re < 5e6 && Pr > 0.5 && Pr < 2000
    %"Turbulent"
    Nu = Gnielinski(Re,Pr,f);
elseif Re < 3000 && Pr > 5
    %"Laminar Pr>5"
    Nu = Hausen(Re,Pr,r,l);
elseif Re < 3000 && Pr < 5
    %"Laminar Pr<5"
    Nu = SiederTate(Re,Pr,r,l);
else
    %Tubulent "Dittus"
    Nu = Dittus_boelter(Re,Pr);
end
h = abs(Nu)*k/(2*r); % specify absolute value to avoid breaking second law of Thermodynamics

end

% Calculates the Nusslet Number from the Gnielinski correlation
% https://en.wikipedia.org/wiki/Nusselt_number
function Nu = Gnielinski(Re,Pr,f)
    Nu = (f/8)*(Re-1000)*Pr/(1 + 12.7*(f/8)^0.5*(Pr^(2/3)-1));
    if Re < 3000 || Re > 5e6 || Pr < 0.5|| Pr > 2000
        warning("Heat transfer coefficient is outside applicability bounds.\n" + ...
            "Re = %f, should be in range: [3000,5e6].\n" + ...
            "Pr = %f, should be in range: [0.5,2000].\n",Re,Pr)
    end
end

% Calculates the Nusslet Number from the Dittus Boelter equation
% Formally used in the earlier versions
% https://en.wikipedia.org/wiki/Nusselt_number
function Nu = Dittus_boelter(Re,Pr)
    Nu = 0.021*Re^0.8*Pr^0.43;
    if Re < 10000 || Pr < 0.6 || Pr > 160
        warning("Heat transfer coefficient is outside applicability bounds.\n" + ...
            "Re = %f, should be > 10000.\n" + ...
            "Pr = %f, should be > 160.\n", Re,Pr)
    end
end

% Calculate Nusslet Number using Hausen Equation presented in 
% https://www.sciencedirect.com/science/article/pii/S0375650517301955?via%3Dihub#eq0085
% suitable for a laminar regime in both non-fully developed and fully
% developed regions - accoriding to https://web2.clarkson.edu/projects/subramanian/ch330/notes/Heat%20Transfer%20in%20Flow%20Through%20Conduits.pdf
% 1. A.F. Mills, Heat Transfer, Second Edition, Prentice-Hall, New Jersey, 1999.
function Nu = Hausen(Re,Pr,r,l)
    Nu = 3.66 + (0.0688*(2*r/l)*Re*Pr)/(1+0.04*((2*r/l)*Re*Pr)^(2./3));
    if Re > 3000 || Pr < 5
        warning("Heat transfer coefficient is outside applicability bounds.\n" + ...
            "Re = %f, should be < 3000.\n" + ...
            "Pr = %f, should be > 5.\n", Re,Pr)
    end
end

% Calculate Nusslet Number using Sieder Tate Eq Equation presented in 
% https://www.sciencedirect.com/science/article/pii/S0375650517301955?via%3Dihub#eq0085
% suitable for a laminar regime in both non-fully developed and fully
% developed regions where Hansen cannot be applied
function Nu = SiederTate(Re,Pr,r,l)
    Nu = 1.86 * (Re*Pr*2*r/l)^(1/3); % * (mu_b/mu_w)^0.14 - ratios of kinematic viscosities assumed to be 1
    if Re > 3000 || Pr > 5
        warning("Heat transfer coefficient is outside applicability bounds.\n" + ...
            "Re = %f, should be < 3000.\n" + ...
            "Pr = %f, should be < 5.\n", Re,Pr)
    end
end

