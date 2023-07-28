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

clear

iterations = 1;

output = ones(iterations,7)*-999;
% Any unfilled output spaces will equal -999


for i=1:iterations
    % Produce random values for variables within their respective ranges
    [kr_in, Cp_in, rho_in, diameter_in, Tr_in] = random_input_mod();
    
    % Input these into mine_geothermal_ensemble
    [Tout, Qout] = mine_geothermal_ensemble_mod(kr_in, Cp_in, rho_in,...
        diameter_in, Tr_in);
    
    % Save outflow temperature and flow rate
    output(i,:) = [Tout, Qout, kr_in, Cp_in, rho_in,...
        diameter_in, Tr_in];
end

disp('')
disp('1st col: Output temperature, deg C')
disp('2nd col: Flow rate, L/s')
disp('3rd col: Rock heat conductivity (W/m,K)')
disp('4th col: Rock specific heat (J/kg,K)')
disp('5th col: Rock density (kg/m^3)')
disp('6th col: Pipe diameters, m')
disp('7th col: Rock initial temperature, degC')
disp(output)
    