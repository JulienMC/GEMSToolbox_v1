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
function Outputs = step_mineflow_pipes(PhysicalProperties, NumericalProperties, Outputs)
% Description:
% This routine solves for the flow in the pipes. The method follows the Todi
% EPANET manual (Appendix D and p.30/Table 3.1), which follows (Todini&Paliti, 1987)
% The 'minor loss coefficient m in EPANET App D is assumed 0
%
% Parameters:
%   PhysicalProperties: Structure containing physical properties of the mine network.
%     - eps: Darcy-Weisbach roughness (meters).
%     - rho_f: Density of the fluid (kg/m³).
%     - nu_f: Kinematic viscosity of the fluid (m²/s).
%
%   NumericalProperties: Structure containing numerical properties for the flow computation.
%     - sumdQrel_th: Threshold for the relative head residuals to achieve convergence.
%     - maxdQrel_th: Threshold for the relative flow residuals to achieve convergence.
%     - a: Oscillation damping parameter (0 < a < 1, 'a' and 'b' together should equal 1).
%     - b: Oscillation damping parameter (0 < b < 1, 'a' and 'b' together should equal 1).
%     - num_doubles_th: Threshold to detect stagnating computation, stops after that number of repeat values.
%
%   Outputs: Structure containing outputs from previous iterations.
%     - nn: Number of unknown nodes.
%     - np: Number of pipes.
%     - A12: Connectivity matrix of unknown nodes.
%     - A10: Connectivity matrix of fixed head nodes.
%     - q: flow demand vector (typically contains the well flow rates in m^3.s^-1).
%     - L: Pipe lengths (meters).
%     - d: Pipe hydraulic diameters (meters).
%     - Ho: Known nodal heads.
%
% Returns:
%   Outputs: Updated structure containing the results of the flow calculation.
%     - H: Head distribution in the network after the flow calculation.
%     - Ho: Known nodal heads.
%     - Q: Flow distribution in the pipes after the flow calculation.
%     - Qdiff: Difference in flow values between iterations.
%     - f1: Residual vector for head calculation.
%     - f2: Residual vector for flow calculation (including a zero error for fixed head nodes).


%%% Internal parameters to solve for flow:
% Set friction coef arguments
eps = PhysicalProperties.eps;
rho_f = PhysicalProperties.rho_f;
nu_f = PhysicalProperties.nu_f;
% Set numerical error thresholds
sumdQrel_th = NumericalProperties.sumdQrel_th;
maxdQrel_th = NumericalProperties.maxdQrel_th;
% Set Oscillation Damping params
a = NumericalProperties.a; % a+b must equal 1
b = NumericalProperties.b; % a+b must equal 1
% Set stagnation threshold
num_doubles_th = NumericalProperties.num_doubles;
% Set outputs from previous steps
nn = Outputs.nn;
np = Outputs.np;
A12= Outputs.A12;
A10= Outputs.A10;
q  = Outputs.q;
L  = Outputs.L;
d  = Outputs.d;
Ho = Outputs.Ho;

%%% Internal parameters to solve for flow:
% Set initial guess for heads H and flow Q
H    = zeros(nn,1);
Q    = ones(np,1);
% Fluid resistance coefficient (called R in
r    = zeros(np,1);
% flow powerlaw exponent ( called B in EPANET manual, p.30, and n in
%                          (Todini&Paliti, 1987))
B    = 2;  % Not sure if eqns below are valid if B is not 2!

A21  = A12'; % Make sure both A21, A12 and A10 are all sparse matrices for speed!

% Start iterative calculation using Newton method:
f1_rel=1e10; % sets very large error to start with
f2_rel=1e10; % sets very large error to start with
counter = 1;
sumdQrel_arr = zeros(1,1);
sign_flip = 0; % counts the number of sign changes in the error
ips = 1:np; % this is the i and j vector to build the sparse matrices A11 and Dinv
A11v  = zeros(1,np); % this is the v vector to build the sparse matrix A11
Dv    = zeros(1,np);
Dinvv = zeros(1,np); % this is the v vector to build another sparse matrix Dinv
while (f1_rel>sumdQrel_th || f2_rel>maxdQrel_th) %original option was 1e-8
    % XXXXXXX  Todini and Rossman 2013 - Q and H computation XXXXX
    for ip=1:np
        %% Open Pipes
        % calc resistance coeff using Darcy-Weisbach formula
        %    (EPANET manual p.30 & Table 3.1)
        r(ip) = res_coef_Brown03(Q(ip),L(ip),d(ip), eps, rho_f, nu_f);
        % Todini & Rossman, 2013, Eq 17
        A11v(1,ip) = (r(ip)*abs(Q(ip))^(B-1));
        % Todini & Rossman, 2013, Eq 40
        Dv(1,ip) = (B*r(ip)*abs(Q(ip))^(B-1));
        Dinvv(1,ip) = Dv(1,ip)^-1;
    end
    % Generate matrices from the vectors
    A11  = sparse(ips,ips,A11v);
    Dinv = sparse(ips,ips,Dinvv);
    D    = sparse(ips,ips,Dv);
    %Todini & Rossman 2013, Eq 41 & 42 (note q is negative for us as we use
    %a different convention
    A = A21*Dinv*A12;
    F = A21 * Dinv * ((D-A11)*Q - A10*Ho) - q;
    Hnew = A\F;
    Qnew = Q - Dinv*(A11*Q+A12*Hnew+A10*Ho);
    %% XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    % Replace zero flow values with very small nonzero values
    % If Q contains 0 flows this results in divide-by-0 errors
    Qnew(Qnew==0) = min(Qnew(Qnew>0))/10;

    % Calculate difference in Q between iterations:
    avQ = sum(abs(Qnew))/length(Qnew);
    dQrel = abs((Q-Qnew)./avQ);
    sumdQrel = sum(dQrel);

    % Stores the change in flow for the iteration
    Qdiff = abs(Q-Qnew);

    % prepare for next iteration:
    Q = a*Q+b*Qnew; % damping oscillation - Q reset

    % Residual
    f1 = A11*Q + A12*Hnew + A10*Ho;
    f2 = A21*Q -q;

    % Relative Error
    f1_rel = sum(abs(f1))/(max(Q)-min(Q));
    f2_rel = sum(abs(f2))/(max(Hnew)-min(Hnew));

    counter = counter + 1;
    sumdQrel_arr(counter,1) = sumdQrel;

    % check solution - if stuck in infginite loop the model quits
    num_doubles = 0; % count the number of double numbers in sumdQrel_arr
    for icounter = 1:counter-1
        if round(sumdQrel_arr(counter),5,'significant') == round(sumdQrel_arr(icounter),5,'significant')
            num_doubles = num_doubles + 1;
            if num_doubles > num_doubles_th % value to ensure double numbers are not due to chance
                error('Error in flow calculation. Stuck in infinite loop. Consider changing a and b in Q reset in mineflow.m. Decrease a and increase b.');
            end
        end
    end

    % checks for oscillating convergence
    if counter > 2 && sign(sumdQrel_arr(counter)-sumdQrel_arr(counter-1)) ~= sign(sumdQrel_arr(counter - 1)-sumdQrel_arr(counter - 2))
        sign_flip = sign_flip + 1;
        if sign_flip > 10 % value to ensure sign flips numbers are not due to chance
            %error('Error in flow calculation. Stuck in oscillations.');
            abs_err  = abs((Hnew-min(Hnew))-(H-min(H)));
            warning('Flow convergence is oscillating. Computation stopped after %d iterations, at a sum(Head_residual)/Head_range = %e, and um(Flow_residual)/Flow_range = %e, results might be innacurate.\n', counter, f2_rel, f1_rel);
            break;
        end
    end

end
% Normalise the final heads so that no negative heads exist (minimum head
% set to 0)
H=Hnew-min(Hnew);
Ho=Ho-min(Hnew);

% Bundle outputs variable to struct
Outputs.('H') = H;
Outputs.('Ho')= Ho;
Outputs.('Q') = Q;
Outputs.('Qdiff') = Qdiff;
Outputs.("f1") = f1;
Outputs.("f2") = [f2; 0]; % 0 error at fixed head node
 end
