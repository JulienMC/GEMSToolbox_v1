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

%% DESCRIPTION
% This function will generate a multi-level grid according to input file
% specifications to be used in GEMSToolbox

function Outputs = step_grid_geometry(PhysicalProperties, NumericalProperties, Outputs)
    % Unpacking inputs
    q_out = PhysicalProperties.q_out;
    q_in = PhysicalProperties.q_in;
    qset = str2num(PhysicalProperties.qset);
    fhn = NumericalProperties.fhn;
    fhn_info = str2num(fhn);
    fhns = fhn_info(:,1); % extract nodes
    fhs  = fhn_info(:,2); % extract corresponding fixed head values to assign

    n  = NumericalProperties.n;   % grid width (number of nodes wide)
    m  = NumericalProperties.m;   % grid height (number of nodes high)
    s  = NumericalProperties.s;   % number of seams
    l1 = NumericalProperties.l1;  % length of horizontal pipes
    l2 = NumericalProperties.l2;  % length of vertical pipes
    h  = NumericalProperties.h;   % interval between seams
    cnx = str2num(NumericalProperties.cnx);%list of connexions between seams
    
    % generate the multi grid geometry
    [A120, xtotal, N, np] = geometry_multigrid(n , m, s, l1, l2, h, cnx);

    % Specify the number of desired fixed heads
    no = length(fhns);
    nn = N - no;

    % creates a map between GIS indices and Matlab code indices
    internalState.Init(nn,no);

    % Set fixed hydraulic heads
    Ho = zeros(no,1); % creats the internal array of fixed head of size no
    for i = 1:no
        Ho = internalState.SetAsFixedHead(fhns(i), fhs(i), Ho); % sets the GIS node fhns(i) to the internal fixed head node, with a value of fhs(i) m
    end

    % Adjust A12, xo, x and A10
    [A12, A10] = internalState.MatSetup(A120);
    [x, xo] = internalState.MatSetup(xtotal, 'inv');

    q = setupwells(qset, q_in, q_out, nn);
    idiagn = nn;
    
    % Bundling variables into struct
    Outputs.('nn') = nn;
    Outputs.('no') = no;
    Outputs.('np') = np;
    Outputs.('A12')= A12;
    Outputs.('A10')= A10;
    Outputs.('xo') = xo;
    Outputs.('x')  = x;
    Outputs.('Ho') = Ho;
    Outputs.('q')  = q;
    Outputs.('idiagn')=idiagn;
    Outputs.('q_out')= q_out;
    Outputs.('q_in')= q_in;
end