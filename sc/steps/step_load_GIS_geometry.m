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
% This is the generic geometry function used to load a specific ArcGIS
% shapefile into GEMSToolbox

function Outputs = step_load_GIS_geometry(PhysicalProperties, NumericalProperties, Outputs)
    % Unpacking inputs
    ArcGIS_file_name = NumericalProperties.ArcGIS_file_name;
    q_out = PhysicalProperties.q_out;
    q_in = PhysicalProperties.q_in;
    qset = str2num(PhysicalProperties.qset);
    fhn = NumericalProperties.fhn;
    fhn_info = str2num(fhn);
    fhns = fhn_info(:,1); % extract nodes
    fhs  = fhn_info(:,2); % extract corresponding fixed head values to assign
    cnx = str2num(NumericalProperties.cnx); %list of connexions between model nodes
    depths = PhysicalProperties.depths; %list of depths to apply to different GIS files (if no POINT_Z value is specified)
    reuse_geom = NumericalProperties.reuse_geom;
    
    if reuse_geom == 0
        [A120, xtotal, N, np] = build_arc_geometry({ArcGIS_file_name}, cnx, depths);
    else
        warning("Re-using previous geometry in run no %d", Outputs.irun);
        A120 = Outputs.A120_toreuse;
        xtotal = Outputs.xtotal_toreuse;
        N = Outputs.N_toreuse;
        np = Outputs.np_toreuse;
    end

 
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

    % added to be on the safe side as otherwise mineflow takes ages!
    A12 = sparse(A12);
    A10 = sparse(A10);
    
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

    % store variable for potential reuse of the geometry
    Outputs.('A120_toreuse') = A120;
    Outputs.('xtotal_toreuse') = xtotal;
    Outputs.('N_toreuse') = N;
    Outputs.('np_toreuse') = np;
end