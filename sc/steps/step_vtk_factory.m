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

% A function which is used to generate VTK output files compatible with
% Paraview.
% Is able to handle both Point and Element data for 2D elements (i.e. pipe)
% Can handle scalars and vector fields.
% Currentl all data is converted into a double.
% Julien Mouli-Castillo 
% v20222202 - initial version
function step_vtk_factory(PhysicalProperties, NumericalProperties, Outputs)
    verbose = 0;

    filename = "results/"+Outputs.sc_name+"/"+Outputs.vtkf;
    timestep = Outputs.irun;
    time = PhysicalProperties.nyrs;
    pipe_nodes = Outputs.pipe_nodes;
    xtotal = Outputs.xtotal;
    n_data_struct = Outputs.n_outputs;
    e_data_struct = Outputs.e_outputs;

    % Reshapes the input cell arrays 
    n_data = {n_data_struct{:,1}};
    n_dat_names = string({n_data_struct{:,2}});
    e_data = {e_data_struct{:,1}};
    e_dat_names = string({e_data_struct{:,2}});
    
    % takes in a filename, the geometry data as pipe_nodes, an array of node data, and an array of element
    % data.

    fname = strcat(filename, sprintf("_%05d%", timestep), ".vtk");
    fid = fopen(fname, 'w');

    %%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Define Headers %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%

    %File header
    header_fmt = "# vtk DataFile Version 3.0\nUnstructured Grid from GEMSToolbox\n" + ...
        "ASCII\n" + ...
        "DATASET UNSTRUCTURED_GRID\n" + ...
        "FIELD TimesAndCycles 2\n" + ...
        "TIME 1 1 double\n" + ...
        "%.12e\n" + ... % Time as 0.000000000000e+00
        "CYLCE 1 1 long\n" + ...
        "%d\n" + ... % Timestep as int
        "POINTS %d double\n";

    % point coords in format: 0.000000000000e+00 0.000000000000e+00 0.000000000000e+00

    cell_head_fmt = "CELLS %d %d\n"; %the cell header with number of cells and size = cells * 3

    % then print the cells "cell_number input_node output_node\n"

    cell_fmt = "2 %d %d\n"; %the cell format

    cell_type_head_fmt = "\nCELL_TYPES %d\n"; %the cell type code for each cell
    
    % number of points, data_name
    point_data_head_fmt = "\nPOINT_DATA %d\n"; % number of points
    
    %either scalars
    scalar_data_fmt = "SCALARS %s double 1\nLOOKUP_TABLE default\n"; % followed by fmt 7.677728988728e+15
    % or vectors
    vector_data_fmt = "VECTORS %s_ double\n"; % followed by vector value in folowing fmt "-1.357478178396e-05 -1.932580868139e-07 0.0 ""

    % Cell data followed by similar formating for Scalars or Vectors
    cell_data_head_fmt = "CELL_DATA %d\n"; % number of cells

    %%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Write to file %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf("writing VTK output to file: %s ...\n", fname)
    num_points = length(xtotal);
    num_ele = length(pipe_nodes);

    fprintf(fid, header_fmt, time, timestep, num_points);

    % Point coords
    for n = 1:num_points
        x = full(xtotal(n,1));
        y = full(xtotal(n,2));
        s = size(xtotal);
        if s(2) > 2
             z = full(xtotal(n,3));
        else
            z = 0;
        end
        
        point_coord_fmt = '%.12e %.12e %.12e\n';
        fprintf(fid, point_coord_fmt, x, y, z);
    end

    %Cell header
    fprintf(fid, cell_head_fmt, num_ele, num_ele*3);
    %Cell definition
    for p = 1:length(pipe_nodes)
        % -1 because ParaView starts indexing at 0
        in = pipe_nodes(p,1)-1;
        out = pipe_nodes(p,2)-1;
    
        fprintf(fid, cell_fmt, in, out);
    end

    %Cell type header
    fprintf(fid, cell_type_head_fmt, num_ele);

    % Cell type
    for p = 1:length(pipe_nodes)       
        cell_type_fmt = '3\n'; % 1: VERTEX, 2: VERTEX POLY, 3: LINE, 4: POLY LINE, 5: TRIANGLE, 9: QUAD (more at https://www.princeton.edu/~efeibush/viscourse/vtk.pdf Figure 2) 
        fprintf(fid, cell_type_fmt);
    end
    
    %Point & Cell Data
    datapacks = {{n_data, n_dat_names, point_data_head_fmt, num_points}, {e_data, e_dat_names, cell_data_head_fmt, num_ele}};
    for i = 1:2
        dp = datapacks{i};
        data = dp{1};
        dat_names = dp{2};
        data_head_fmt = dp{3};
        points = dp{4};
        if ~isempty(data)
            fprintf(fid, data_head_fmt, points);
            for d = 1:length(data)
                vals = cell2mat(data(d));
                if verbose == 1 fprintf("...writing %s %s ...", strcat(string(size(vals))), dat_names(d)); end
                % determines wether data is vector or scalar type
                type = 0;
                if min(size(data{d})) == 1
                    type = 1;
                    fprintf(fid, scalar_data_fmt, dat_names(d));
                elseif min(size(data{d})) == 2
                    type = 2;
                    fprintf(fid, vector_data_fmt, dat_names(d));
                else
                    type = 3;
                    fprintf(fid, vector_data_fmt, dat_names(d));
                end
    
                %vals(vals <= 1e-15) = 0;                
                for v = 1:length(vals)      
                    if type == 1
                        fprintf(fid, "%.12e\n", vals(v));
                    elseif type == 2
                        fprintf(fid, "%.12e %.12e 0.000000000000e+00\n", vals(v,:));
                    else
                        fprintf(fid, "%.12e %.12e %.12e\n", vals(v,:));
                    end
                end
                if verbose == 1 fprintf("complete\n"); end
            end
        end
    end
    fclose(fid);
    fprintf("... VTK file successfully completed.\n")
