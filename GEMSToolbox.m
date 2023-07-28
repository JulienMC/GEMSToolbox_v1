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

%% GEMSToolbox Main script to run mine geothermal model
%
% This is the main script that executes the key steps of GEMSToolbox
% including loading inputs, generating geometry, 
% calculating flow, computing temperatures, post-processing outputs,
% exporting results, and saving well temperatures.
%
% Outputs are bundled into a struct 'Outputs' that is passed between 
% functions. Helper functions exist for vectorizing data, assigning 
% visualizations, and dynamically calling model steps.
%
% The model is run over multiple tests defined in a scenario input file.

clear
clc

% More (debug?) output? --> Set verbose to 1.
verbose = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% PRE-PROCESSING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('sc'); % adds the source code folder to the workspace
addpath('sc/steps'); % adds the source code folder to the workspace
% Load the scenario
fid = fopen('inputfiles/_InputFileName.csv', 'r');
scenarioName = textscan(fid, '%s', 1);
scenarioPath = "inputfiles/"+cellstr(scenarioName{1});
%scenarioPath = "inputfiles/DU_1stPaper_opSA.csv";
[sc_path, sc_name, sc_ext] = fileparts(scenarioPath);
sc_cell = readcell(scenarioPath);
sc = cell2table(sc_cell(2:end,:),"VariableNames",sc_cell(1,:));
ntests = height(sc);
% Load default input values
tableD = readtable('inputfiles/DefaultScenario.csv');

% loading sound for successfull completion
playSound = 0;
[y, Fs] = audioread('success_run_sound.mp3');

% create output directory if it doesn't already exist
if ~exist("results/"+sc_name, 'dir')
       mkdir("results/"+sc_name)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MINEGEOTHERMAL MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
Touts = zeros(10,ntests)*nan;

% Build the object in which all the outputs of function will be stored
Outputs = table();
Outputs = PropertyStruct(Outputs);

for irun = 1:ntests
%      clearvars -except atTout irun fid Fs ntests playSound sc sc_cell sc_ext sc_name sc_path scenarioName scenarioPath Touts verbose y tableD
    %% Assign Run Parameters
    % bundle parameters into propeties
    Outputs.('irun') = irun;
    [PhysicalProperties, NumericalProperties] = bundleProperties(sc, tableD, irun);
    vtkf = sc_name+"_"+char(string(NumericalProperties.tag)); % output VTK file name
    Outputs.('sc_name') = sc_name;
    Outputs.('vtkf') = vtkf;
    fprintf("\n\nXXXX Running "+vtkf+" XXXX\n");

    % Read in geometries and return arrays required by mineflow
    disp('Running: geometries.m - generating problem geometry')
    %Outputs = geometries(PhysicalProperties, NumericalProperties, Outputs);
    Outputs = call_step(PhysicalProperties.igeom, PhysicalProperties, NumericalProperties, Outputs);
    disp('Call complete: geometry generation sucessful, flow model variable initialised')
    disp('-------')

    % Random seed for the RNG
    rng(NumericalProperties.rng_seed);

    %% RUN
    % Calculate pipe lengths, initialise pipe diameters as a vector
    Outputs = mine_array_setup('1', PhysicalProperties, NumericalProperties, Outputs);

    %%% MATERIAL PROPERTIES
    % Assigns the material properties
    Outputs = material_properties(PhysicalProperties, NumericalProperties, Outputs);

    %%% INITIAL CONDITIONS
    % Apply rock temperatures
    Outputs = source_terms(5, PhysicalProperties, Outputs);

    %%% SOURCE TERMS
    % Apply Inflow Temperatures
    Outputs = source_terms(3, PhysicalProperties, Outputs);
    % Apply Groundwater Source Terms
    Outputs = source_terms(4, PhysicalProperties, Outputs);

    % Calculate flow through pipe system:
    disp('Running: mineflow.m - calculating hydraulic heads and pipe flow volumes')
    Outputs  = call_step(NumericalProperties.flow_model,PhysicalProperties, NumericalProperties, Outputs);
    disp('Call complete: geometry generation sucessful, flow model variable initialised')
    disp('-------')

    %%% Setup pipe flow arrays:
    Outputs = mine_array_setup('2', PhysicalProperties, NumericalProperties, Outputs);

    %%% Calculate temperature of pipe system:
    disp('Running: mine_heat.m - calculating nodal and pipe temperatures')
    Outputs = call_step(NumericalProperties.heat_model, PhysicalProperties, NumericalProperties, Outputs);
    fprintf('\nCall complete: thermal calculation sucessful')
    disp('-------')

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% POST-PROCESSING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Add useful outputs
    Outputs.('head') = [Outputs.H; Outputs.Ho];
    % convert Tn to the relative range (Tr-Tin)
    internal_nodes = internalState.nodes;
    Outputs.('T_frac') = (Outputs.Tn-min(Outputs.Tn(internal_nodes(PhysicalProperties.q_in))))/(max(Outputs.Tr_st)-min(Outputs.Tn(internal_nodes(PhysicalProperties.q_in))));
    % Calculate Reynolds number through the system:
    Outputs.('Re') = 4 * abs(Outputs.Q) ./ (pi*Outputs.d * PhysicalProperties.nu_f);
    % Outputs the GIS node numbers
    Outputs.('node_id') = internalState.nodes';
    % Outputs the pipe numbers
    Outputs.('pipe_id') = [1:Outputs.np]';
    Outputs.Tp = Outputs.Tp.';
    % Create output to show injection and abstraction points
    Outputs = assign_well_fhn_visualisation(NumericalProperties, Outputs);
    % Convert the Q magnitudes to a vector, for easy glyph plotting in
    % paraview
    Outputs.('Qv') = vectorise(Outputs.Q, Outputs.pipe_nodes, Outputs.xtotal, Outputs.L).';
    
    % Required to extract the data from the model based on the user input
    Outputs = build_outputs(PhysicalProperties, NumericalProperties, Outputs);

    % This step is required to merge porous zone and pipe zone data - For
    % now comment out if you do not wish to run a porous zone model
    % Outputs = step_recombine_model_data(PhysicalProperties, NumericalProperties, Outputs);
    
    %% Generate VTK file compatible with ParaVIEW
    % step_vtk_factory_mixed(PhysicalProperties, NumericalProperties, Outputs); %for porous flow
    step_vtk_factory(PhysicalProperties, NumericalProperties, Outputs); % for pipe only

    %Play 2 beeps on completion / 1 beep would indicate a failed run
    if playSound == 1
        try % in case the computer audio is innacessible or has been changed
        sound(y, Fs, 16);
        catch
            %%
        end
    end

    % store the outflow well temperatures for future saving in a single
    % file
    Touts(1:length(PhysicalProperties.q_out),irun) = Outputs.T_frac(internal_nodes(PhysicalProperties.q_out));


    % Save all outflow temperatures from the wells to a single file.
    % Moved here to avoid loosing data if one scenario crashes.
    writetable(table(Touts), "results/"+sc_name+"/"+sc_name+"_output.csv")

    fclose('all'); % makes sure all the open files are closed to release any memory
end
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% HELPER FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function v = vectorise(x, pipe_nodes, xtotal, L)
%   v = vectorise(x, pipe_nodes, xtotal, L) takes input vectors x, pipe_nodes,
%   xtotal, and L to compute the vector components and returns a matrix v.
%
% Input arguments:
%   x          - Vector containing values corresponding to each pipe node.
%   pipe_nodes - Matrix specifying pipe connectivity with two columns: 
%                [n_in, n_out], where n_in and n_out are node indices.
%   xtotal     - Matrix containing node coordinates in the form [x, y, z] 
%                (for 3D) or [x, y] (for 2D) for each node.
%   L          - Vector containing the length of each pipe corresponding 
%                to the nodes in 'pipe_nodes'.
%
% Output:
%   v          - Matrix of size np-by-3, where np is the number of pipes. 
%                The columns represent the vector components [vx, vy, vz].
%                Note: vz will be zero for 2D data (when size(xtotal, 2) == 2).
%
    np = size(pipe_nodes,1);
    v = zeros(length(x),min(size(xtotal)));
    for ip = 1:np
        n_in = pipe_nodes(ip,1);
        n_out = pipe_nodes(ip,2);
        x_in = xtotal(n_in,1);
        x_out = xtotal(n_out,1);
        y_in = xtotal(n_in,2);
        y_out = xtotal(n_out,2);
        ratioxL = x(ip)/L(ip);

        v(ip,1) = (x_out - x_in) * ratioxL ;% vx
        v(ip,2) = (y_out - y_in) * ratioxL ;% vy
        if min(size(xtotal)) == 3
            z_in = xtotal(n_in,3);
            z_out = xtotal(n_out,3);
            v(ip,3) = (z_out - z_in) * ratioxL ;% vz
        else
            v(ip,3) = 0;
        end
    end
end

function Outputs = assign_well_fhn_visualisation(NumericalProperties, Outputs)
% assign_well_fhn_visualisation - Assign well and FHN visualization values
%
%   Outputs = assign_well_fhn_visualisation(NumericalProperties, Outputs)
%   takes input PhysicalProperties, NumericalProperties, and Outputs, and 
%   assigns well and FHN (fixed-head node) visualization values to the 
%   'Outputs' structure.
%
% Input arguments:
%   NumericalProperties - Structure containing numerical properties data.
%                         It must have the following field:
%       .fhn  - string of 2D array representing the fixed head data (format: [node1, fixed_head_value1; node2, fixed_head_value2).
%   Outputs             - Structure containing output data that will be modified.
%                         It must have the field:
%       .N      - Total number of nodes.
%       .q_in  - Vector of node indices representing wells' inflow points.
%       .q_out - Vector of node indices representing wells' outflow points.
%
% Output:
%   Outputs - Modified structure containing the following additional field:
%       .wells - Vector of size N-by-1 representing the well and FHN 
%                visualization values for each node. Possible values:
%                1 - Represents a well inflow node.
%               -1 - Represents a well outflow node.
%                0 - Represents an FHN node.
%                NaN - Represents regular nodes without wells or FHN data.
%    
    nn = Outputs.nn;
    no = Outputs.no;
    q_in = Outputs.q_in;
    q_out = Outputs.q_out;
    internal_nodes = internalState.nodes;

    wells = nan(nn+no,1);
    for i = 1:length(q_in)
        wells(internal_nodes(q_in(i))) = 1;
    end
    for i = 1:length(q_out)
        wells(internal_nodes(q_out(i))) = -1;
    end
    fhns = str2num(NumericalProperties.fhn);
    for i = 1:size(fhns,1)
        wells(internal_nodes(fhns(i,1))) = 0;
    end
    Outputs.('wells') = wells;
end

function Outputs = call_step(stepName, PhysicalProperties, NumericalProperties, Outputs)
% call_step - Call a specified step function for data processing
%
%   Outputs = call_step(stepName, PhysicalProperties, NumericalProperties, Outputs)
%   takes input 'stepName' specifying the name of the step function to be
%   executed, and 'PhysicalProperties', 'NumericalProperties', and 'Outputs'
%   containing relevant data for the specified step. The function calls the 
%   specified step function and returns the updated 'Outputs' structure.
%
% Input arguments:
%   stepName           - Name of the step function to be executed. It should
%                        be a valid function name as a character array.
%   PhysicalProperties - Structure containing physical properties data.
%                        This data is passed as an input to the specified step.
%   NumericalProperties - Structure containing numerical properties data.
%                         This data is passed as an input to the specified step.
%   Outputs             - Structure containing output data. This data is 
%                         passed as an input to the specified step, and the
%                         function will return an updated version of it.
%
% Output:
%   Outputs - Updated structure containing the output data after the 
%             specified step function is executed.
%     
    if stepName == ""
        return
    end
    fh = str2func(stepName); % creates a function handle from a char array
    Outputs = fh(PhysicalProperties, NumericalProperties ,Outputs);
end