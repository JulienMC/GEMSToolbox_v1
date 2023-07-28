# Installation
Install MATLAB R2021b or above
Install ParaView v5.10.1
Download the GEMSToolbox archive.

To install the code you need to unzip the archived into a directory.
The location of the `GEMSToolbox.m` file is considered the home directory.

# How to run the code?
Open MATLAB and make sure you are located in the home directory.

To run the code, run the `GEMSToolbox.m file`

You can find some input files in the `inputfiles` directory

To use one of these files simply write their name in the provided  `_InputFileName.csv` located in the `inputefiles` directory.

## View outputs
The outputs of each Scenario input file are saved in the `results/[scenario_filename]` folder.
The outputs consit of:
1. A .vtk file for each of the runs in that scenario named: `[scenario_filename]_[tag]_[runindex]``
2. A .csv file which records the relative temperature of each outflow node in each run. The relative temperature is given by [temperature_at_outflow_node]/abs(max(rock_temperature)-min(water_temperature))

## Making your own input file
Here are the properties which can be set in the input file. Any property not setup in the input file will be read from the 'DefaultScenario.csv' file.
**/!\ Only edit 'DefaultScenario.csv' if you know what you are doing.**

Each line of an input file represents a distinct run.

1. Copy the DefaultScenario.csv, or create a blank csv file.
2. Edit the following parameters as you see fit:
- **Tag**: Allows you to specify a tag that will be added to the output .vtk filename.
- **node_outputs**: Specifies the desired nodal outputs to be saved in the vtk file for each run.
    List of valid node outputs:
    - node_id: the id of the node found in the underlying GIS files or/and ID to be used to set inflow, outflow, or fixed head nodes.
    - Tn: the nodal temperatures of the model.
    - head: the nodal hydraulic heads of the model.
    - tt_n: the time it takes water to flow from the closest net inflow point to this node.
    - wells: indicates the inflow (1), outflow (-1) and fixed head (0) nodes. All other nodes have a NaN value.
    - n_tree: the position of nodes in the tree used to sequentially solve the heat transfer.
- **pipe_outputs**: Specifies the desired pipe outputs to be saved as cell values in the vtk file for each run.
    List of valid pipe outputs:
    - pipe_id: the id of the pipes as held within the GEMS toolbox.
    - Re: the Reynolds number.
    - v: the water velocity inside the pipe (m/s) assumed uniform.
    - Qv: the vectorised flow rate through the pipes (m^3/s). Saved as a vector in VTK with an x, y and z component.
    - Q: the absolute flow rate through the pipes. Is equal to the magnitude of Qv.
    - poro: the porosity of each pipe.
    - rp: the radial distance of the thermal front inside an infinite homogeneous porous rock.
    - Tp: a vectorised value holding two components labelled x and y. x corresponding to the inflow temperature of the pipe and y to its outflow temperature.
    - pmod: indicates which case of the if statement is used in the heat calculation code.
    - Qdiff: the flow difference in the pipe between the last and penultimate iteration in the flow computation.
- **igeom**: Specifies the geometrical function file to be used by the code. Note the `.m` extension is not included.
    List of valid geometrical functions:
    - step_load_GIS_geometry: will read the `ArcGIS_file_name` input shape file and will converted it into a GEMSToolbox compatible geometry.
    - step_grid_geometry: will generate a set of uniform grid seams, as per the specified `grid_width`, `grid_height`, `no_seams`, `h_pipe_length`, `v_pipe_length`, `seam_spacing` and `connections` input keywords.
- **ArcGIS_file_name**: the path to the shapefiles to be used to generate the geometry with `step_load_GIS_geometry`. The format is `[[directoryPath]/[fileName1].shp, [directoryPath]/[fileNameN].shp]`.
- **nyrs**: The number of years the simulaton is to be run for.
- **qset**: The flow rate of the injection and abstraction wells in m^3.s^-1.
    qset can be set as follows:
        float value -> will be applied to all the wells. format: `float`. e.g. `0.0186`
        array of foats -> will be applied to the individual wells. format: `[float float ... float]`. The values are applied to the injection wells starting from the first value, and then to the output wells. e.g. 0.04 0.01 0.02 0.03 for 2 injection well and 2 abstraction wells. Note that the total sum of injected volumes is equal to 0 (= (0.04 + 0.01) * (-1 for injection) + (0.02 + 0.03) * (1 for abstraction)).
- **q_in**: The injection nodes. Specified as integers not exceeding the total number of nodes in the selected geometry. format: `[int int ... int]`.
- **q_out**: the abstraction nodes. Specified as integers. (Note the total number of points specified in q_in and q_out should match the number of values in qset, if the array input is used). format: `[int int ... int]`.
- **k_r**: thermal conductivity of the rockmass surrounding the mine. W.m^−1.K^−1
- **Cp_r**: specific heat capacity of the rockmass surrounding the mine. J.kg^-1.K^-1
- **rho_r**: density of the rockmass surrounding the mine. kg.m^-3
- **k_f**: thermal conductivity of the water inside the mine.
- **Cp_f**: specific heat capacity of the water inside the mine.
- **rho_f**: density of the mine water. kg.m^-3
- **nu_f**: dynamic viscosity of the water in the mine. Pa.s
- **Tf_ini**: Injection temperature of the water (C). Array of floats. If a single value is specified it will be applied to all injection wells. Otherwise a value for each well is requried.
- **Tr**: Geothermal gradient used to compute the starting rock temperature. The first value is the gradient in C.m^-1, the second value is the temperature at 0 meters in C. format `[float float]`. e.g. [0.0376 17].
- **d_set**: The diameter of the mine galleries (m). Constant value for all the mine.
- **eps**: The absolute roughness of the mine galleries (m).
- **heat_model**: The heat model function file name to be used. Can select between: `mine_heat_radial` for the Rodriguez and Diaz (2009) reworked model.  `mine_heat_MMF` for a model which combines R&D with a planar heat transfer model to account for thermal interference between galleries, fast and RECOMMENDED. `mine_heat_FDprofile` for a model computing heat in the most conceptually representative way (most computationally expensive). `mine_heat_FDreload` allows the user to run consecutive runs whilst retaining the rockmass thermal state (computationally expensive).
- **flow_model**: The flow model function file name to be used. `step_mineflow_pipes` to be used for models composed of pipe networks only. `step_mineflow_mixed` to be used for models composed of a mixture of pipe networks and porous zones (in dev).
- **Qth_tree**: the cutoff value for the flow comparison check used when disregarding a node for the computation of the flow tree required for the sequential heat calculations. This might need adjusting if when adding porous pipes sections of the models are not added to the tree.
- **int_d**: Number of integration point for the computation of the weighting factor for the `MixedModel` - 100 varies by about 3%, 1000 by less than 1%. Significantly affect the computational time of the `mine_heat_FDprofile` model.
- **B**: The exponent from Todini and Pilati (1988) for the flow computation.
- **sumdQrel_th**: The threshold against which the average change in flow between iterations is compared to determine convergence.
- **maxdQrel_th**: The threshold against which the maximum change in flow between successive iterations is compared to determine convergence.
- **a**: The damping oscillation factor applied to the flow at iteration n-1.
- **b**: The damping oscillation factor applied to the flow at iteration n.
- **num_doubles**: The number of consecutive instances for which the average change in flow hasn't change to 5 significant figures. Used to exist stagnating run.
- maxtraintime: Not yet implemented
- AI_modelpath: Not yet implemented
- **timesteps**: The minimum number of timesteps used per year to compute the FD schemes.
- **max_sim_time**: (in years) the maximum time it should take for water to reach the end of any given pipe. Beyond which the outflow temperature of a pipe will be considered to be that of the initial rock temperature.
- **n_steps**: The target amount of steps that the water should take to flow through any gallery in the model.
- **nz**: The amount of steps the rockmass should be discretised by in the FD scheme in the radial direction.
- **nitermax**: The maximum number of iterations used to determine the radius of the thermal front in the R&D method.
- **niter**: The minimum number of iterations used to determine the radius of the thermal front in the R&D method.
- VF_speed_th: The fluid velocity threshold below which the R&D model will apply the rock temperature as the outflow temperature. (m.s^-1) (Not yet implemented)
- **grid_width**: The number of nodes in the x dimension of the `step_grid_geometry` model (integer)
- **grid_height**: The number of nodes in the y dimension of the `step_grid_geometry` model (integer)
- **no_seams**: The number of seams to be created in the the `step_grid_geometry` model (integer)
- **h_pipe_length**: The length of individual pipes in the x direction for the `step_grid_geometry` model (meters)
- **v_pipe_length**: The length of individual pipes in the y direction for the `step_grid_geometry` model (meters)
- **seam_spacing**: The distance between seams for the `step_grid_geometry` model (meters)
- **connections**: The connections to impose in the `step_grid_geometry` model. Given as a 2D array of node pairs e.g `[1 43; 5 10]` would connect nodes 1 to 43 with a pipe, and 5 to 10 with a pipe.
- **fixed_head_node**: The nodes at which a specified fixed hydraulic head is applied. **These node should be different from the injection and abstraction nodes.**. An array is provided containing the node and the fixed head value for each node. format `[int float; int float; ...; int float]`. e.g. [10 0.0] fixing node 10 at a value of 0 meters.
- testbank: disused
- **reloadTrock**: Set to 0 to use the initial rock temperature Tr, set to 1 to use the rockmass temperature values form the previous run.
- **init_time**: Set to 0 when reloadTrock is 0, set to the amount of years corresponding to the time at which the rockmass temperatures were saved. (in years)
- **gw_darcy_velocity**: the maximum darcy velocity at which ground water flows from the surrounding rockmass into the mine. format: `float` (in dev)
- **mat_prop_case**: the switch case code used to apply the material properties to the model (currently porosity and hydraulic conductivity). 1: sets porosity to 1 for all pipes (i.e. all pipes are considered open). 2: sets all the pipes to the porosity and hydraulic conductivity specified in the `porosity` and `Hydraulic_conductivity` input keywords respectively. 3: uses a Perlin Noise map to set the porosity values of the pipes, and beyond a certain threshold (specified in `open_th` keyword) all the pipes are considered open with a porosity of 1. The value specified in the `Hydraulic_conductivity` input keyword is used as the hydraulic conductivity in all pipes. 4: like 3, but around a radius specified by the `well_exclusion_th` input keyword, where all pipes are considered open. 5: like 4 but instead of forcing open pipes around the wells, we force porous pipes.
- **Hydraulic_conductivity**: the value of the hydraulic conductivity to be used in porous pipes.
- **Porosity**: the porosity value of all pipes to be used for the `mat_prop_case` number 2.
- **rng_seed**: the random seed value used for all randomly generated data in the model. This is used for repeatability. format: `int`.
- **mat_d_prop**: An array of values `[float ... float]` indicating the new diameter (in meters) to apply to the pipes specified under the keyword `mat_pipe_ids`.
- **mat_K_prop**: An array of values `[float ... float]` indicating the new Hydraulic conductivity (meters/second) to apply to the pipes specified under the keyword `mat_pipe_ids`. A value of `0` will be ignored and the pipe forced to open.
- **mat_pipe_ids**: An array of values `[int ... int]` indicating the ID of the pipes to that will have their Hydraulic conductivity and diameter overwritten. Note that the following parameters should therefore all have the same number of values: `mat_pipe_ids`, `mat_d_prop` and `mat_K_prop`.
- **depths**: An array of values `[float ... float]` indicating the depths (in meters) to apply to any GIS shapefile specified in the `ArcGIS_file_name` keyword. Note that if a `POINT_Z` is specified in the shapefile then the matching value specified in this keyword will be ignored. For example if the user specifies 3 shapefile names, but only the first and last one have been assigned `POINT_Z` depth values in GIS. Then the user needs to specify three values in this keyword (e.g. `[50 100 300]`) but only the 100 m depth will be assigned as it corresponds to the second file specified, the only file that doesn't already specify a depth. Effectively, the 50 and 300 m depth values will be ignored and the values specified by the shapefile `POINT_Z` attribute will be used instead.
