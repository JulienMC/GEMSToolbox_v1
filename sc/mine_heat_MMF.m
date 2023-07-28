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

%% DESCRIPTION
% This function computes the temperature distribution in a mine network 
% considering water and heat flow. It uses an analytical and semi-analytical 
% approach combining the planar and radial methods to account for thermal
% interference between galleries. The function considers fixed inflow 
% temperatures and can handle various thermal properties of rocks and fluids.
%
% INPUTS:
%   PhysicalProperties: Structure containing physical properties of the mine network.
%     - nyrs: The number of years of the simulation.
%     - k_f: Fluid thermal conductivity [W m^-1 K^-1].
%     - nu_f: Fluid kinematic viscosity [m^2 s^-1].
%     - rho_f: Fluid density [kg m^-3].
%     - Cp_f: Fluid heat capacity [J kg^-1 K^-1].
%     - eps: The roughness of the walls [m].
%     - k_r: Rock thermal conductivity [W m^-1 K^-1].
%     - Cp_r: Rock heat capacity [J kg^-1 K^-1].
%     - rho_r: Rock density [kg m^-3].
%
%   NumericalProperties: Structure containing numerical properties for the flow computation.
%     - int_d: The number of points used to evaluate the thermal interference distance around a pipe and determine alpha.
%     - max_sim_time: The maximum number of years after which the outflow temperature is assumed to be that of the rock temperature.
%
%   Outputs: Structure containing outputs from previous iterations (initially empty).
%     - d: Pipe diameters [m].
%     - L: Pipe lengths [m].
%     - Q: Pipe flow rates [m^3 s^-1].
%     - v: Pipe fluid velocity [m s^-1].
%     - np: The number of pipes.
%     - nn: The number of unknown head nodes.
%     - no: The number of known fixed head nodes.
%     - npipes: Matrix storing the number of pipes flowing in and out of each node size(2,nn).
%     - node_pipes_out: Matrix storing the ids of pipes flowing out of each node size(max(npipes),nn).
%     - node_pipes_in: Matrix storing the ids of pipes flowing into each node size(max(npipes),nn).
%     - pipe_nodes: Matrix storing the ID of the inflow and outflow node of each pipe size(nn,2).
%     - xtotal: Matrix storing the spacial coordinates of each node.
%     - Tr_st: Matrix storing the rock temperature at each node size(nn+no,1).
%     - T_st: Matrix storing the inflowing fluid temperature at each node size(nn+no,1).
%     - n_tree: A tree organizing all the nodes in the order in which their temperature needs solving.
%
% OUTPUTS:
%   Outputs: Updated structure containing the results of the flow calculation.
%     - Tn: The nodal temperatures [°C] size(nn+no,1).
%     - Tp: The pipe inflow and outflow temperatures. size(np,2).
%     - rp: The radial analytical thermal front distance of each pipe. size(np,1).
%     - tt_n: The transit time, the shortest time it takes for water to reach this node from any inflow point. size(nn+no,1).
%     - pmod: The rule applied to the temperature calculation.  
%           1: No special case - the main model is applied
%           2: Set to Tr because it take over the max_sim_time years for water to flow through the pipe
%           3: Set to Tr because inflow temperature equals Tr, so no heat
% exchange happens
%     - alphas: The weighting factor to account for the thermal interference. size(np.1).
%     - Tn_max: The maximum potential nodal temperature [°C]. size(nn+no,1).
%     - Tn_max: The minimum potential nodal temperature [°C]. size(nn+no,1).
%

function Outputs = mine_heat_MMF(PhysicalProperties, NumericalProperties, Outputs)
    % More output needed? --> Set verbose to 1.
    verbose = 0;

    % Loading Physical Properties for sub function
    nyrs = PhysicalProperties.nyrs;
    k_f = PhysicalProperties.k_f;
    nu_f = PhysicalProperties.nu_f;
    rho_f = PhysicalProperties.rho_f;
    Cp_f = PhysicalProperties.Cp_f;
    eps = PhysicalProperties.eps;
    k_r = PhysicalProperties.k_r;
    Cp_r = PhysicalProperties.Cp_r;
    rho_r = PhysicalProperties.rho_r;

    % Numerical Properties
    int_d = NumericalProperties.int_d;
    max_sim_time = NumericalProperties.max_sim_time; % years

    % Loading previous outputs
    d = Outputs.d;
    L = Outputs.L;
    Q = Outputs.Q;
    v = Outputs.v;
    np= Outputs.np;
    nn= Outputs.nn;
    no = Outputs.no;
    npipes = Outputs.npipes;
    node_pipes_out = Outputs.node_pipes_out;
    node_pipes_in  = Outputs.node_pipes_in;
    pipe_nodes = Outputs.pipe_nodes;
    xtotal = Outputs.xtotal;
    Tr_st = Outputs.Tr_st;
    Tr_st_max = max(Tr_st);
    T_st  = Outputs.T_st;
    n_tree = Outputs.n_tree;

    % local variables pre-allocation
    Tn = zeros(nn+no,1);
    Tn_max = zeros(nn+no,1);
    Tn_min = zeros(nn+no,1);
    Tp = zeros(np, 2);
    rp = nan(np, 1); %radius from which heat is extracted around the pipe
    alphas = zeros(np, 1); % the alpha parameter of each pipe
    % Impose inflow temperature on nodes that have external inflow
    % (which will be updated later with T from inflowing pipes if applicable)
    weighted_flow_in = zeros(np,1);

    % stores the type of handling performed on the pipe for debugging
    pmod = ones(np,1)*nan;

    % local variables assigment
    tree_length = length(n_tree(n_tree>0));
    progressStep = fix(tree_length/1000);
    reverseStr = ''; % For tracking percentage progress in command line

    % calculate the shortest transit time for a potential heat disturbance
    % to reach every node.
    tt_n = CalcTransitTimes(npipes,node_pipes_in,pipe_nodes,n_tree,L,v,Tr_st);
    % calculate the pipe center coordinates and appropriate data
    [pipe_centres, xs1, ys1, zs1, np] = CalcPipeCentres(pipe_nodes, xtotal);

    t = 3600*24*365*nyrs;
    r=d./2;
    Qmax = max(abs(Q));
    % First, set T of all incoming nodes:
    for inode = 1:nn+no
        sumQin = 0.;
        sumQout = 0.;

        Tn(inode) = Tr_st(inode); % Sets all nodes to the background rock temp

        % Flow balance
        external_inflow=0.;
        npipe_in = npipes(1,inode);
        % If node has inflow pipes, sum that inflow
        if npipe_in>0
            sumQin = sum(Q(node_pipes_in(1:npipe_in,inode)));
        end

        npipe_out = npipes(2,inode);
        % If node has outflow pipes, sum outflow
        if npipe_out>0
            sumQout = sum(Q(node_pipes_out(1:npipe_out,inode)));
        end
        % compute difference between in and outflow
        sumQ = sumQin - sumQout;
        if (sumQ/Qmax<-1e-6) % significant external inflow at node:
            external_inflow = -sumQ;
            Tn(inode)=T_st(inode)*external_inflow/(sumQin+external_inflow);  % Give all these nodes initial inflow T at first
        end

        % If more than 1 pipe flows into node, then T of node will be a weighted
        % average of the T of those inflowing pipes.
        % Those weights are stored in weighted_flow_in
        weighted_flow_in(node_pipes_in(1:npipe_in,inode)) = Q(node_pipes_in(1:npipe_in,inode))/(sumQin+external_inflow);
    end

    % Solves for the system temperature using the n_tree which ensures
    % that when every node is solved all the required significant inflow Temps have already
    % been calculated or are known boundaries.
    for i = 1:length(n_tree)
        in = n_tree(i);

        % checks if end of tree reached
        if in == 0
            break;
        end

        Tr = Tr_st(in); % set Tr using the specified value in the source term array

        % Because all the nodes are set to Tr initially, if a node is still
        % set to Tr (i.e. not overriden by inflow T), then set it to 0
        % because its full temperature will be determined in this iteration
        if Tn(in) == Tr
            Tn(in) = 0;
        end

        % iterates through the current node's incoming pipes
        pipes = transpose(node_pipes_in(1:npipes(1,in), in));
        if weighted_flow_in(pipes) > 1
            warning("At node %d, weighted flows are greater than 1.\n",in)
        end

        for j = 1:length(pipes)
            ip = pipes(j);
            fn = pipe_nodes(ip,1); % the inflow node at the other end of the pipe
            %  (start of pipe ip) fn ---ip---> in (end of pipe ip)
            Tp(ip,1) = Tn(fn); % this will work because any edge nodes will be at Tr
            tbt = L(ip)/v(ip); % breakthrough time for water to flow through the pipe.
            tt = tt_n(fn); % transit time to the node at the start of the pipe.
            tBt = tbt + tt;% the time for water to flow to the pipe, and through it.

            t_in = max([t-(tt+tbt),0]); % We offset the time for the computation of the solution to account for the time it takes for the closest temperature distrubance to reach the node

            %% A case to handle the computation of heat by combining the Cylindrical approach
            % from Rodriguez and Diaz, and the planar approach from
            % Pruess and bodvarsson, using a weighting factor
            if Tn(fn) == Tr || tBt > max_sim_time*(3600*24*365) || t_in == 0
                pmod(ip) = 2;
                % if the tBt is greater than max_sim_time the water has
                % not had yet time to flow through the pipe, so
                % Tout is set to Tr as breakthrough time is far
                % away
                Tp(ip,2) = Tr;
                rp(ip,1) = nan;

                if Tn(fn) == Tr
                    pmod(ip) = 3;
                end

                Ttemp = ones(2,1)*Tr;
            else
                % determine all the radii related to the pipe
                [~, rii_mode_frac, rii_minmax, r0] = calcRadii(v,t_in,r,L,ip,pipe_nodes,xtotal,d, int_d, k_r, rho_r, Cp_r, k_f,nu_f,rho_f,Cp_f,eps, pipe_centres, xs1, ys1, zs1, np);
                Ttemp = zeros(2,1);
                for ir = 1: 2
                    if rii_minmax(ir) >= r0 % then we use R&D
                        [Ttemp(ir,1), ~ ]= pipeheat(r(ip), L(ip), Tn(fn), v(ip), t_in, Tr, k_r, rho_r, Cp_r, k_f, nu_f,rho_f, Cp_f, eps);
                    else % if not we are radius limited and need to use P&B
                        pmod(ip) = 1;
                        ts = heatplate_PB(Tp(ip,1), Tr, L(ip), 2*pi*r(ip), Q(ip)*rho_f, Cp_f, k_r, rho_r, Cp_r, t_in+tbt, pi*r(ip)^2, rho_f);
                        Ttemp(ir,1) = ts(end,end);
                    end
                end
                
                % Determine the outflow temperature as the weighted
                % average of the temperature obtained for the min
                % and max critical radii found.
                % Because max(rii_c) is now always == r0, for ir ==
                % 2, we will always use R&D. For ir it will depend
                % on whether interference happens.
                Tp(ip,2) = Ttemp(2,1)*rii_mode_frac + Ttemp(1,1)*(1 - rii_mode_frac);
                alphas(ip,1) = rii_mode_frac;
                rp(ip,1) = r0;
            end
            
            if (verbose)
                fprintf('Solved T for pipe %d: T= %f -> %f\n',ip,Tp(ip,1),Tp(ip,2))
            end

            % Clamps calculated Temperatures to physically possible Temp
            Tp(ip,2) = max(min(Tr,Tn(fn)),min(max(Tr,Tn(fn)),Tp(ip,2)));

            % T at end of pipe contributes to T at node at end of pipe:
            Tn(in) = Tn(in) + weighted_flow_in(ip) * Tp(ip,2);
            Tn_max(in) = Tn_max(in) + weighted_flow_in(ip) * Ttemp(1,1);
            Tn_min(in) = Tn_min(in) + weighted_flow_in(ip) * Ttemp(2,1);

            if Tn(in) > max([T_st(in); Tr_st_max; Tn(fn)])*1.001
                fprintf('WARNING Node %d Temperature = %f C, which is greater than maximum possible model temperature of %f C\n',in, Tn(in), max([T_st; Tr]));
            end
        end

        %display % progress
        if mod(i,progressStep) == 0
            percentDone = i/tree_length* 100;
            msg = sprintf('Percent done: %3.1f', percentDone); %Don't forget this semicolon
            fprintf([reverseStr, msg]);
            reverseStr = repmat(sprintf('\b'), 1, length(msg));
        end
    end

    % Bundle outputs in struct
    Outputs.('Tn') = Tn;
    Outputs.('Tp') = Tp;
    Outputs.('rp') = rp;
    Outputs.('tt_n') = tt_n;
    Outputs.('pmod') = pmod;
    Outputs.('alphas') = alphas;
    Outputs.('Tn_max') = Tn_max;
    Outputs.('Tn_min') = Tn_min;
end

