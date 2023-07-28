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
%%%%% General puprpose book-keeping script called multiple times by
%%%%% minegeothermal.m

%%% Version 
%%%     20220112 - Function modified to accept a variable number of input
%%%     arguments, callCase added to determine which part of function is
%%%     actually running. Since mine_array_set up is a generic
%%%     'house-keeping' script, we use it to initialise major book-keeping
%%%     arrays and excecute ad-hoc calculations.

function Outputs = mine_array_setup(callCase,varargin)

%%% Extra output information?
verbose=0;

%%% Set up index arrays
switch callCase
    
    case '1'
        %%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% Define varargin %%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%
        [PhysicalProperties, NumericalProperties, Outputs] = varargin{1:end};
        
        d_set = PhysicalProperties.d_set;
        np = Outputs.np; 
        A12 = Outputs.A12;
        A10 = Outputs.A10;
        xo = Outputs.xo;
        x = Outputs.x;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% Begin computations %%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Array with start & end node of each pipe:
        pipe_nodes  = zeros(np,2);  % array to store start & end node of each pipe
        % note the 'global' node numbering, including
        %    fixed and free head nodes
        A102 = [A12 A10];           % merges table of fixed (A10) and unknown (A12) nodes
        [row,col] = find(A102==-1); % node with -1 indicates one end point of pipe
        pipe_nodes(row,1) = col;    % store that node nr in pipe_nodes array
        [row,col] = find(A102==1);  % node with 1 indicates the other end point
        pipe_nodes(row,2) = col;    % store that node nr in pipe_nodes array

        % Pipe lengths:
        switch size(x,2)
            case 2 % 2D case
                
                L   = zeros(np,1);
                xtotal = [x;xo];
                p1 = pipe_nodes(:, 1);
                p2 = pipe_nodes(:, 2);
                x1 = xtotal(p1,:);
                x2 = xtotal(p2,:);
                dx = x2-x1;
                L  = sqrt(dx(:,1).^2+dx(:,2).^2);
                
            case 3 % 3D case
                
                L = zeros(np,1);
                xtotal = [x;xo];
                p1 = pipe_nodes(:, 1);
                p2 = pipe_nodes(:, 2);
                x1 = xtotal(p1,:);
                x2 = xtotal(p2,:);
                dx = x2-x1;
                L  = sqrt(dx(:,1).^2+dx(:,2).^2+dx(:,3).^2);
                
                
        end
        
        d = d_set*ones(np,1);
        %%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Assign varargout %%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%
        
        Outputs.('pipe_nodes') = pipe_nodes;
        Outputs.('xtotal') = xtotal;
        Outputs.('L') = L;
        Outputs.('d') = d;
        
    case '2'
        %%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% Define varargin %%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%
       [PhysicalProperties, NumericalProperties, Outputs] = varargin{1:end};

        np = Outputs.np;
        nn = Outputs.nn;
        no = Outputs.no;
        A10 = Outputs.A10;
        A12 = Outputs.A12;
        pipe_nodes = Outputs.pipe_nodes;
        q = Outputs.q;
        Q = Outputs.Q;
        d = Outputs.d;
        Qth = NumericalProperties.Qth_tree;
                     
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% Begin computations %%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Swap indices if flow is from end to start node:
        flow_dir = Q<=0;
        for ipipe=1:np
            if flow_dir(ipipe)==1
                dummy               = pipe_nodes(ipipe,1);
                pipe_nodes(ipipe,1) = pipe_nodes(ipipe,2);
                pipe_nodes(ipipe,2) = dummy;
                Q(ipipe)            = -Q(ipipe);
            end
        end
       
       % Compute pipe flow velocity from volumetric flow rate
       v = Q ./ (pi*(d/2).^2);
        
        
        if (verbose)
            disp('Flow direction of pipes:')
            for ipipe=1:np
                fprintf('  pipe %d has Q = %f from pipe %d to %d\n', ipipe, Q(ipipe), pipe_nodes(ipipe,1), pipe_nodes(ipipe,2))
            end
        end
        
        % Store nr of incoming and outgoing pipes for every node;
        npipes = zeros(2,nn+no); % nr of incoming (first row) and outgoing (2nd row) pipes
        for ipipe=1:np
            % pipe ipipe enters node inode:
            inode = pipe_nodes (ipipe,2);
            % So update nr of p
            npipes(1,inode) = npipes(1, inode)+1;
            % pipe ipipe exits node inode:
            inode = pipe_nodes (ipipe,1);
            npipes(2,inode) = npipes(2, inode)+1;
        end
        
        %%% INITIALIZE
        % Create arrays to list all pipes entering/exiting every node:
        maxnpipes_in   =(max(npipes(1,:)));
        node_pipes_in  = zeros(maxnpipes_in+1,nn+no); % add 1 to maxnpipes_in to allow for potential external inflow (stored in array q)
        maxnpipes_out  =(max(npipes(2,:)));
        node_pipes_out = zeros(maxnpipes_out,nn+no);
        
        % Populate 2-D array of all pipes exiting every node:
        for ipipe=1:np
            % node inode at start of pipe ipipe:
            inode = pipe_nodes (ipipe,1);
            % add this pipe to the table of all pipes exiting node inode:
            % pipe not yet added to node_pipes_out array:
            pipe_reported = 0;
            % pipe to be added in one of the reserved places in array:
            for ipipe2=1:maxnpipes_out
                if pipe_reported == 0
                    if node_pipes_out(ipipe2,inode) == 0
                        % location in array still free, so enter pipe nr here:
                        node_pipes_out(ipipe2,inode) = ipipe;
                        % and flag pipe as entered
                        pipe_reported = 1;
                        % leave the ip loop
                        break;
                    end
                end
            end
        end
        
        % Create 2-D array of all pipes entering every node:
        for ipipe=1:np
            % node at end of pipe ipipe:
            inode = pipe_nodes (ipipe,2);
            % add this pipe to the table of all pipes entering node in:
            % pipe not yet added to node_pipes_in array:
            pipe_reported = 0;
            % pipe to be added in one of the reserved places in array:
            for ipipe2=1:maxnpipes_in
                if pipe_reported == 0
                    if node_pipes_in(ipipe2,inode) == 0
                        % location in array still free, so enter pipe nr here:
                        node_pipes_in(ipipe2,inode) = ipipe;
                        % and flag pipe as entered
                        pipe_reported = 1;
                        % leave the ip loop
                        break;
                    end
                end
            end
        end

        % JMC BUILDING NODE TREE FOR HEAT COMPUTATION
        % Array to keep track of node index in tree
        n_tree_idx = zeros(nn+no,1);
        % Array representing the flow-tree in which to solve for Temp
        n_tree = zeros(nn+no,1); % note that the end of the tree might
        % be 0 padded if not all the nodes are connected to the mine
        % network.
        np = 1; % the current node tree position free to add a new node
        for inode = 1:nn+no
            if npipes(1,inode) == 0 && npipes(2,inode) > 0% external inflow into (free-pressure) node
                % will create the tree branche starting from cnode   
                % we pass -1 as the previous pipe argument as there are none
                % fprintf("Starting node: %d\n",inode);
                [n_tree_idx, n_tree, np] = tree_branche(inode, n_tree_idx, n_tree, np, node_pipes_out, node_pipes_in, npipes, pipe_nodes, Q, Qth); 
            end
        end
        
        % Some nodes have external inflow, which needs recording too:
        for inode = 1:(nn)
            if q(inode)<0    % external inflow into (free-pressure) node
                % add this info to the node_pipes_in array
                external_inflow_reported = 0;
                % external inflow to be added in one of the reserved places in array
                % find first 'empty' location, and add a -1 there (to indicate external inflow):
                for index=1:maxnpipes_in
                    if external_inflow_reported == 0
                        if node_pipes_in(index,inode) == 0
                            % location in array still free, so enter a -1 here:
                            node_pipes_in(index,inode) = -1;
                            % and flag pipe as entered
                            external_inflow_reported = 1;
                            % leave the index loop
                            break;
                        end
                    end
                end
            end
        end    

            
        %%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Assign varargout %%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%
        Outputs.('pipe_nodes') = pipe_nodes;
        Outputs.('npipes') = npipes;
        Outputs.('node_pipes_out') = node_pipes_out;
        Outputs.('node_pipes_in') = node_pipes_in;
        Outputs.('Q') = Q;
        Outputs.('v') = v;
        Outputs.('n_tree') = n_tree;
        Outputs.('n_tree_idx') = n_tree_idx;
        %%%%%%%%%%%%%%%%%%%%%%%%%
end



