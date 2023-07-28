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
function [n_tree_idx, n_tree, np] = tree_branche(cnode, n_tree_idx, n_tree, np, node_pipes_out, node_pipes_in, npipes, pipe_nodes, Q, Qth)
    % This is a recursive function which creates a tree by connecting the
    % nodes that are connected together by flowing pipes. It allows to
    % reduce the heat computation by avoiding nested for loops in mineheat.m. It
    % could also be used to establish the mine network's connectivity by
    % checking which nodes are not in the tree using the n_tree_idx.
    verbose = 0;

    % Prints out VTK files to highlight the tree formation every time 5
    % nodes are added to it
    global xtotal;
    if verbose ==1 && mod(np-1,5) == 0
        inTree = n_tree_idx;
        inTree(inTree > 0) = 1; 
        vtk_factory("tree_formation_", np, np, pipe_nodes, xtotal, {inTree, npipes}, ["InTree", "Nconnections"], {}, {});
    end

    [cnode, n_tree_idx, n_tree, np, go_to_next] = tree_add(cnode, n_tree_idx, n_tree, np, npipes, node_pipes_in, pipe_nodes, Q, Qth);
    if verbose == 1
        fprintf('%d pipes leave node %d, and %d flow into it. Next-> %d\n', npipes(2,cnode), cnode, npipes(1,cnode),go_to_next);
    end
    if go_to_next == 1   
        % iterate through all the pipes coming out of cnode
        pipes = transpose(node_pipes_out(1:npipes(2,cnode), cnode));
        for i = 1:length(pipes)
            ip = pipes(i);
            nnode = pipe_nodes(ip,2); % the next node to consider
            if verbose == 1
                fprintf('next node is %d\n', nnode);
            end
            [n_tree_idx, n_tree, np] = tree_branche(nnode, n_tree_idx, n_tree, np, node_pipes_out, node_pipes_in, npipes, pipe_nodes, Q, Qth);
        end
    end
    % end of the tree branch reached, go back up to the next valid
                % fork in the tree
end