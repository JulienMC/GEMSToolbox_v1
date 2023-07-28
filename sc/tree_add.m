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

function [cnode, n_tree_idx, n_tree, tp, go_to_next] = tree_add(cnode, n_tree_idx, n_tree, tp, npipes, node_pipes_in, pipe_nodes, Q, Qth)
% This functions adds the node to the tree if all the upstream nodes have
% been added to the tree already or if the flow from an incoming pipe is so
% negligeable that it can be ignored. Starting nodes are not added to the
% tree, as their temperature is known, but the go_to_next flag is set to true
% to keep recursing from the starting node.
%             cnode = the current node that is being considered
%             % Array to keep track of node index in tree
%             n_tree_idx = zeros(1,nn+no);
%             % Array representing the flow-tree in which to solve for Temp
%             n_tree = zeros(1,nn+no); % note that the end of the tree might
%             be 0 pgo_to_next if not all the nodes are connected to the mine
%             network.
%             % tp is the index of the tree's next vacant position
%             go_to_next = 1 or 0, to determine whether to keep recursion
%             going in branch_add.m

    verbose = 0;
    % Have we tried adding the node to the tree an amount of
    % time equal to the number of incoming pipes (i.e. is this
    % is now the last time chance to add this node and all
    % incoming nodes have been go_to_next).
    go_to_next = 0;

    % Determines whether the upstream nodes have been added to the tree
    % and/or if incoming pipes can be ignored as very low contribution to
    % the inflow to cnode
    count = 0;
    if npipes(1, cnode) > 0
        incoming_pipes = node_pipes_in(1:npipes(1,cnode),cnode);
        for i = 1:length(incoming_pipes)
            p = incoming_pipes(i);
            % fetch previous node from the incoming pipe
            n = pipe_nodes(p, 1);
            % check if the incoming node is already in the tree, or a
            % starting point
            if n_tree_idx(n) > 0 || npipes(1,n) == 0
                count = count + 1;
            % or if its flow is negligeable (less than 1e-10)
            % compared to the maximum model flow
            elseif Q(p)/max(Q) < Qth || Q(p) / sum(Q(incoming_pipes)) < 1e-3
                count = count + 1;
            end
        end
    end    

    if npipes(1,cnode) == 0
        go_to_next = 1; % here we set the flag but don't actually add the node to the tree as it is a starting node, hence its T in known

    elseif (count == npipes(1,cnode)) && n_tree_idx(cnode) == 0
        % final time the node will be encountered or other inflow is negligeable.
        n_tree(tp) = cnode; % adds the node to the tree
        tp = tp+1; % moves the tree position marker to the next free position in the tree
        n_tree_idx(cnode) = tp; % marks the node pos in the tree 
        go_to_next = 1; % to inform the recursive algo that the node has been go_to_next
    end
end
