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

% a Function which reads shape files, build the geometry and joins various 
% geometries using the specified connextions (cnx).
% cnx is a specified as a user input as 2D array [[x y],[v w]] will connect node x to y and v to w
% it will be parsed as a 1D array [x y v w], and it will be reshaped back
% to a 2D array here. This is because it needs to be stored in the Table of
% PhysicalProperties.
function [A120, xtotal, N, np] = build_arc_geometry(filenames, cnx, depths)
ncnx = size(cnx,1); % number of connections
cnx = reshape(cnx,[ncnx 2]);
depths = str2num(depths);

T = {};
npoints = {};
flow_dir = {};
nps = {}; % the total number of pipes without the connections
Ns = {}; % the maximum node ID

% Removes parenthesis from input string
% and splits the string at ','  or ';'
filenames = regexprep(filenames,'[[](){} ]','');
filenames = regexp(filenames, '[,;]', 'split');
filenames = filenames{1};

for i = 1:length(filenames)
    [T{i}, npoints{i}, flow_dir{i}, nps{i}, Ns{i}] = build_file_geometry(filenames{i});
end


% size of problem:
no  = 1;       % nr of fixed head nodes

% Parameters to be solved in this function:
np = sum(cell2mat(nps)); % file pipes + the connections
N = sum(cell2mat(Ns));
A120 = zeros(np,N);
xtotal   = zeros(N,3);

% locations of nodes
prev_pts = 0;
prev_ps  = 0;
for i = 1:length(filenames)
    T1 = T{i};
    hasZ = ismember('POINT_Z',fieldnames(T1(1)));
    for irow = 1:npoints{i}
        % locate points
        nodeid = prev_pts + T1(irow).NODEID;
        if hasZ == 1
            xtotal(nodeid,:) = [T1(irow).POINT_X,T1(irow).POINT_Y,T1(irow).POINT_Z];
        else
            xtotal(nodeid,:) = [T1(irow).POINT_X,T1(irow).POINT_Y,depths(i)];
        end
        % locate pipes
        pipeid = prev_ps  + T1(irow).PIPEID;
        fd = flow_dir{i};
        A120(pipeid,nodeid) = fd(irow);
    end
    prev_pts = prev_pts +  Ns{i}; % adds the points from the current file to the count
    prev_ps  = prev_ps  + nps{i}; % adds the pipe count from the current file
end

% In linear pipe configuration, for each node in, pipe in feeds into 
% node in, and pipe in+1 leaves node in:

% % handles the connections
% for i = 1:ncnx
%     A120(np-ncnx+i,cnx(i,1))=1;
%     A120(np-ncnx+i,cnx(i,2))=-1;
% end

% clean the A120 matrix for errors
[A120, xtotal] = cleanNetworkMatrix(A120, xtotal);
%update the number of nodes and pipes
[np, N] = size(A120);

% handles the connections
for i = 1:ncnx
    A120(np+i,:) = 0;
    A120(np+i,cnx(i,1))=1;
    A120(np+i,cnx(i,2))=-1;
end

% update the number of nodes and pipes
[np, N] = size(A120);
A120 = sparse(A120);
end

function [T1, nshprows, flow_dir1, np1, N1] = build_file_geometry(filename)
    T1 = shaperead(filename);
    sizeofT1 = size(T1);
    nshprows = sizeofT1(1); 
    flow_dir1 = zeros(nshprows,1);
    
    for irow = 1:nshprows
        if (T1(irow).POINT_X == T1(irow).EXT_MIN_X && T1(irow).POINT_Y == T1(irow).EXT_MIN_Y)
           flow_dir1(irow) = -1; 
        end
        if (T1(irow).POINT_X == T1(irow).EXT_MAX_X && T1(irow).POINT_Y == T1(irow).EXT_MAX_Y)
           flow_dir1(irow) = 1;   
        end
    end

    % Counts the number of nodes
    max_node_id1 = 0;
    for irow = 1:nshprows
        if T1(irow).NODEID > max_node_id1
            max_node_id1 = T1(irow).NODEID;
        end
    end
    
    N1 = max_node_id1;
    
    % Counts the number of pipes
    max_pipe_id1 = 0;
    for irow = 1:nshprows
        if T1(irow).PIPEID > max_pipe_id1
            max_pipe_id1 = T1(irow).PIPEID;
        end
    end

    np1 = max_pipe_id1;
end