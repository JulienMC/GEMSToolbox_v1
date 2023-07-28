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

function [A120, xtotal, N, np] = geometry_multigrid(n,m,s,l1,l2,h,cnx)

% with: 
%   - = pip
%   o = unknown head node
%   * = fixed head node
% 
% size of problem:
% n, m, l1 & l2 are imported from ensemble run
% n      % grid width (number of nodes wide)
% m      % grid height (number of nodes high)
% l1     % length of horizontal pipes
% l2     % length of vertical pipes
% s      % the number of seams
% h      % the spacing of each seam
% cnx    % additional connections between two nodes in two different seams % of dimension (cnx_number, 2)  

N  = n*m*s;     % nr of unknown head nodes
np  = ((n-1)*m + (m-1)*n)*s;     % nr of pipes
A120 = sparse(np,N); % global incidence matrix

%% Assigning Node Coordinates
xtotal = zeros(N,3);
for in = 1:N
    so = n*m; % seam offset
    cs = floor((in-1)/(n*m)); %current seam
    xtotal(in,:) = [((in-cs*so-1)-floor((in-cs*so-1)/n)*n)*l1, floor((in-cs*so-1)/n)*l2, cs*h];
end

%% locations of boundary unknown-nodes:
for is = 1:s
    A12 = seam_incidence_mat(n,m,np/s);
    A120(1+(np/s)*(is-1):(np/s)*(is),1+(N/s)*(is-1):(N/s)*(is)) = A12;
end

%% Applying the connections between seams only
cnx_size = size(cnx);
ACnx = sparse(cnx_size(1),N);
for ic = 1:cnx_size(1)
    ACnx(ic,cnx(ic,:)) = [-1, 1]; % assigns outflow to the first node in cnx, and inflow to the second node
end
A120 = [A120; ACnx]; % adds the new connecting pipes to the A120 matrix
np = np + cnx_size(1); % adds the number of new pipes to the total pipe number
end

function A12 = seam_incidence_mat(n,m,np)
    % This function generates the incidence matrix for each seam
    % n is the number of nodes in the x direction for the seam
    % m is the number of nodes in the y direction for the seam
    % np is the number of pipes in the seam
    A12 = sparse(np,n*m);
    l1 = 10;
    l2 = 10;
    for in = 1:n*m
        if in==1 % corner1 - check
            A12(1,in) = -1;
            A12((n-1)*m+1,in) = -1;
        elseif in==n % corner2 - check
            A12((n-2)*m+1,in) = 1;
            A12((n-1)*m+n,in) = -1;
        elseif in==(m-1)*n+1 % corner3 - check
            A12((n-1)*m+(m-2)*n+1,in) = 1;
            A12(m,in) = -1;
        elseif in==n*m % corner4 - check
            A12((n-1)*m,in) = 1;
            A12((n-1)*m+(m-1)*n,in) = 1;
        elseif in>=2 && in<=n-1 % bottom boundary of grid - check
            A12((in-2)*m+1,in) = 1;
            A12((n-1)*m+in,in) = -1;
            A12((in-1)*m+1,in) = -1;
        elseif (in/n)>1 && (in/n)<m-1 && mod(in-1,n)==0 % left boundary of grid - check
            A12(in-1+(n-1)*(m-1),in) = 1;
            A12((in-1)/n+1,in) = -1;
            A12(in-1+n+(n-1)*(m-1),in) = -1;
        elseif (in/n)>1 && (in/n)<=m-1 && mod(in,n)==0% right boundary of grid - check
            A12(in-n+(n-1)*m,in) = 1; 
            A12(in/n+(n-2)*m,in) = 1;
            A12(in+(n-1)*m,in) = -1;
        elseif in>(m-1)*n+1 && in<m*n % top boundary of grid - 
            A12((in-(m-1)*n-1)*m,in) = 1;
            A12((n-1)*m+in-n,in) = 1;
            A12((in-(m-1)*n)*m,in) = -1;
        else % the middle of grid -  
            A12((n-1)*m+in-n,in) = 1;
            A12((in-floor(in/n)*n-2)*m+floor(in/n)+1,in) = 1;
            A12((in-floor(in/n)*n-2)*m+floor(in/n)+m+1,in) = -1;
            A12((n-1)*m+in,in) = -1;
        end
    end
end