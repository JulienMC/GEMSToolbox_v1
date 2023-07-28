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

function r_crit = critical_radii_dot(pipe_nodes, xtotal, A12, A10, d, PhysicalProperties, NumericalProperties)
    tic
    verbose = 0;
    max_radius = PhysicalProperties.max_radius; % default is nothing is found
    int_d = NumericalProperties.int_d; % integration density
    np = size(pipe_nodes,1);
    r_crit = size(np,1);

    xs1 = xtotal(pipe_nodes(:,1),1); % xs of starting nodes
    ys1 = xtotal(pipe_nodes(:,1),2); % ys of starting nodes
    if size(xtotal,2) < 3 % handles 2D case
        zs1 = zeros(size(pipe_nodes(:,1))); % zs of starting nodes
    else
        zs1 = xtotal(pipe_nodes(:,1),3); % zs of starting nodes
    end

    xs2 = xtotal(pipe_nodes(:,2),1); % xs of end nodes
    ys2 = xtotal(pipe_nodes(:,2),2); % ys of end nodes
    if size(xtotal,2) < 3
        zs2 = zeros(size(pipe_nodes(:,1))); % zs of starting nodes
    else
        zs2 = xtotal(pipe_nodes(:,2),3); % zs of end nodes
    end

    pipe_centres = zeros(np,3);
    pipe_centres(:,1) = mean([xs1,xs2],2);
    pipe_centres(:,2) = mean([ys1,ys2],2);
    pipe_centres(:,3) = mean([zs1,zs2],2);

    for ip = 1:np
        ipc = pipe_centres(ip,:); % ip's centre coordinates
        norm_ipc = [xs1(ip),ys1(ip),zs1(ip)] - ipc; % vector from centre point of ip to one of its ends
        li = [[1:ip-1],[ip+1:np]]; % all indices except ip's
        
        pc = pipe_centres(li,:); % pipe centres except ip's
        norm_pc = pc - ipc; % normalised pipe centres

        % distances between all the pipe_centres and ip's centre
        dists = sqrt((pc(:,1) - pipe_centres(ip,1)).^2 + ...
                     (pc(:,2) - pipe_centres(ip,2)).^2 + ...
                     (pc(:,3) - pipe_centres(ip,3)).^2);
        
        % calculate angles between the pipe_centres  
        u = repmat(norm_ipc,np-1,1);
        v = norm_pc;
        w = cross(u,v,2);
        n = null([norm_ipc(:).'])'; % [rand() rand() rand()];
        n = repmat(n(end,:),np-1,1);
        s = sign(dot(w,n,2));
        N = s .* vecnorm(w, 2, 2);% s .* 
        D = dot(u, v, 2);
        angles = atan2d(N,D);
        pos_a = angles; neg_a = angles;
        pos_a(angles<0) = 1000; % these values can never be reached in the for loop below
        neg_a(angles>0) = -1000; % these values can never be reached in the for loop below

        % get minimum distance and their index where the angles are between
        % 45 and 135 degrees
        %size(dists((angles > 80) + (angles < 100) == 2));
        cone = 60;
        minA = 90-(cone/2);
        maxA = 90+(cone/2);
        found = [0 0];
        [m_pos, idx_pos] = min(dists(((pos_a > minA) + (pos_a < maxA)) == 2));
        [m_neg, idx_neg] = min(dists(((neg_a > -maxA) + (neg_a < -minA)) == 2));
        if ~isempty(m_pos)
            found(1,1) = 1;
        end
        if ~isempty(m_neg)
            found(1,2) = 1;
            %break
        end

        minA = minA - 5;
        maxA = maxA + 5;
        minA = max(0, minA);
        maxA = min(180, maxA);
        
        % assign found value angles
        % filter out values outside range of interest so that the idx_
        % match the array size on which the comparision is done in the for
        % loop when m_pos and m_neg are obtained.
        pos_a = pos_a(((pos_a > minA) + (pos_a < maxA)) == 2);
        neg_a = neg_a(((neg_a > -maxA) + (neg_a < -minA)) == 2);
        ca_pos = pos_a(idx_pos);
        ca_neg = neg_a(idx_neg)+360;

        % if nothing is found one side or the other, then assign a lot of
        % space to it
        if found(1,1) == 0
            m_pos = max_radius*2;
            if found(1,2) == 0 % if neg is also not found assume angles are 180 opposits
                ca_pos = 0;
            else % if the neg angle is found, then assume we are opposit to it
                ca_pos = ca_neg-180;
            end
        end
        if found(1,2) == 0
            m_neg = max_radius*2;
            if found(1,1) == 0 % if pos is also not found assume angles are 180 opposits
                ca_neg = 180;
            else % if the pos angle is found, then assume we are opposit to it
                ca_neg = ca_pos+180;
            end
        end
        %r_crit(ip,1) = (m_pos+m_neg)/(length(m_neg)+length(m_pos)) .* 0.5;
        r_crit(ip,1) = radial_approximation([m_pos/2,m_neg/2],deg2rad([ca_pos,ca_neg]),max_radius,int_d);
        % makes sure that the r_crit is not lower than the pipe's own radius
        % and if so adds +10%!
        if r_crit(ip,1) <= d(ip)/2
            r_crit(ip,1) = d(ip)/2*1.1;
        end    
    end
    toc
end
 