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

function rs_crit = critical_radii(ip, pipe_nodes, xtotal, d, max_radius, int_d, pipe_centres, xs1, ys1, zs1, np)
    % scaling factor to account for interference of pipes from all sides
    % in effect, the ratio of a dual inverted cone to that of the cylinder
    scl = 1/2;%1/3; % formally was 1/2 as assumed interference started at half the distance to opposite pipe
    
    % Calculate the radii for the pipe ip
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
    N = s .* vecnorm(w, 2, 2);
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

    rs_crit = radial_integration([m_pos*scl,m_neg*scl],deg2rad([ca_pos,ca_neg]),max_radius,int_d);
    % makes sure that the r_crit is not lower than the pipe's own radius
    % and if so adds +10%!
    if rs_crit <= d(ip)/2
        rs_crit = d(ip)/2*1.1;
    end    

end

% % Moved to external file
% % A scaling function to approximate the critical radius of a pipe in 3D
% % r_cs = (2,1); contains the min and max critical radii for that pipe
% % r_ca = (2,1); contains the matching angles at which these critical radii
% % are observed.
% % Density = int; is the density of the integration a density of 100 offers
% % a stable solution.
% % r_as = array; all the radii at 360 angle around the pipe
% function r_as = radial_integration(r_cs, a_cs, r_inf, density)
% %density = 1000; for accurate results
% verbose = 0;
% r_cmin = r_cs(1);
% r_cmax = r_cs(2);
% a_cmin = a_cs(1);
% a_cmax = a_cs(2);
% a_mean = mean(a_cs);
% a_means = [a_mean, a_mean+pi];
% angles = linspace(0,2*pi,density)+ a_cmin;
% r_as = nan(density,1); % to store the representative radii at various integraiton angles
% 
%     for i = 1:density
%         a = angles(i);
%         angles(i) = a;
%         if (a >= a_cmin && a <= a_means(1)) || (a > a_means(2) && a <= a_cmin + 2*pi)
%             % region where to use scaled r_cmin
%             r_as(i) = abs(r_cmin / cos(a - a_cmin));
%         elseif (a > a_means(1) && a <= a_cmax) || (a > a_cmax && a <= a_means(2)) 
%             % region where to use scaled r_cmax
%             r_as(i) = abs(r_cmax / cos(abs(a_cmax - a)));
%         end
%         if (a > a_cmax + pi/2 && a < a_cmin + 3/2 * pi && a_cmax - a_cmin < pi) || (a > a_cmin + pi/2 && a < a_cmax - pi/2 && a_cmax - a_cmin > pi) 
%             % define the infinite heat supply region when the interference
%             % angles are not at 180 degrees.
%             r_as(i) = r_inf;
%         end
%         r_as(i) = min(r_as(i),r_inf);
%     end
% end


 