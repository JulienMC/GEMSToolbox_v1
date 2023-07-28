function [q] = groundwater_st(q, x, V_gw, A_n)
    % A function to generate the source terms from the host rock
    % it takes in the source array, the coordinates of the nodes
    % and the darcy velocity of the regional groundwater into the mine

    xmax = max(x(:,1));
    xmin = min(x(:,1));
    ymax = max(x(:,2));
    ymin = min(x(:,2));
    
    % A_n m^2 nodal area of influence - it is a function of the mine working's geometry.
    q = q + V_gw .* A_n .* lin_scaler(xmax, xmin, ymax, ymin, x);

    fprintf("Sum q = %e\n", sum(q));
end

% A function which scales the applied groundwater flow by multiplying it by
% a value in the interval -1 to 1.
function mod = lin_scaler(xmax, xmin, ymax, ymin, x)
    medx = median(x(:,1));
    r = max([medx - xmin, xmax - medx]);
    rmin = medx - r;
    rmax = medx + r;
    normx = (x(:,1) - rmin) / (rmax - rmin);

    b = zeros(length(normx),1)+1;
    n = sum(b(x(:,1)<median(x(:,1))));
    
    mod = linear_scale(normx, -1);%sig_scale(normx, -0.5, 11);% xmod to be multiplied by projected x
    %sum(mod)
    if sum(mod) > 0
       mod(mod > 0) = mod(mod > 0) + sum(mod)/n; 
    elseif sum(mod) < 0
        mod(mod < 0) = mod(mod < 0) - sum(mod)/n; 
    end
    
    %mod(:,2) = (cos(normy*pi-pi/2)*0.75-cos(normy*pi*5-pi/2)/4)*0.5; % ymod to be multiplied by projected y area
end

function mod = sin_scale(normx)
    mod = sin(normx*pi+pi/2);
end

function mod = sig_scale(normx, offset, mult) % offset = -0.5, mult=11
    mod = 1./(1+exp(-(normx+offset)*mult));
end

function mod = linear_scale(normx, dir)
    mod = (normx * 2 - 1) * dir;
end