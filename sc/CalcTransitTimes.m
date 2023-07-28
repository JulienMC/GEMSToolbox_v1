function [tt] = CalcTransitTimes(npipes,node_pipes_in,pipe_nodes,n_tree,L,v,Tr_st)
    ntot = size(npipes,2);
    tt = zeros(ntot,1); % stores the transit times for each pipe (i.e. shortest time from closest thermal interference to pipe start)
    for i = 1:length(n_tree)
        in = n_tree(i);
        if in == 0 % end of tree reached
            break
        end
        % Calculate the minimum total time for water to flow from closest
        % point of same initial temperature to the start of the pipe.
        % find the shortest time from incoming pipes 
        pipes = transpose(node_pipes_in(1:npipes(1,in), in));
        min_t = inf; % sets a very large starting minimum value
        for j = 1:length(pipes)
            ip = pipes(j);
            fn = pipe_nodes(ip,1); % the inflow node at the other end of the pipe

            % determine if an earlier thermal interference can arise from
            % initial water temperature differences
            reset = ~(Tr_st(fn) == Tr_st(in));

            % transit time to in
            if reset
                % if we need to reset then it is the time it takes to flow
                % through the pipe linking fn and in.
                tt_in = L(ip)/v(ip);
            else
                % else it is the transit time to fn plus the time it takes
                % to flow through the pipe.
                tt_in = (L(ip)/v(ip) + tt(fn));
            end

            if tt_in < min_t
                min_t = tt_in;
            end
        end
        tt(in) = min_t;
    end
end