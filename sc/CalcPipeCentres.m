function [pipe_centres, xs1, ys1, zs1, np] = CalcPipeCentres(pipe_nodes, xtotal)
    np = size(pipe_nodes,1);

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
end