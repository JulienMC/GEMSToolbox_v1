function [stable_network, stable_xtotal] = cleanNetworkMatrix(A120, xtotal)

    prev_size = 0;
    current_size = size(A120,1);
    while ~(current_size == prev_size)
        current_size = size(A120,1);
        % find no-flow pipes (i.e. pipes connected to a single node)
        sn_idx = findColumnsBySum(abs(A120), 1);
        matchingRows = any(abs(A120(:,sn_idx)) == 1,2);
        A120 = A120(~matchingRows,:);
        prev_size = size(A120,1);
    end

    % remove duplicate pipes
    A120 = removeDuplicateRows(A120);
    % find loose nodes column indices
    ln_idx = findColumnsBySum(abs(A120), 0);
    % Remove the columns with the specified indices
    A120(:,ln_idx) = [];
    xtotal(ln_idx,:) = [];

    stable_network =  A120;
    stable_xtotal = xtotal;
end