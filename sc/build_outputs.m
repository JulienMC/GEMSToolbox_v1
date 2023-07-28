


function Outputs = build_outputs(PhysicalProperties, NumericalProperties, Outputs)

    n_pairs = split(NumericalProperties.node_outputs,',');
    e_pairs = split(NumericalProperties.pipe_outputs,',');

    % Adding additional node outputs
    n_outputs = cell(length(n_pairs),1);
    for i = 1:length(n_pairs)
        vars = split(n_pairs{i,:});
        % ensure consitent output data dimensions
        if size(Outputs.(vars{1}),1) < size(Outputs.(vars{1}),2)
            Outputs.(vars{1}) = Outputs.(vars{1})';
        end
        A = Outputs.(vars{1});
        if ispc % If on windows Paraview cannot read nans, so replace nans with -999
            A(isnan(A)) = -999;
        end
        n_outputs{i} = {A string(vars{2})};
    end

    % Adding additional pipe outputs
    e_outputs = cell(length(e_pairs),1);
    for i = 1:length(e_pairs)
        vars = split(e_pairs{i,:});
        % ensure consitent output data dimensions
        if size(Outputs.(vars{1}),1) < size(Outputs.(vars{1}),2)
            Outputs.(vars{1}) = Outputs.(vars{1})';
        end
        A = Outputs.(vars{1});
        if ispc % If on windows Paraview cannot read nans, so replace nans with -999
            A(isnan(A)) = -999;
        end
        e_outputs{i} = {A string(vars{2})};
    end

    Outputs.('n_outputs') = cat(1, n_outputs{:});
    Outputs.('e_outputs') = cat(1, e_outputs{:});
end