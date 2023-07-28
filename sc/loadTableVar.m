function vals = loadTableVar(table,rows,vars, tableD, varargin)
    try
        vals = table{rows,vars};
    catch exception
        try
        vals = tableD{1,vars};
        msg = exception.message;
        info = split(msg,"'");
        missTag = string(info(2));
        warning("Scenario Input Variable '%s' not found. Using default value of %s instead.", missTag, string(vals))
        catch exception
            msg = exception.message;
            info = split(msg,"'");
            missTag = string(info(2));
            error("Scenario input variable '%s' does not exist in the specified file, nor in the default. Please check the spelling of your variables in your Scenario file.", missTag)
        end
    end
    forceType = "";
    if length(varargin) > 0
        forceType = varargin{1};
    end
    if forceType == "str"
        vals = vals;
    else if iscell(vals)
            if input_eval_check(vals{1}) == 1
                vals = str2num(cell2mat(vals));
            else
                error("The %s input contained some characters which are not allowed. Only numbers, '.', ';', '_','[' and ']', are allowed. Please check your input file.", vars)
            end
    end
end