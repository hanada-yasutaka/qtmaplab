function [varargout] = readheader(path)

of = fopen(path);

line = fgetl(of);
line = line(find(~isspace(line))); % remove white space
keys = strsplit( line(2:end), '=') ;
dtype = cell2mat(keys(2));
basis = '';
eigenvalue = '';

while ~feof(of)
    line = fgetl(of);
    if line(1) == '#' & contains(line, '=')
        line = line(find(~isspace(line))); % remove white space
        keys = strsplit( line(2:end), '=') ;
        var = cell2mat(keys(1));
        val = cell2mat(keys(2));
        if strcmp(dtype, 'mp')
            try 
                expr = sprintf("%s = mp('%s');", var, val);
                eval(expr);
            catch 
                expr = sprintf("%s = %s;", var, val);
                eval(expr); 
            end
        else
            expr = sprintf("%s = %s;", var, val);
            eval(expr);
        end

    elseif line(1) ~= '#'
        break
    end
end
fclose(of);

dtype = dtype;
dim = dim;
domain = [qmin qmax; pmin pmax];
hbar = hbar;
basis = basis;
eigenvalue = eigenvalue;
varargout = {dtype; dim; domain; hbar; basis; eigenvalue};
end

