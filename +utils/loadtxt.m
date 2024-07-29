function data = loadtxt(path, dsize, dtype)

    of = fopen(path);
    data = zeros(dsize(1), dsize(2), dtype);
    i = 1;
    while ~feof(of)
        line = fgetl(of);
        if line(1) ~= '#'
            if strcmp(dtype, 'mp')
                strexp = regexprep(line, '\t', ',');
                str2mp = mp(sprintf('[%s]', strexp));
                data(i, :) = str2mp;
            else
                strexp = regexprep(line, '\t', ',');
                data(i, :) = str2num(strexp);
            end
            i = i + 1;
        end
    end
end

