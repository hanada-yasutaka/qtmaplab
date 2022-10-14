function data = readdata(path, dsize)

of = fopen(path);

line = fgetl(of);
line = line(find(~isspace(line))); % remove white space
keys = strsplit( line(2:end), '=') ;
dtype = cell2mat(keys(2));

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

