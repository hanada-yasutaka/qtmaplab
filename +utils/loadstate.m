function state = loadstate(path)
%LOADSTATE この関数の概要をここに記述
%   詳細説明をここに記述
    impd = importdata(path);
    length(impd.textdata);
    of = fopen(path, 'r');
    tline = fgetl(of);
    eval(strcat(tline(2:end), ';') ); % dtype
    %tline = fgetl(of); % date
    for i=1:length(impd.textdata)-2
        tline = fgetl(of);
        if tline(1) ~= '#'
            break
        end
    
        if strcmp(dtype, 'mp')
            strsp = strsplit(tline(2:end), ' = ');
            try
                eval(sprintf('%s = mp(''%s'');', strsp{1}, strsp{2} ));
            catch ME
                eval(strcat(tline(2:end), ';'));
            end
        else
            eval(strcat(tline(2:end), ';'));
        end
    end

    %strdata = strings(size(f.data));
    if strcmp(dtype, 'double')
        data = impd.data;
    else
        data = zeros(size(impd.data), dtype);
        i=1;
        while ~feof(of)
            tline = fgetl(of);
            if tline(1) ~= '#'
                strsp = regexprep(tline, '\t', ',');
                str2mp = mp(sprintf('[%s]', strsp));
                data(i,:) = str2mp;
                i = i + 1;
            end
        end
    end
    
    domain = [qmin qmax; pmin pmax];
    sysinfo = SystemInfo(dim, domain);
    vec = data(:,3) + 1i*data(:,4);
    state = FundamentalState(sysinfo, basis, vec, eigenvalue);
end

