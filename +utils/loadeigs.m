function states = loadeigs(flist)
%LOADEIGS この関数の概要をここに記述
%   詳細説明をここに記述
    [~, reindex] = sort( str2double( regexp( {flist.name}, '\d+', 'match', 'once' )));
    list = flist(reindex);
    if length(list) == 0
        error("files are not found.")
    end
    path = sprintf("%s/%s", list(1).folder, list(1).name);
    [dtype, dim, domain, hbar, basis, eval] = utils.readheader(path);
    states = [];
    sysinfo = SystemInfo(dim, domain);
    for i = 1:length(list)
        path = sprintf("%s/%s", list(i).folder, list(i).name);
        data = utils.readdata(path, [dim, 4]);
        vec = data(:,3) + 1i*data(:,4);
        s = State(sysinfo, basis, vec, eval);
        states = [states; s];
    end
end

