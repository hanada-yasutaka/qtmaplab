function evalus = loadeigenvalues(path)
%LOAD_EIGENVALUES この関数の概要をここに記述
%   詳細説明をここに記述
of = fopen(path, "r");
str = strsplit(fgetl(of), " : ");
dtype = cell2mat(str(end) );

for i = 1:4
    str = strsplit(fgetl(of), " : ");
    
    if strcmp(dtype, 'double')
        var = str2num( cell2mat( str(end) ) );
    elseif strcmp(dtype, 'mp'
        var = mp(sprintf('%s', cell2mat( str(end) ) ) )
    else
        error(sprintf("data type : %s does not support", dtype);
    end
end

fclose(of)
end

