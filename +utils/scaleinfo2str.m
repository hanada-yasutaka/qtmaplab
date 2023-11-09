function y = scaleinfo2str(obj)
%GENERATE_SAVEFILE_HEADER この関数の概要をここに記述
%   詳細説明をここに記述
y = "";
y = y + sprintf('# dtype = %s\n', class(obj.domain));
y = y + sprintf('# eps = %.18e\n', obj.eps);
y = y + sprintf('# dim = %d\n', obj.dim);
y = y + sprintf('# qmin = %e\n', obj.domain(1,1));
y = y + sprintf('# qmax = %e\n', obj.domain(1,2));
y = y + sprintf('# pmin = %e\n', obj.domain(2,1));
y = y + sprintf('# pmax = %e\n', obj.domain(2,2));
y = y + sprintf('# hbar = %.18f\n', obj.hbar);
end

