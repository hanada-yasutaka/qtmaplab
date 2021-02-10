function y = scaleinfo2str(obj)
%GENERATE_SAVEFILE_HEADER この関数の概要をここに記述
%   詳細説明をここに記述
y = "";
y = y + sprintf("# data type : %s\n", class(obj.domain));
y = y + sprintf("# qmin : %e\n", obj.domain(1,1));
y = y + sprintf("# qmax : %e\n", obj.domain(1,2));
y = y + sprintf("# pmin : %e\n", obj.domain(1,1));
y = y + sprintf("# pmax : %e\n", obj.domain(1,2));
y = y + sprintf("# dim : %d\n", obj.dim);
y = y + sprintf("# hbar : %.18f\n", obj.hbar);
y = y + sprintf("# machine eps : %.18e\n", obj.eps);
y = y + sprintf("# date : %s\n", date);
end

