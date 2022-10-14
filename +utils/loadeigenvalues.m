function evals = loadeigenvalues(path)
%LOAD_EIGENVALUES この関数の概要をここに記述

    [dtype, dim] = utils.readheader(path);
    data = utils.readdata(path, [dim, 3]);
    evals = data(2) + 1i*data(3);

end

