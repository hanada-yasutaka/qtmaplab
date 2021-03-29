function states = eigs2states(obj, evecs, evals)
    states = [];
    dim = obj.dim;
    system = SystemInfo(obj.dim, obj.domain, class(obj) );
    if exist('evals', 'var')
        if isequal( size(evals), [dim dim] )
            evals = diag(evals);
        end
    end
    [sr,sc] = size(evecs);
    for i = 1:sc
        if exist('evals', 'var')
            stat = FundamentalState(system, obj.basis, evecs(:, i), evals(i) );
        else
            stat = FundamentalState(system, obj.basis, evecs(:, i));            
        end
        states = [states; stat];
    end
end

