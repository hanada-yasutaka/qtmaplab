function saveeigs(obj, evals, states, varargin)
    par = inputParser;
    addRequired(par, 'obj', @isobject);
    parse(par, obj);
    obj = par.Results.obj;
    addRequired(par, 'evals', @isnumeric);    
    addRequired(par, 'states', @(x) isnumeric(x) | isa(x, 'FundamentalState') );
    addOptional(par, 'basis', obj.basis, @(x) isstring(x) | ischar(x) );
    addOptional(par, 'savedir', 'Data', @(x) isstring(x) | ischar(x) );
    addOptional(par, 'trim', false, @islogical);
       
    parse(par, obj, evals, states, varargin{:});
    %obj    = par.Results.obj;
    states  = par.Results.states;
    basis = par.Results.basis;
    savedir = par.Results.savedir;
    trim = par.Results.trim;
    
    if ~isa(states, 'FundamentalState')
        states = eigs2states(obj, evals, states);
    end
    
    utils.save_eigenvalues(obj, evals, 'savedir', savedir, 'trim', trim);
        
    for i=1:obj.dim
        s = states(i);
        filename = sprintf('eigen_%srep_n%d.dat', basis, i);        
        if strcmp(basis, 'q')
            utils.savestate(obj, s.qrep(), filename, varargin{:});
        else
            utils.savestate(obj, s.prep(), filename, varargin{:});
        end
        disp(filename);
    end
end

