function saveeigs(obj, evals, states, varargin) 
    par = inputParser;
    addRequired(par, 'obj', @isobject);
    addRequired(par, 'evals', @isnumeric);    
    addRequired(par, 'states', @(x) isnumeric(x) | isa(x, 'FundamentalState') );    
    %parse(par, obj, evals,states);
    %obj = par.Results.obj;
    
    addOptional(par, 'header', 'eigen' ); %, @(x) isstring(x) | ischar(x) );     
    addOptional(par, 'basis', obj.basis);
    addOptional(par, 'savedir', 'Data');%, @(x) isstring(x) | ischar(x) );
    addOptional(par, 'trim', false);
    addOptional(par, 'varbose', true);
       
    parse(par, obj, evals, states, varargin{:});            
    header = par.Results.header;
    basis = par.Results.basis;
    savedir = par.Results.savedir;    
    trim = par.Results.trim;
    varbose = par.Results.varbose;
    
    if ~isa(states, 'FundamentalState')
        states = eigs2states(obj, evals, states);
    end
    
    fname = sprintf('%s_evals.dat', header);
    utils.save_eigenvalues(obj, evals, 'savedir', savedir, 'filename', fname, 'trim', trim, 'varbose', varbose);
    for i=1:obj.dim
        s = states(i);
        fname = sprintf('%s_%srep_n%d.dat', header, basis, i);
        if strcmp(basis, 'q')
            utils.savestate(obj, s.qrep(), fname, 'savedir', savedir, 'trim', trim, 'basis', basis, 'varbose', varbose);
        else
            utils.savestate(obj, s.prep(), fname, 'savedir', savedir, 'trim', trim, 'basis', basis, 'varbose', varbose);
        end
    end
end

