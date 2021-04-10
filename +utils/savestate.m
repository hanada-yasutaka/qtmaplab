function savestate(obj, state, filename, varargin)
  par = inputParser;
  addRequired(par, 'obj', @isobject);
  parse(par, obj);
  obj = par.Results.obj;
  
  addRequired(par, 'state', @(x) isa(x,'FundamentalState'));  
  addRequired(par, 'filename', @(x) isstring(x) | ischar(x) );
            
  addOptional(par, 'savedir', 'Data') %@(x) isstring(x) | ischar(x) );
  addOptional(par, 'trim', false);
  addOptional(par, 'basis', obj.basis, @(x) isstring(x) | ischar(x) );    
            
  parse(par, obj, state, filename, varargin{:} );
  state    = par.Results.state;
  savedir  = par.Results.savedir;
  filename = par.Results.filename;
  trim     = par.Results.trim;
  basis    = par.Results.basis;
  
  if ~exist(savedir, 'dir')
      mkdir(savedir);
  end
  
  dtype = class(state.dtype);
  
  path = strcat(savedir, '/', filename);
  header = utils.scaleinfo2str(obj);

  of = fopen(path, 'w');
  fprintf(of, '%s', header);

  if strcmp(obj.dtype, 'mp')
      fprintf(of, '# eigenvalue = %s\n', state.eigenvalue);
  else
      fprintf(of, '# eigenvalue = %.18e\n', state.eigenvalue);
  end
  fprintf(of, '# basis = ''%s''\n', basis);
  fprintf(of, '# x,\t|<x|psi>|^2,\tRe[<x|psi>],\tIm[<x|psi>]\n');
  data = [state.x.'; abs2(state.y).'; real(state.y).'; imag(state.y).'];  

  if strcmp(obj.dtype, 'double')
      fmt = '%.18f \t%.18e \t%.18e \t%.18e\n';      
  elseif strcmp(obj.dtype, 'mp')
      if (~trim)
          fmt = '%s \t%s \t%s \t%s\n';
      else
          fmt = sprintf("%%.%ds \\t%%.%ds \\t%%.%ds \\t%%.%ds\\n", trim, trim, trim, trim);
      end
  else
      error("something wrong");
  end
  fprintf(of, fmt, data);    
  fclose(of);
end
