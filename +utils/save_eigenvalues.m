function save_eigenvalues(obj, evals, varargin)
  
  par = inputParser;
  addRequired(par, 'obj', @isobject);
  addRequired(par, 'evals', @isnumeric);  
            
  addOptional(par, 'savedir', 'Data', @(x) isstring(x) | ischar(x) );
  addOptional(par, 'filename', 'eigen_evals.dat', @(x) isstring(x) | ischar(x) );
  addOptional(par, 'trim', false);
            
  parse(par, obj, evals, varargin{:} );
  obj      = par.Results.obj;
  evals    = par.Results.evals;
  savedir  = par.Results.savedir;
  filename = par.Results.filename;
  trim     = par.Results.trim;
  
  if ~exist(savedir, 'dir')
      mkdir(savedir);
  end
  
  [row, col] = size(evals);
  
  if row == col
      evals = diag(evals);
  end

  dtype = class(evals);
  
  path = strcat(savedir, '/', filename);
  header = utils.scaleinfo2str(obj);

  of = fopen(path, 'w');
  fprintf(of, '%s', header);
  fprintf(of, "# date : %s\n", date);
  fprintf(of, "# n, \treal(eval),\t imag(eval)\n");
  data = [1:length(evals); real(evals).'; imag(evals).'];  

  if strcmp(dtype, 'double')
      fmt = '%d\t%.18e\t%.18e\n';      
  elseif strcmp(dtype, 'mp')
      if (~trim)
          fmt = '%d\t%s\t%s\n';
      else
          fmt = sprintf("%%d\\t%%.%ds\\t%%.%ds\\n", trim, trim);
      end
  else
      error("something wrong");
  end
  fprintf(of, fmt, data);    
  fclose(of);
end


