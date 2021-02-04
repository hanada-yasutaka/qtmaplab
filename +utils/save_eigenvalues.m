function save_eigenvalues(obj, evals, dirname, filename, trim)
  arguments
      obj;
      evals {mustBeNumeric};
      dirname = '.';
      filename = 'eigen_evals.dat';
      trim {mustBeInteger} = false;
  end
  
  if ~exist(dirname, 'dir')
      mkdir(dirname);
  end
  
  [row, col] = size(evals);
  
  if row == col
      evals = diag(evals);
  end

  dtype = class(evals);
  
  path = strcat(dirname, '/eigen_evals.dat');
  header = utils.scaleinfo2str(obj);

  of = fopen(path, 'w');
  fprintf(of, '%s', header);
  fprintf(of, "# n, \treal(eval),\t imag(eval)\n");
  data = [0:length(evals)-1; real(evals).'; imag(evals).'];  

  if strcmp(dtype, 'double')
      fmt = '%d\t%.18e\t%.18e\n';      
  elseif strcmp(dtype,'mp')
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


