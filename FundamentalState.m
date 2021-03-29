classdef FundamentalState < matlab.mixin.SetGet & SystemInfo
    % FundamentalState provide a class of a 1-dim. wavefunction with fundamental basis e.g, position $q$ and momentum $p$ basis.
    properties (SetAccess = protected)
        sysinfo % instance of the class SystemInfo        
        basis {mustBeMember(basis,{'q','p'})} = 'p' % basis (representation) of wavefunction                 
        x % coordinate: $q$ or $p$
        y (:,1) {mustBeNumeric} = [] % <x|¥psi> where x is $q$ or $p$
        eigenvalue {mustBeNumeric} = NaN % eigenvalue data        
    end
    
    
    methods
        function obj = FundamentalState(sysinfo, basis, vec, eigenvalue)
            if ~isa(sysinfo, 'SystemInfo')
                error("scl must be SystemInfo class")
            end
            
            obj@SystemInfo(sysinfo.dim, sysinfo.domain)     
            obj.sysinfo = sysinfo;
            
            obj.basis = basis;            
            
            if strcmp(basis, 'q')
                obj.x = obj.q ;
            elseif strcmp(basis, 'p')
                obj.x = obj.p;
            end
                                    
            if exist('vec', 'var')
                obj.y = vec;
            end
            
            if exist('eigenvalue', 'var')
                obj.eigenvalue = eigenvalue;
            end
        end
        
        function obj = tostate(obj, vec)
            % convert to vector array to state             
            obj = FundamentalState(obj.sysinfo, obj.basis, vec);
        end
        
        function obj = double(obj)
            % truncate to double 
            sysinfo = SystemInfo(obj.dim, double(obj.domain) );
            obj = FundamentalState(sysinfo, obj.basis, double(obj.y));            
        end
        
        function obj = q2p(obj)
            assert( strcmp(obj.basis, 'q'), 'this state basis is not "q"');
            if abs(obj.domain(2,1)) == abs(obj.domain(2,2))
                vec = fftshift( fft(obj.y) );
            else
                vec = fft(obj.y);
            end
            
            obj = FundamentalState(obj.sysinfo, 'p', vec/sqrt(obj.dim));
        end
        
        function obj = p2q(obj)
            assert( strcmp(obj.basis, 'p'), 'this state basis is not "p"');            
            if abs(obj.domain(2,1)) == abs(obj.domain(2,2))
                vec = ifft( ifftshift(obj.y) );                
            else
                vec = ifft(obj.y);                                

            end
            %vec = ifft(obj.y);                
            obj = FundamentalState(obj.sysinfo, 'q', vec*sqrt(obj.dim));
        end
        
        function obj = qrep(obj)
            if strcmp(obj.basis, 'q')
                obj = obj;
            else
                obj = obj.p2q();
            end
        end
        
        function obj = prep(obj)
            if strcmp(obj.basis, 'p')
                obj = obj;
            else
                obj = obj.q2p();
            end
        end

        function [X, Y, mat] = hsmrep(varargin)
            
            par = inputParser;
            addRequired(par, 'obj', @isobject);
            parse(par, varargin{1});
            obj = par.Results.obj;
                                    
            addOptional(par, 'gridnum',       50 )%, @(x) all(x>0 & x == floor(x)) );
            addOptional(par, 'vrange', obj.domain)%, @(x) isnumeric(x) & all(size(x) == [2 2], 'all') );
            addOptional(par, 'ismp',       false )%, @islogical);
                                        
            parse(par, varargin{:} );
            vrange  = par.Results.vrange;
            gridnum  = par.Results.gridnum;                                   
            ismp = par.Results.ismp;
            
            assert(vrange(1,2) > vrange(1,1), 'vrange(1,1) >= vrange(1,2)');
            assert(vrange(2,2) > vrange(2,1), 'vrange(2,1) >= vrange(2,2)');
            
                                    
            if length(gridnum) == 1
                gridnum = [gridnum gridnum];
            end                        
                        
            if strcmp(obj.basis, 'q')
                v = obj.y;                
            elseif strcmp(obj.basis, 'p')
                vec = obj.p2q();
                v = vec.y;
            end
            
            v = v / sqrt( dot(v,v) );
            
            if ismp
                x = obj.q;
                if ~strcmp( class(vrange), 'mp')
                    msg = {"The computation precision is set as mp, but your the input arguments"
                        sprintf("vrange are %s", class(vrange))
                        "This is not safe to keep preciison."
                        };
                    warning(sprintf("%s\n",msg{:}))
                end
                hbar = obj.hbar;
                qc = linspace(vrange(2,1), vrange(2,2), gridnum(1));
                pc = linspace(vrange(2,1), vrange(2,2), gridnum(2));
                
                [X,Y] = meshgrid(qc,pc);              
                mat = zeros(gridnum(1), gridnum(2), class(v));
                if obj.periodic
                    csfun = @(qc, pc) csp(x, hbar, obj.domain(1,:), qc, pc, 5);
                else
                    csfun = @(qc, pc)  cs(x, hbar, qc, pc);
                end                
            else
                x = double(obj.q);              
                %x = linspace(double(vrange(1,1)), double(vrange(1,2)), obj.dim);
                v = double(v);
                vrange = double(vrange);
                
                hbar = double(obj.hbar);
                
                %dq = double( ( vrange(1,1) - vrange(1,2) )/ gridnum(1) );
                %dp = double( ( vrange(2,1) - vrange(2,2) )/ gridnum(2) );               
                
                %qc = vrange(1,1): dq : vrange(1,2);
                %pc = vrange(2,1): dp : vrange(2,2);               
                
                qc = linspace(vrange(1,1), vrange(1,2), gridnum(1));
                pc = linspace(vrange(2,1), vrange(2,2), gridnum(2));
                
                [X,Y] = meshgrid(qc, pc);
                mat = zeros(gridnum(1), gridnum(2), 'double');
                
                if obj.periodic
                    csfun = @(qc, pc) csp(x, hbar, double(obj.domain(1,:)), qc, pc, 5);
                else
                    csfun = @(qc, pc)  cs(x, hbar, qc, pc);
                end
                
            end
                                   
            for i = 1:gridnum(1)
                for j = 1:gridnum(2)
                    csv = transpose( csfun(qc(i), pc(j)) );
                    % dot(c, v) = abs2( sum( conj(c) .* v) )
                    mat(j,i) = abs2( dot(csv, v) );
                end
            end
        end
        
        
        function y = abs2(obj)
            % return |<x|psi>|^2
            y = abs( obj.y .* conj(obj.y) );
        end
        
        function y = norm(obj)
            % overload built-in norm function to  <¥psi|¥psi>
            y = norm( obj.y );
        end
        
        function display(obj)
            % overload built-in display function to display <x|¥psi>
            display( obj.y )
        end
        
        function disp(obj)
            %overload built-in disp function to display <x|¥psi>            
            disp( obj.y );
        end
        
        function info(obj)
            builtin('display', obj);
        end
        
        
        function obj = plus(v1, v2)
            % overload built-in plus ( + ) function to |v1> + |v2>
            basisconsistency(v1, v2);            
            if strcmp(class(v2), 'FundamentalState')
                v2 = v2.y;
            end
            v3 = v1.y + v2;
            obj = FundamentalState(v1.sysinfo, v1.basis, v3);
        end
        
        function obj = minus(v1, v2)
            % overload built-in minus ( - ) function to |v1> - |v2>            
            basisconsistency(v1,v2);            
            if strcmp( class(v2), 'FundamentalState')
                v2 = v2.y;
            end
            v3 = v1.y - v2;
            obj = FundamentalState(v1.sysinfo, v1.basis, v3);
        end
                    
        function scalar = inner(v1, vm)
            % return < v1 | vm > if ins isvector
            % if vm is matrix such that
            % vm = [a, b, c, ...]
            % then inner returns vector for the inner product
            % [ <v1|a>, <v1|b>, <v1|c> ]
            
            if strcmp( class(vm), 'FundamentalState')
                vm = vm.y;
            end
            
            if isvector(vm)
                scalar = dot(v1.y, vm);
            elseif ismatrix(vm)
                scalar = conj(v1.y).' * vm;
            else
                error("mv must be vector or matrix");
            end
        end
        
        function obj = times(v1, a)
            % overload built-in times ( .* ) (element wise multiple) function to a * |v1> where "a" must be scalar
            if ~isscalar(a)
                error('"a" must be scalar')
            end            
            obj = FundamentalState(v1.sysinfo, v1.basis, v1.y .* a);
        end
        
        
        function obj = mrdivide(v1, a)
            % return v1/a
            if ~isscalar(a) && a ~= 0
                error('"a" must be non-zero scalar')
            end
            obj = FundamentalState(v1.sysinfo, v1.basis, v1.y/a);
        end            
            
            
                
        %function v = inners(v1, vecs)
            % return [<v1 |vecs(1)>, <v1|vecs(2)>, ... ]
        %    v = transpose( conj(v1.y) ) * vecs
        %end
            
        
        function obj = normalize(obj)
            v = obj.y;
            a = dot(v, v);
            obj = FundamentalState(obj.sysinfo, obj.basis, v/sqrt(a));
        end
        
        function obj = coherent(obj, qc, pc, varargin)
            
            par = inputParser;
            addRequired(par, 'obj', @isobject);                        
            addRequired(par, 'qc',  @isreal);
            addRequired(par, 'pc',  @isreal);

            default = true;
            addOptional(par,'periodic', default, @islogical);
            addOptional(par,'verbose',  default, @islogical);
            pnum = 5;
            addOptional(par,'pnum',  pnum, @(x) x>0 & x==floor(x) );
                        
            parse(par, obj, qc, pc, varargin{:});
            qc  = par.Results.qc;
            pc  = par.Results.pc;
            periodic = par.Results.periodic;
            verbose = par.Results.verbose;            
            pnum = par.Results.pnum;
            
            if strcmp(obj.dtype, 'mp') && verbose
                if abs(mp(num2str(pc)) - pc) ~=0 || abs(mp(num2str(qc)) - qc) ~=0
                %if strcmp(class(qc), 'double') || strcmp(class(pc), 'double') 
                    st = dbstack;                   
                    funcname =st.name;
                    str = sprintf( "qc(%s):%s\npc(%s):%s", class(qc), num2str(mp(qc)), class(pc), num2str(mp(pc)) );
                    precisionWarning(str, funcname);
                end
            end
            
            basis = obj.basis;            
            
            if periodic
                v = csp(obj.q, obj.hbar, obj.domain(1,:), qc, pc, pnum);
                a = dot(v, v);
                obj = FundamentalState(obj.sysinfo, 'q', v / sqrt(a));
                if strcmp(basis, 'p')
                    obj = obj.q2p();
                end
            elseif ~periodic
                v = cs(obj.q, obj.hbar, qc, pc);
                a = dot(v, v);
                obj = FundamentalState(obj.sysinfo, 'q', v / sqrt(a));
                obj.periodic = false;
                if strcmp(basis, 'p')
                    obj = obj.q2p();
                end
            else
                error("periodic must be true/false");
            end
        end
                        
        
        function obj = delta(obj, xc, varargin)
            
            par = inputParser;
            addRequired(par, 'obj', @isobject);                        
            addRequired(par, 'xc',  @isreal);

            default = true;
            addOptional(par,'verbose', default, @islogical);
                        
            parse(par, obj, xc, varargin{:});
            xc  = par.Results.xc;
            verbose = par.Results.verbose;            
                        
            [~, ind] = min( abs(obj.x - xc) );
            y = zeros(obj.dim, 1, class(obj.domain) );
            y(ind) = 1;
            
            if xc ~= obj.x(ind) && verbose
                msg = {
                    sprintf('input argument \nxc = %.16f but the ''xc'' is set to ', xc)
                    sprintf('xc = %.16f (obj.y(ind) = 1 where index = %d) due to the discreteness', obj.x(ind), ind)
                    'If you want to be slient on this warning, input arguments add as'
                    sprintf('delta(xc, ''verbose'', false)')
                    };
                warning('%s\n', msg{:});
            end
            obj = FundamentalState(obj.sysinfo, obj.basis, y);            
        end

        
        function save_state()
        end
        
        function save_eigs()
        end
        
    end
end

function y = cs(x, hbar, qc, pc)
  arguments
      x {mustBeNumeric};      
      hbar {mustBeReal};
      qc {mustBeReal};
      pc {mustBeReal};
  end
  if isrow(x)
      x = x.';
  end  
  y = exp( -(x - qc) .^ 2 /(2 * hbar) + 1j*(x - qc) * pc / hbar);
end

function y = csp(x, hbar, qdomain, qc, pc, partition)
  arguments
      x {mustBeNumeric};
      hbar {mustBeReal};
      qdomain (1,2) {mustBeReal};
      qc {mustBeReal};
      pc {mustBeReal};
      partition {mustBeInteger} = 5; % mustbe odd number
  end
  
  if mod(partition, 2) == 0
      error("partition must be odd integer");
  end
  
  if isrow(x)
      x = x.';
  end
  
  dim = length(x);
  width = qdomain(1,2) - qdomain(1,1);
  blk = fix(partition/2);
  
  qmin = qdomain(1,1) - blk*width;
  qmax = qdomain(1,2) + blk*width;
  q = linspace(qmin, qmax, partition*dim);
  %q = q(1:end-1);
  
  vec = cs(q, hbar, qc, pc);
  vec = reshape(vec, dim, partition);
  y = zeros(dim, 1, class(x));
  for i = 1:partition
      y = y + vec(:,i);
  end
end

function y = abs2(x)
  y = abs( conj(x) .* x);
end 
   
function basisconsistency(v1, v2)
  if ~strcmp(v1.basis, v2.basis)
      error("basis must be same");
  end
end

function precisionWarning(args, funcname)
format longG;
if length(args) == 1
    argstr = sprintf("%s", args);
else
    argstr = sprintf("%s", args{:});
end

msg = {"The computation precision is set as mp, but your the input arguments"
    sprintf("%s", argstr)
    "are double. This is not safe to keep preciison."
    "If you want to be slient on this warning, input arguments add as"
    sprintf("%s(..., verbose, false)", funcname)
    };
warning(sprintf("%s\n",msg{:}))
end


