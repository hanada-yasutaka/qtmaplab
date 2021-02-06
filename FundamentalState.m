classdef FundamentalState < matlab.mixin.SetGet
    properties (SetAccess = protected)
        dim {mustBeInteger, mustBePositive};
        domain(2,2) {mustBeReal} = [0,0; 0,0];
        hbar {mustBeReal};
        q {mustBeNumeric};
        p {mustBeNumeric};  
        basis {mustBeMember(basis,{'q','p'})} = 'q';
        eps {mustBeNumeric};
        twopi {mustBeNumeric};
        system;
        eigenvalue {mustBeNumeric};
        y (:,1) {mustBeNumeric} = [];
        periodic {mustBeNumericOrLogical} = true;
    end
    
    properties (Access = private)
        isshift;
    end
    
    methods
        function obj = FundamentalState(system, vec, eigenvalue)
            if ~isa(system, 'SystemInfo')
                error("scl must be SystemInfo class")
            end
            %obj@SystemInfo(scl.dim, scl.domain, scl.basis)     
            obj.system = system;
            
            props = properties(system);
            for i = 1:length(props)
                p = get(system, props{i});
                if isprop(obj, props{i})
                    set(obj, props{i}, p);
                end
            end
            
            if exist('vec', 'var')
                obj.y = vec;
            end
            if exist('eigenvalue', 'var')
                obj.eigenvalue = eigenvalue;
            end
        end
        
        function y = toarray(obj)
            obj.y
        end
                
        function y = q2p(obj)
            y = 1;
        end
        
        function y = p2q(obj)
            y = 1;
        end
        
 
       function obj = double(obj)
            scl = SystemInfo(obj.dim, double(obj.domain), obj.basis);
            obj = FundamentalState(scl, double(obj.y));            
        end
         
        function y = zeros(obj)
            data = zeros(1, obj.dim, class(obj.domain));
            y = FundamentalState(obj.system, data);
        end
        
        function y = abs2(obj)
            y = abs( obj.y .* conj(obj.y) );
        end
        
        function y = norm(obj)
           y = norm( obj.y );
        end
        
        function z = inner(x, y)
            % return < x | y >
            xx = x.data;
            if strcmp( class(y), 'FundamentalState')
                y = y.data;
            elseif size(xx) ~= size(yy)
                error("size must be same");
            end
            % dot(x, y) = abs2( sum( conj(x) .* y) )
            z = dot(xx , y);
        end
        
        function obj = normalize(obj)
            x = obj.y;
            a = dot(x, x);
            obj = FundamentalState(obj.system, x/sqrt(a) );
        end
        
        function obj = coherent(obj, qc, pc, periodic)
            arguments
                obj;
                qc {mustBeReal};
                pc {mustBeReal};
                periodic {mustBeNumericOrLogical} = true;
            end
            
            if periodic
                y = csp(obj.q, obj.hbar, obj.domain(1,:), qc, pc, 5);
                a = dot(y, y);
                obj = FundamentalState(obj.system, y / sqrt(a) );
            elseif ~periodic
                y = cs(obj.q, obj.hbar, qc, pc);
                a = dot(y, y);
                obj = FundamentalState(obj.system, y / sqrt(a) );
                obj.periodic = false;
            else
                error("periodic must be true/false");
            end
            
        end
                
        function [X, Y, mat] = hsmrep(obj, num, range)
            arguments
                obj;
                num {mustBeInteger, mustBePositive} = 50;                
                range {mustBeReal} = [obj.domain(1,1) obj.domain(1,2); obj.domain(2,1), obj.domain(2,2)];
            end
            vec = obj.y;
            x = transpose(obj.q);
            qc = linspace(obj.domain(1,1), obj.domain(1,2), num);
            pc = linspace(obj.domain(2,1), obj.domain(2,2), num);
            [X,Y] = meshgrid(qc,pc);
            mat = zeros(num, class(vec));
            
            if obj.periodic
                csfun = @(qc, pc) csp(x, obj.hbar, obj.domain(1,:), qc, pc, 5);                
            else
                csfun = @(qc, pc)  cs(x, obj.hbar, qc, pc);
            end
            
            for i = 1:length(qc)
                for j = 1:length(pc)
                    c = transpose( csfun(qc(j), pc(i)) );
                    % dot(c, vec) = abs2( sum( conj(c) .* vec) )                    
                    mat(i,j) = abs2( dot(c, vec));

                end
            end
        end
        
        
        function obj = delta(obj)
            
            obj.y=1;
        end
        
        
        function save_state()
        end
        
        function save_eigs()
        end
        
       %% operator overload
        
        %function display(obj)
        %    display(obj.y);
        %end
        
        %function disp(obj)
        %    disp(obj.y);
        %end
            
       
        %function y = plus(obj1, obj2)
        %    isstate(obj1, obj2);
        %    data = obj1.data + obj2.data;
        %    y = FundamentalState(obj1.dim, obj1.domain, obj1.basis, data);
        %end
        %
        %function isstate(obj1, obj2)
        %    if ~isa(obj1, 'FundamentalState')
        %        error(sprintf('object must be FundamentalState' ));
        %    elseif exist('obj2', 'var') && ~isa(obj2, 'FundamentalState')
        %        error(sprintf('object must be FundamentalState', getobjname(obj1) ));
        %    end
        %end
        %%%        
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
   
    

