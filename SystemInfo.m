classdef SystemInfo < matlab.mixin.SetGet
    %SCALEINFO このクラスの概要をここに記述
    % define a domain of system
    
    properties (SetAccess = protected)
        dim {mustBeInteger, mustBePositive};
        domain(2,2) {mustBeReal} = [0,0; 0,0];
        hbar {mustBeReal};
        q {mustBeNumeric};
        p {mustBeNumeric};
        eps {mustBeNumeric};
        pi {mustBeNumeric};
        twopi {mustBeNumeric};
        basis {mustBeMember(basis,{'q','p'})} = 'p';
        stype = '';
        dtype {mustBeMember(dtype,{'mp','double'})} = 'double';
    end
    
    methods
        function obj = SystemInfo(dim, domain, basis, stype)
            
            assert(domain(1,2) > domain(1,1), 'domain(1,1) >= domain(1,2)');
            assert(domain(2,2) > domain(2,1), 'domain(2,1) >= domain(2,2)');
            
            if class(domain) == "double"
                twopi = 2*pi;
                obj.pi = pi;
            elseif class(domain) == "mp"
                twopi = mp('2*pi');
                obj.pi = mp('pi');
            else
                str = sprintf('class(domain)=%s\nclass of domain must be "double" or "mp"\n', class(domain));
                error(str);
            end
            
            obj.dim = dim;
            obj.domain = domain;
            obj.twopi = twopi;
            obj.basis = basis;
            obj.dtype = class(domain);
            
            if exist('stype', 'var')
                obj.stype = stype;
            end
            
            w = (domain(1,2) - domain(1,1) ) * (domain(2,2) - domain(2,1));
            planck = w/obj.dim;
            obj.hbar = planck/twopi;
            obj.eps = eps(obj.hbar);

            q = linspace(domain(1,1), domain(1,2), obj.dim+1);
            p = linspace(domain(2,1), domain(2,2), obj.dim+1);
            obj.q = q(1:end-1);
            obj.p = p(1:end-1);
        end
        
        
        function bool = isfftshift(basis)
            if basis == 'q'
                bool = true;
            else
                bool = true;
            end
            
        end
    end
end

