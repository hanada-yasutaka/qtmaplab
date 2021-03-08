classdef SystemInfo < matlab.mixin.SetGet
    % SystemInfo the domain of the system
    
    properties (SetAccess = protected)
        dim {mustBeInteger, mustBePositive} % HilbertSpace
        domain(2,2) {mustBeReal} = [0,0; 0,0] % domain of phase space: must be [qmin, qmax; pmin, pmax]
        hbar {mustBeReal} % effective Planck constant / (2*pi)
        %basis {mustBeMember(basis,{'q','p'})} = 'p' % basis (representation) of wavefunction         
        q {mustBeNumeric} % position coordinate
        p {mustBeNumeric} % momentum coordinate
        eps {mustBeNumeric} % machine epsilon
        pi {mustBeNumeric}
        twopi {mustBeNumeric}
        stype = '' % system type such as  'Hamiltonian', 'Unitary' ...
        dtype {mustBeMember(dtype,{'mp','double'})} = 'double' % data type: 'double' or 'mp'
        periodic {mustBeNumericOrLogical} = true % coordinate is periodic or not
        isfftshift = [false false]
    end
    
    methods
        function obj = SystemInfo(dim, domain, stype)
            
            assert(domain(1,2) > domain(1,1), 'domain(1,1) >= domain(1,2)');
            assert(domain(2,2) > domain(2,1), 'domain(2,1) >= domain(2,2)');
            
            if (domain(1,1) == 0) || (domain(1,2) == 0)
                obj.isfftshift(1) = false;
            elseif abs( domain(1,1)) == abs(domain(1,2) )
                obj.isfftshift(1) = true;
            else
                error('non-symmetric domain is not supported yet');
            end
            
            if (domain(2,1) == 0) || (domain(2,2) == 0)
                obj.isfftshift(2) = false;
            elseif abs(domain(2,1)) == abs(domain(2,2))
                obj.isfftshift(2) = true;
            else
                error('non-symmetric domain is not supported yet');
            end
            

            
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
            %obj.basis = basis;
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
        
        
    end
end

