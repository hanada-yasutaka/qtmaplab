classdef SystemInfo < matlab.mixin.SetGet
    % SystemInfo the domain of the system
    
    properties (SetAccess = protected)
        dim {mustBeInteger, mustBePositive} % Hilbert space dimension
        domain(2,2) {mustBeReal} = [0,0; 0,0] % domain of phase space: must be [qmin, qmax; pmin, pmax]
        hbar {mustBeReal} % effective Planck constant / (2*pi)
        q (:,1) {mustBeNumeric} % a column vector for the position coordinate (q-direction).
        p (:,1) {mustBeNumeric} % a column vector for the momentum coordinate (p-direction).
        eps {mustBeNumeric} % machine epsilon
        pi {mustBeNumeric} % pi
        twopi {mustBeNumeric} % 2*pi
        stype = '' % system type such as  'Hamiltonian', 'Unitary', ...
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
            obj.dtype = class(domain);
            
            if exist('stype', 'var')
                obj.stype = stype;
            end
            
            w = (domain(1,2) - domain(1,1) ) * (domain(2,2) - domain(2,1));
            planck = w/obj.dim;
            obj.hbar = planck/twopi;
            obj.eps = eps(obj.hbar);
            
            if mod(obj.dim, 2) == 0
                dq = (domain(1,2) - domain(1,1))/obj.dim;
                dp = (domain(2,2) - domain(2,1))/obj.dim;
                q = domain(1,1): dq: domain(1,2) - dq;
                p = domain(2,1): dp: domain(2,2) - dp;               
            else
                dq = (domain(1,2) - domain(1,1))/obj.dim;
                dp = (domain(2,2) - domain(2,1))/obj.dim;
                q = domain(1,1) + dq/2: dq: domain(1,2) - dq/2;
                p = domain(2,1) + dp/2: dp: domain(2,2) - dp/2;                
            end
            obj.q = q.';
            obj.p = p.';
        end
    end
end

