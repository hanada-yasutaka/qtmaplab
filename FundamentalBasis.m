classdef FundamentalBasis < handle
    properties
        dim {mustBeInteger, mustBePositive};
        domain(2,2) {mustBeReal} = [0,0; 0,0];
        hbar {mustBeReal};
        q {mustBeNumeric};
        p {mustBeNumeric};       
        rep {mustBeMember(rep,{'q','p'})} = 'p';
        eps {mustBeNumeric};
    end
    
    methods
        function obj = FundamentalBasis(dim, domain, rep)
            assert(domain(1,2) > domain(1,1) );
            assert(domain(2,2) > domain(2,1) );
            
            obj.dim = dim;
            obj.domain = domain;
            obj.rep = rep;
            
            w = (domain(1,2) - domain(1,1) ) * (domain(2,2) - domain(2,1));
            planck = w/obj.dim;
            eps(planck)
            class(planck);
            
            if class(planck) == "double"
                twopi = 2*pi;
            elseif class(planck) == "mp"
                twopi = mp('2*pi');
            else
                error("something wrong.");
            end
            
            obj.hbar = planck/twopi;
            obj.eps = eps(obj.hbar);
            
            q = linspace(domain(1,1), domain(1,2), obj.dim+1);
            p = linspace(domain(2,1), domain(2,2), obj.dim+1);            
            obj.q = q(1:end-1);
            obj.p = p(1:end-1);
        end
        
        function mat = matT(obj, funcT)
            % return <x|T(p)|x> where x = 'q' or 'p'
            if obj.rep == 'p'
                mat = obj.pBaseMatT(funcT);
            elseif obj.rep == 'q'
                mat = obj.qBaseMatT(funcT);
            else
                error("rep must be 'q' or 'p'");
            end
        end
        
        function mat = matV(obj, funcV)
            % return <x|V(p)|x> where x = 'q' or 'p'
            if obj.rep == 'p'
                mat = obj.pBaseMatV(funcV);
            elseif obj.rep == 'q'
                mat = obj.qBaseMatV(funcV);
            else
                error("rep must be 'q' or 'p'");
            end
        end
            
  
        function mat = pBaseMatT(obj, funcT)
            % retun <p'|T(p)|p>  = \delta_{p,p'}T(p)
            x = obj.p;
            mat = diag(funcT(x));
        end
        
        function mat = pBaseMatV(obj, funcV)
            % retun <p'|V(p)|p>  = \sum_q <p'|V(q)|q><q|p>
            x = obj.q;
            mat = zeros(length(x), class(x));
            for i = 1:length(x)
                pvec = zeros(1, length(x), class(x));
                pvec(i) = 1;
                qvec = ifft(pvec);
                if x(1)*x(end) < 0
                    qvec = fftshift(funcV(x)) .* qvec;
                else
                    qvec = funcV(x) .* qvec;
                end
                mat(:,i) = fft(qvec);
            end
        end
        
        function mat = qBaseMatV(obj, funcV)
            % retun <q'|V(p)|q>  = \delta_{q,q'}V(q)
            x = obj.q;
            mat = diag(funcV(x));
        end
        
        function mat = qBaseMatT(obj, funcT)
            % retun <q'|T(p)|q>  = \sum_p <q'|T(p)|p><p|q>            
            x = obj.p;
            mat = zeros(length(x), class(x));
            for i = 1:length(x)
                qvec = zeros(1, length(x), class(x));
                qvec(i) = 1;
                pvec = fft(qvec);
                if x(1)*x(end) < 0
                    pvec = fftshift(funcT(x)) .* pvec;
                else
                    pvec = funcT(x) .* pvec;
                end
                mat(:,i) = ifft(pvec);
            end
        end
     
    end
end

