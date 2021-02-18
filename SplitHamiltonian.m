classdef SplitHamiltonian < matlab.mixin.SetGet & SystemInfo
    properties (SetAccess = protected)
        basis {mustBeMember(basis,{'q','p'})} = 'p';
        tau {mustBeNumeric} = 1;
        funcT;
        funcV;
    end
    
    methods
        function obj = SplitHamiltonian(dim, domain, basis)
            obj@SystemInfo(dim, domain, 'SplitHamiltonian');
            obj.basis = basis;
        end    
        
        function mat = matT(obj, funcT)
            % return <x|T(p)|x> where x = 'q' or 'p'
            if obj.basis == 'p'
                mat = obj.pBaseMatT(funcT);
            elseif obj.basis == 'q'
                mat = obj.qBaseMatT(funcT);
            else
                error("representation error");
            end
        end
        
        function mat = matV(obj, funcV)
            % return <x|V(p)|x> where x = 'q' or 'p'
            if obj.basis == 'p'
                mat = obj.pBaseMatV(funcV);
            elseif obj.basis == 'q'
                mat = obj.qBaseMatV(funcV);
            else
                error("representation error...");                
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
        
        function mat = qBaseMatT_old(obj, funcT)
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
        
        function mat = qBaseMatT(obj, funcT)
            % retun <q'|T(p)|q>  = \sum_p <q'|T(p)|p><p|q>            
            x = obj.p;
            if x(1)*x(end)<0
                f = @(x) fftshift( funcT(x) );
            else
                f = @(x) funcT(x);
            end
            qvecs = diag( ones(1, obj.dim, class(x) )); 
            qvecs = transpose( f(x) ) .* fft(qvecs);
            mat = ifft(qvecs);
        end
        
        
        function state = nullstate(obj)
            system = SystemInfo(obj.dim, obj.domain, class(obj) );            
            y = zeros(1, obj.dim, class(obj.domain));
            state = FundamentalState(system, obj.basis, y);            
        end
        
        function state = vec2state(obj, data)
            system = SystemInfo(obj.dim, obj.domain, class(obj) );            
            state = FundamentalState(system, obj.basis, data);                        
        end
        
        function states = eigs2states(obj, evals, evecs)
            states = [];
            dim = obj.dim;
            system = SystemInfo(obj.dim, obj.domain, class(obj) );
            
            if isequal( size(evals), [dim dim] )
                evals = diag(evals);
            end
            
            for i = 1:dim
                stat = FundamentalState(system, obj.basis, evecs(:, i), evals(i) );
                states = [states; stat];
            end
        end
        
        
        
    end
end


