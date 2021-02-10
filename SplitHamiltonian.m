classdef SplitHamiltonian < matlab.mixin.SetGet
    properties (SetAccess = protected)
        dim {mustBeInteger, mustBePositive};
        domain(2,2) {mustBeReal} = [0,0; 0,0];
        hbar {mustBeReal};
        q {mustBeNumeric};
        p {mustBeNumeric};       
        basis {mustBeMember(basis,{'q','p'})} = 'p';
        eps {mustBeNumeric};
        tau {mustBeNumeric} = 1;
        funcT;
        funcV;
        system;
        state;
    end
    
    methods
        function obj = SplitHamiltonian(dim, domain, basis)
            obj.system = SystemInfo(dim, domain, basis, class(obj));
            obj.state = FundamentalState(obj.system);

            props = properties(obj.system);
            for i = 1:length(props)
                p = get(obj.system, props{i});
                if isprop(obj, props{i})
                    set(obj, props{i}, p);
                end
            end
        end    
%     properties
%         dim {mustBeInteger, mustBePositive};
%         domain(2,2) {mustBeReal} = [0,0; 0,0];
%         hbar {mustBeReal};
%         q {mustBeNumeric};
%         p {mustBeNumeric};       
%         basis {mustBeMember(basis,{'q','p'})} = 'p';
%         eps {mustBeNumeric};
%     end
%     
%     methods
%         function obj = SplitHamiltonian(dim, domain, basis)
%             obj.basis = basis;
%             obj = scaleinfo(obj, dim, domain);
%         end
        
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
            data = zeros(1, obj.dim, class(obj.domain));
            state = FundamentalState(obj.dim, obj.domain, obj.basis, data)
        end
        
        function states = eigs2states(obj, evals, evecs)
            states = [];
            dim = obj.dim;
            if isequal( size(evals), [dim dim] )
                evals = diag(evals);
            end
            
            for i = 1:dim
                stat = FundamentalState(obj.system, evecs(:, i), evals(i));
                states = [states; stat];
            end
        end
        
        
        
    end
end


