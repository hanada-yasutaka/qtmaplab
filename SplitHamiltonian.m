classdef SplitHamiltonian < matlab.mixin.SetGet & SystemInfo
    properties (SetAccess = protected)
        basis {mustBeMember(basis,{'q','p'})} = 'p';
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

            %if obj.isfftshift(1)
            %    x = fftshift(obj.q);
                %fv = fftshift( funcV(obj.q) );
            %else
            x = obj.q;
                %fv = funcV(obj.q);
            %end
            
            mat = zeros( obj.dim, class(obj.domain));
            pvecs = diag( ones(1, obj.dim, class(obj.domain) ));
            %pvecs = transpose( fv ) .* ifft(pvecs);
            pvecs = transpose( funcV(x) ) .* ifft(pvecs);            
            mat = fft(pvecs);
        end
        
        function mat = qBaseMatV(obj, funcV)
            % retun <q'|V(p)|q>  = \delta_{q,q'}V(q)
            mat = diag( funcV(obj.q) );
        end
        
        function mat = qBaseMatT(obj, funcT)
            % retun <q'|T(p)|q>  = \sum_p <q'|T(p)|p><p|q>                        
            if obj.isfftshift(2)
                %fv = fftshift( funcT(obj.p) );
                x = fftshift(obj.p);
            else
                x = obj.p;
                %fv = funcT(obj.p);
            end
            qvecs = diag( ones(1, obj.dim, class(obj.domain) )); 
            qvecs = transpose( funcT(x) ) .* fft(qvecs);
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
        
%         function states = eigs2states(obj, evals, evecs)
%             states = [];
%             dim = obj.dim;
%             system = SystemInfo(obj.dim, obj.domain, class(obj) );
%             
%             if isequal( size(evals), [dim dim] )
%                 evals = diag(evals);
%             end
%             
%             for i = 1:dim
%                 stat = FundamentalState(system, obj.basis, evecs(:, i), evals(i) );
%                 states = [states; stat];
%             end
%         end
        
        
        
    end
end


