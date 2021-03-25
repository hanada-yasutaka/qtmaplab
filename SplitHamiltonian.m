classdef SplitHamiltonian < matlab.mixin.SetGet & SystemInfo
    % SplitHamiltonian provides the matrix representation for a Hamiltonian such a split (separate) form
    % 
    % $H(p,q) = T(p) + V(q)$ 
    %
    % by using fundamental basis ('q' or 'p').
    %
    % Example:
    %
    %     dim = 50;
    %     domain = [-pi pi; -pi pi];         % or mp('[-pi pi; -pi pi'])
    %     basis = 'p';                       % or 'p'
    %     sH = SplitHamiltonian(dim, domain, basis);
    %     matT = sH.matT(@(x) x.^2/2);       %<p'|T(p)|p> where T(p) = p^2/2
    %     matV = sH.matV(@(x) cos(x));       %<p'|V(p)|p> where V(q) = cos(q)
    %     matH = matT + matV;                %<p'|H|p>
    %
    % or 
    %
    %     T = @(x) x.^2 /2;
    %     matT = sH.matT(T);
    %
    % or 
    %
    %     matT = sH.matT(@T);
    %     ...
    %     function y = T(x)
    %         y = x.^2/2;
    %     end
    %
    
    properties (SetAccess = protected)
        basis {mustBeMember(basis,{'q','p'})} = 'p'; % representation basis which must be 'q' or 'p'
    end
    
    methods
        function [obj, state] = SplitHamiltonian(dim, domain, basis)
            obj@SystemInfo(dim, domain, 'SplitHamiltonian');
            obj.basis = basis;
            if strcmp(basis, 'q')
                fprintf('Worning: basis "q" has not fully tested yet')
            end
            state = FundamentalState(obj, obj.basis);            
        end
        
        function state = vec2state(obj, vec, basis)
            % create FundamentalState from vector array
            if isrow(vec)
                vec = vec.';
            end
            
            if ~exist('basis', 'var')
                basis = obj.basis;
            end

            sysinfo = SystemInfo(obj.dim, obj.domain, 'SplitHamiltonian');
            state = FundamentalState(sysinfo, basis, vec);
        end
        
        
        function mat = matT(obj, funcT)
            % return <x|T(p)|x> where x is obj.basis ('q' or 'p')
            msg = {'argument must be function_handle'
                'usage:'
                '    matT(@(x) x)'
                'or'
                '    T=@(x) x'
                '    matT(T),'
                'or'
                '    matT(@T);'
                '    ...'
                '    function y = T(x)'
                '        y = x'
                '    end'
                };
            %fprintf('%s\n', msg{:})
            %assert( strcmp(class(funcT), 'function_handle'), sprintf('%s\n', msg{:}) );
            if obj.basis == 'p'
                mat = obj.pBaseMatT(funcT);
            elseif obj.basis == 'q'
                mat = obj.qBaseMatT(funcT);
            else
                error("representation error");
            end
        end
        
        function mat = matV(obj, funcV)
            % return <x|V(p)|x> where x is obj.basis ('q' or 'p')
            if obj.basis == 'p'
                mat = obj.pBaseMatV(funcV);
            elseif obj.basis == 'q'
                mat = obj.qBaseMatV(funcV);
            else
                error("representation error...");                
            end
        end                 
        
        %function state = nullstate(obj)
        %    system = SystemInfo(obj.dim, obj.domain, class(obj) );            
        %    y = zeros(1, obj.dim, class(obj.domain));
        %    state = FundamentalState(system, obj.basis, y);            
        %end                     
    end
    
    methods(Access = private)
        
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
            %pvecs = transpose( funcV(x) ) .* ifft(pvecs);
            pvecs = funcV(x) .* ifft(pvecs);            
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
            %qvecs = transpose( funcT(x) ) .* fft(qvecs);
            qvecs = funcT(x)  .* fft(qvecs);
            mat = ifft(qvecs);
        end
    end
end


