classdef SplitUnitary  < matlab.mixin.SetGet & SystemInfo
    %SplitUnitary このクラスの概要をここに記述
    %   詳細説明をここに記述
    properties (SetAccess = protected)
        %dim {mustBeInteger, mustBePositive};
        %domain(2,2) {mustBeReal} = [0,0; 0,0];
        %hbar {mustBeReal};
        %q {mustBeNumeric};
        %p {mustBeNumeric};       
        basis {mustBeMember(basis,{'q','p'})} = 'p';
        %eps {mustBeNumeric};
        tau {mustBeNumeric} = 1;
        funcT;
        funcV;
        %system;
        %state;
    end
    
    methods
        function obj = SplitUnitary(dim, domain, basis)
            obj@SystemInfo(dim, domain, 'SplitUnitary');
            obj.basis = basis;
        end
        
        function qmat = expTV(obj, funcT, funcV, tau)
            % return <q'|exp(-i/hbar *funcT(p)*tau) exp(-i/hbar*funcV(q)*tau)|q>            
            arguments
                obj
                funcT
                funcV
                tau = 1
            end
                        
            if obj.isfftshift(2)
                p = fftshift(obj.p);
            else
                p = obj.p;
            end
            
            if strcmp(obj.basis, 'q')
                qvecs = diag( ones(1, obj.dim, class(obj.domain) ) );
            elseif
                pvecs = diag( ones(1, obj.dim, class(obj.domain) ) );
                qvecs = ifft(pvecs)/sqrt(obj.dim);
            end
            
            qvecs = exp( -1j * funcV(obj.q) * tau / obj.hbar ) .* qvecs;
            pvecs = transpose( exp( -1j * funcT(p) * tau / obj.hbar ) ) .* fft( qvecs );
            qmat = ifft(pvecs);
        end
        
        function expVT(obj, funcT, funcV)
            
        end
        
        %function y = nullstate(obj)
        %    y = FundamentalState(obj.system);
        %end
        
        %function y = systeminfo(obj)
        %    y = obj.system
        %end
        
        function op = evolutionOP(obj)
            op = @(x) EvolveTV(x);
        end
        
        function setop(obj, funcT, funcV, tau)
            if exist('tau', 'var')
                obj.tau = tau;
            end           
            obj.funcT = funcT; % = set(obj, 'funcT', funcT);
            obj.funcV = funcV; % = set(obj, 'funcV', funcV);
        end
        
        function op = getop(obj)
            op = @(x) evolveTV(obj, x)
        end

        
        function outvec = evolveTV(obj, invec)            
            arguments
                obj
                invec (:,1) {mustBeNumeric}
                %funcT % mustbefunction @funcT
                %funcV % mustbefunction @funcV
            end
            
            if obj.basis == 'p'
                invec = ifft(invec);
            end
            
            q = transpose(obj.q);
            p = transpose(obj.p);
            
            size(q)
            size(invec)
            isfftshift = true;
           
            if isfftshift
                V = @(x) obj.funcV( fftshift(x) );
                T = @(x) obj.funcT( fftshift(x) );
            else
                V = @(x) obj.funcV(x);
                T = @(x) obj.funcT(x);               
            end
            
            vec = exp( -1j * V(q) / obj.hbar ) .* invec;
            vec = fft(vec);
            vec = exp( -1j * T(p) / obj.hbar) .* vec;
            outvec = ifft(vec);

            if obj.basis == 'p'
                outvec = fft(outvec);
            end            
            
        end
        
        function outvec = EvolveVT(obj, invec, funcT, funcV)
            arguments
                obj
                invec (1,:) {mustBeNumeric}
                funcT % mustbefunction @funcT
                funcV % mustbefunction @funcV
            end
            
            if obj.basis == 'q'
                invec = q2p(invec);
            end            
            
            vec = exp( -1j .* funcT(obj.q) / obj.hbar ) .* invec;
            vec = fft(vec);
            vec = exp( -1j .* funcV(obj.p) / obj.hbar) .* vec;
            outvec = ifft(vec);
            
            if obj.basis == 'q'
                invec = q2p(invec);
            end                                    
        end
        
        function a = AAA(obj)
            a = 1
        end
        
        
        function y = test(obj, func)
            arguments
                obj
                func
            end
            disp(isa(@func, 'function_handle'))
            y = 1;
        end
        
        
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 このメソッドの概要をここに記述
            %   詳細説明をここに記述
            outputArg = obj.Property1 + inputArg;
        end
    end
end

