classdef SplitUnitary  < matlab.mixin.SetGet
    %SplitUnitary このクラスの概要をここに記述
    %   詳細説明をここに記述
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
        function obj = SplitUnitary(dim, domain, basis)
            obj.system = SystemInfo(dim, domain, basis, class(obj));
            obj.state = FundamentalState(obj.system);
            

            props = properties(obj.system);
            for i = 1:length(props)
                p = get(obj.system, props{i});
                if isprop(obj, props{i})
                    set(obj, props{i}, p);
                end
            end
                        
            %obj = scaleinfo(obj, dim, domain);
            %obj.basis = basis;
        end
        
        function y = nullstate(obj)
            y = FundamentalState(obj.system);
        end
        
        function y = systeminfo(obj)
            y = obj.system
        end
        
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

