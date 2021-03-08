classdef SplitUnitary  < matlab.mixin.SetGet & SystemInfo
    %SplitUnitary このクラスの概要をここに記述
    %   詳細説明をここに記述
    properties (SetAccess = protected)
        basis {mustBeMember(basis,{'q','p'})} = 'p';
        tau {mustBeNumeric} = 1;
        funcT;
        funcV;
    end
    
    methods
        function obj = SplitUnitary(dim, domain, basis)
            obj@SystemInfo(dim, domain, 'SplitUnitary');
            obj.basis = basis;
        end
        
        function op = op_expT(obj, funcT, tau, isvec)
            arguments
                obj
                funcT
                tau = 1
                isvec = true
            end

            if strcmp(obj.basis, 'q') & obj.isfftshift(2)
                x = fftshift(obj.p);
            else
                x = obj.p;
            end

            
            if isvec
                op = @(invec) exp( -1j * funcT(x) * tau /obj.hbar) .* invec;
            else
                op = @(invec) transpose( exp( -1j * funcT(x) * tau /obj.hbar) )  .* invec;
            end
        end
        
        function op = op_expV(obj, funcV, tau, isvec)
            arguments
                obj
                funcV
                tau = 1
                isvec = true;
            end
            
            x = obj.q;
            
            if isvec
                op = @(invec) exp( -1j * funcV(x) * tau / obj.hbar) .* invec;
            else
                op = @(invec) transpose(exp( -1j * funcV(x) * tau / obj.hbar) ) .* invec;
            end
        end
                
        function mat = mat_expTV(obj, funcT, funcV, tau)
            arguments
                obj
                funcT
                funcV
                tau = 1
            end
                        
            expT = op_expT(obj, funcT, tau, false);
            expV = op_expV(obj, funcV, tau, false);
            
            iden = diag( ones(1, obj.dim, class(obj.domain) ) );
            if strcmp(obj.basis, 'p')                              
                iden = ifft(iden);
            end
            
            pvecs = fft( expV(iden ) );
            mat  = expT(pvecs);
            
            if strcmp(obj.basis, 'q')
                mat = fft(mat);
            end                        
        end
        
        function mat = mat_expVT(obj, funcT, funcV, tau)
            arguments
                obj
                funcT
                funcV
                tau = 1
            end
                        
            expT = op_expT(obj, funcT, tau, false);
            expV = op_expV(obj, funcV, tau, false);
            
            iden = diag( ones(1, obj.dim, class(obj.domain) ) );
            if strcmp(obj.basis, 'q')                              
                iden = fft(iden);
            end
            
            pvecs = ifft( expT(iden ) );
            mat  = expV(pvecs);
            
            if strcmp(obj.basis, 'p')
                mat = fft(mat);
            end    
        end
        
        function mat = mat_expVTV(obj, funcT, funcV, tau)
            arguments
                obj
                funcT
                funcV
                tau = 1
            end
                        
            expT = op_expT(obj, funcT, tau, false);
            expV = op_expV(obj, funcV, tau/2, false);
            
            iden = diag( ones(1, obj.dim, class(obj.domain) ) );
            if strcmp(obj.basis, 'p')
                iden = ifft(iden);
            end
            
            pvecs =  fft( expV(iden ) );
            qvecs = ifft( expT(pvecs) );
            mat   =  fft( expV(qvecs) );
            
            if strcmp(obj.basis, 'q')
                mat = ifft(mat);
            end                        

        end
        
        function mat = mat_expTVT(obj, funcT, funcV, tau)
            arguments
                obj
                funcT
                funcV
                tau = 1
            end
                        
            expT = op_expT(obj, funcT, tau/2, false);
            expV = op_expV(obj, funcV, tau, false);
            
            iden = diag( ones(1, obj.dim, class(obj.domain) ) );
            if strcmp(obj.basis, 'q')
                iden = fft(iden);
            end
            
            qvecs =  ifft( expT(iden ) );
            pvecs =  fft ( expV(qvecs) );
            mat   =  ifft( expT(pvecs) );
            
            if strcmp(obj.basis, 'p')
                mat = fft(mat);
            end                        

        end
        
        
        
       
            
        
%         function mat = expTV(obj, funcT, funcV, tau)
%             % return <q'|exp(-i/hbar *funcT(p)*tau) exp(-i/hbar*funcV(q)*tau)|q>            
%             arguments
%                 obj
%                 funcT
%                 funcV
%                 tau = 1
%             end
%                         
%             if obj.isfftshift(2)
%                 p = fftshift(obj.p);
%             else
%                 p = obj.p;
%             end
%             
%             if strcmp(obj.basis, 'q')
%                 qvecs = diag( ones(1, obj.dim, class(obj.domain) ) );
%             else
%                 pvecs = diag( ones(1, obj.dim, class(obj.domain) ) );
%                 %qvecs = ifft(pvecs)/sqrt(obj.dim);
%                 qvecs = pvecs;
%             end
%             
%             qvecs = exp( -1j * funcV(obj.q) * tau / obj.hbar ) .* qvecs;
%             mat = transpose( exp( -1j * funcT(p) * tau / obj.hbar ) ) .* fft( qvecs ) / sqrt(obj.dim);
%             
%             if strcmp(obj.basis, 'q')
%                 mat = ifft(mat) / sqrt(obj.dim);
%             end
%         end
%         
%         function mat = expVT(obj, funcT, funcV, tau)
%             arguments
%                 obj
%                 funcT
%                 funcV
%                 tau = 1
%             end
%                         
%             if obj.isfftshift(2)
%                 q = fftshift(obj.q);
%             else
%                 q = obj.q;
%             end
%             
%             idmat = diag( ones(1, obj.dim, class(obj.domain) ) );
%             idmat = fft(idmat);
%             pvecs = exp( -1j * funcT(obj.p) * tau / obj.hbar ) .* idmat;
%             mat = transpose( exp( -1j * funcV(obj.q) * tau / obj.hbar ) ) .* ifft( pvecs ) / sqrt(obj.dim);
%             
%             if strcmp(obj.basis, 'p')
%                 mat = fft(mat) / sqrt(obj.dim);
%             end
%         end
        
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
        
        function index = sort(obj, targetmat, refmat)
            mat = conj(targetmat).' * refmat;
            mat2 = conj(mat) .* mat;
            
            qnumber = 1:obj.dim;
            index = [];
            for i=1:obj.dim
                [ovl, ind] = max(mat2(:,i));
                index = [index, qnumber(ind)];
                qnumber(ind) = [];
                mat2(ind,:) = [];
            end

            if length(unique(index)) ~= obj.dim
                error("duplicate sort index ")
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

