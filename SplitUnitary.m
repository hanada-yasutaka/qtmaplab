classdef SplitUnitary  < matlab.mixin.SetGet & SystemInfo
    % SplitUnitary provides a solver for the eigenvalue/time-evolution problem of the exponential split operator, e.g,
    % $U = \exp(-i*T(p)*\tau/\hbar) \exp(-i*V(p)*\tau/\hbar)$,
    % $U = \exp(-i*V(p)*\tau/\hbar) \exp(-i*T(p)*\tau/\hbar)$,
    % and
    % higher order symplectic integrator construction.
    %
    properties (SetAccess = protected)
        basis {mustBeMember(basis,{'p', 'q'})} = 'p'; % representation basis which must be 'q' or 'p'
        tau {mustBeNumeric} = 1; % time step of the quatnum map.
%        funcT;
%        funcV;
    end
    
    methods
        function [obj, state] = SplitUnitary(dim, domain, basis)
            obj@SystemInfo(dim, domain, 'SplitUnitary');
            obj.basis = basis;
            if strcmp(basis, 'q')
                fprintf('Worning: basis "q" has not fully tested yet')
            end
            state = FundamentalState(obj, obj.basis);
        end
        
        function state = vec2state(obj, vec, basis)
            if isrow(vec)
                vec = vec.';
            end
            
            if ~exist('basis', 'var')
                basis = obj.basis;
            end

            sysinfo = SystemInfo(obj.dim, obj.domain, 'SplitUnitary');
            state = FundamentalState(sysinfo, basis, vec);
        end
        
        function [op, info] = SIevolve(obj, funcT, funcV, varargin)
            % return evolution operator with n-th order symplectic integrator 
            %
            par = inputParser;
            
            addOptional(par, 'tau', 1, @isnumeric)
            addOptional(par, 'SIorder', 1, @(x) all(x>0 & x == floor(x)) );
            addOptional(par, 'OPorder', 'TV', @(x) ismember(x, ["TV" "VT" "TVT" "VTV"]) );
            
            parse(par, varargin{:} );            
            tau     = par.Results.tau;            
            SIorder = par.Results.SIorder;
            OPorder = par.Results.OPorder;
            
            if strcmp(OPorder, 'TVT')
                op = @(invec) obj.expTVTevolve(invec, funcT, funcV, tau);
            elseif strcmp(OPorder, 'VTV')
                op = @(invec) obj.expVTVevolve(invec, funcT, funcV, tau);
            elseif strcmp(OPorder, "TV")
                if SIorder == 1
                    op = @(invec) obj.expTVevolve(invec, funcT, funcV, tau);
                elseif SIorder == 2
                    op = @(invec) obj.expVTVevolve(invec, funcT, funcV, tau);
                else
                    error('siorder > 2 is not yes implement');
                end
            elseif strcmp(OPorder, "VT")
                if SIorder == 1
                    op = @(invec) obj.expVTevolve(invec, funcT, funcV, tau)
                elseif SIorder == 2
                    op = @(invec) obj.expTVTevolve(invec, funcT, funcV, tau)
                else
                    error('siorder > 2 is not yes implement');
                end                              
            else
                error('oporder must be ''TV'' or ''VT'' ');
            end
            info = '';
        end
        
        function [mat, info] = SImatrix(obj, funcT, funcV, varargin)
            % return <x'|U|x>
            %
            % where U is n-th order simplectic integrator 
            % and x = obj.basis
            par = inputParser;
            
            addOptional(par, 'tau', 1, @isnumeric)
            addOptional(par, 'SIorder', 1, @(x) all(x>0 & x == floor(x)) );
            addOptional(par, 'OPorder', 'TV', @(x) ismember(x, ["TV" "VT" "TVT" "VTV"]) );
            
            parse(par, varargin{:} );            
            tau     = par.Results.tau;            
            SIorder = par.Results.SIorder;
            OPorder = par.Results.OPorder;
            
            
            if strcmp(OPorder, 'TVT')
                
                mat = obj.expTVTmat(funcT, funcV, tau);
                
            elseif strcmp(OPorder, 'VTV')
                
                mat = obj.expVTVmat(funcT, funcV, tau);
                
            elseif strcmp(OPorder, 'TV')
                if SIorder == 1
                    
                    mat = obj.expTVmat(funcT, funcV, tau);
                    
                elseif SIorder == 2
                    
                    mat = obj.expVTVmat(funcT, funcV, tau);
                else
                    error('siorder > 2 is not yes implement');
                end
            elseif strcmp(OPorder, "VT")
                if SIorder == 1
                    mat = obj.expVTmat(funcT, funcV, tau);
                elseif SIorder == 2
                    mat = obj.expTVTmat(funcT, funcV, tau);
                else
                    error('siorder > 2 is not yes implement');
                end                              
            else
                error('oporder must be ''TV'' or ''VT'' ');
            end
            info = '';
        end
        
                
    end
    
    methods (Access = private)
        function op = op_expT(obj, funcT, tau)
            arguments
                obj
                funcT
                tau = 1
            end

            if strcmp(obj.basis, 'q') & obj.isfftshift(2)
                x = fftshift(obj.p);
            else
                x = obj.p;
            end
            %x = obj.p;

            op = @(vecs) exp( -1j * funcT(x) * tau /obj.hbar) .* vecs;
        end
        
        function op = op_expV(obj, funcV, tau)
            arguments
                obj
                funcV
                tau = 1
            end
            
            %if strcmp(obj.basis, 'p') & obj.isfftshift(1)
            %    x = ifftshift(obj.q);
            %else
            %    x = obj.q;
            %end            
            x = obj.q;

            op = @(vecs) exp( -1j * funcV(x) * tau /obj.hbar) .* vecs;
        end
        
                
        function mat = expTVmat(obj, funcT, funcV, tau)
            arguments
                obj
                funcT
                funcV
                tau = 1
            end
                        
            expT = op_expT(obj, funcT, tau);
            expV = op_expV(obj, funcV, tau);
            
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
        
        function mat = expVTmat(obj, funcT, funcV, tau)
            arguments
                obj
                funcT
                funcV
                tau = 1
            end
                        
            expT = op_expT(obj, funcT, tau);
            expV = op_expV(obj, funcV, tau);
            
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
        
        function mat = expVTVmat(obj, funcT, funcV, tau)
            arguments
                obj
                funcT
                funcV
                tau = 1
            end
                        
            expT = op_expT(obj, funcT, tau);
            expV = op_expV(obj, funcV, tau/2);
            
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
        
        function mat = expTVTmat(obj, funcT, funcV, tau)
            arguments
                obj
                funcT
                funcV
                tau = 1
            end
                        
            expT = op_expT(obj, funcT, tau/2);
            expV = op_expV(obj, funcV, tau);
            
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
        
        function outvec = expTVevolve(obj, invec, funcT, funcV, tau)
            arguments
                obj
                invec 
                funcT
                funcV
                tau {mustBeNumeric} = 1
            end
            
            expT = op_expT(obj, funcT, tau);
            expV = op_expV(obj, funcV, tau);
            
            if isa(invec, 'FundamentalState')
                dtypecmp(obj, invec.y);                                               
                                
                if strcmp(invec.basis, 'q')
                    qvec = invec.y;
                else
                    qvec = ifft(invec.y);
                end
                                
                pvec =  fft( expV(qvec) );
                qvec = ifft( expT(pvec) );
                                            
                if strcmp(invec.basis, 'q')
                    outvec = vec2state(obj, qvec, 'q');
                else
                    outvec = vec2state(obj, fft(qvec), 'p');
                end
            else
                dtypecmp(obj, invec);                
                
                pvec = fft(expV(invec));
                outvec = ifft(expT(pvec));
            end    
        end
                
        function outvec = expVTevolve(obj, invec, funcT, funcV, tau)
            arguments
                obj
                invec 
                funcT
                funcV
                tau {mustBeNumeric} = 1
            end
                                    
            expT = op_expT(obj, funcT, tau);
            expV = op_expV(obj, funcV, tau);
            
            if isa(invec, 'FundamentalState')                
                dtypecmp(obj, invec.y);
                
                if strcmp(invec.basis, 'p')
                    pvec = invec.y;
                else
                    pvec = fft(invec.y);
                end
                                
                qvec = ifft( expT(pvec) );
                pvec =  fft( expV(qvec) );
                                            
                if strcmp(invec.basis, 'p')
                    outvec = vec2state(obj, pvec, 'p');
                else
                    outvec = vec2state(obj, ifft(pvec), 'q');
                end
            else
                dtypecmp(obj, invec);                
                pvec   =  fft( expV(invec) );
                outvec = ifft( expT(pvec) );
            end    
        end
        
        function outvec = expVTVevolve(obj, invec, funcT, funcV, tau)
            arguments
                obj
                invec 
                funcT
                funcV
                tau {mustBeNumeric} = 1
            end
            
            expT = op_expT(obj, funcT, tau);
            expV = op_expV(obj, funcV, tau/2);
            
            if isa(invec, 'FundamentalState')
                dtypecmp(obj, invec.y);                                               
                                
                if strcmp(invec.basis, 'q')
                    qvec = invec.y;
                else
                    qvec = ifft(invec.y);
                end
                                
                pvec =  fft( expV(qvec) );
                qvec = ifft( expT(pvec) );
                pvec =  fft( expV(qvec) );                
                                            
                if strcmp(invec.basis, 'q')
                    outvec = vec2state(obj, ifft(pvec), 'q');
                else
                    outvec = vec2state(obj, pvec, 'p');
                end
            else
                dtypecmp(obj, invec);                
                
                pvec   =  fft( expV(invec) );
                qvec = ifft( expT(pvec) );
                outvec   =  fft( expV(invec) );                                
            end    
        end
        
        function outvec = expTVTevolve(obj, invec, funcT, funcV, tau)
            arguments
                obj
                invec 
                funcT
                funcV
                tau {mustBeNumeric} = 1
            end
                                    
            expT = op_expT(obj, funcT, tau/2);
            expV = op_expV(obj, funcV, tau);
            
            if isa(invec, 'FundamentalState')                
                dtypecmp(obj, invec.y);
                
                if strcmp(invec.basis, 'p')
                    pvec = invec.y;
                else
                    pvec = ifft(invec.y);
                end
                                       
                qvec = ifft( expT(pvec) );
                pvec =  fft( expV(qvec) );
                qvec = ifft( expT(pvec) );              
                                            
                if strcmp(invec.basis, 'q')
                    outvec = vec2state(obj, qvec, 'q');
                else
                    outvec = vec2state(obj, fft(qvec), 'p');
                end
            else
                dtypecmp(obj, invec);
                pvec = invec;
                qvec   = ifft( expT(invec) );
                pvec   =  fft( expV(qvec) );
                qvec = ifft( expT(pvec) );
                outvec = qvec;
            end    
        end
        
        
        
        %         %function op = evolutionOP(obj)
%             op = @(x) EvolveTV(x);
%         end
%         
%         %function setop(obj, funcT, funcV, tau)
%             if exist('tau', 'var')
%                 obj.tau = tau;
%             end           
%             obj.funcT = funcT; % = set(obj, 'funcT', funcT);
%             obj.funcV = funcV; % = set(obj, 'funcV', funcV);
%         end
%         
%         %function op = getop(obj)
%             op = @(x) evolveTV(obj, x)
%         end
% 
%         
%         %function outvec = EvolveTV(obj, invec, funcT, funcV)
%             arguments
%                 obj
%                 invec %(:,1) {mustBeNumeric}
%                 funcT % mustbefunction @funcT
%                 funcV % mustbefunction @funcV
%             end
%             
%             if ~strcmp(obj.dtype, class(invec))
%                 %st = dbstack;
%                 %funcname =st.name;
%                 msg = {"data type error:"
%                     sprintf("%s.dtype=%s, but class(invec)=%s", class(obj), obj.dtype, class(invec) )
%                     };
%                 error(sprintf("%s", msg{:}));
%             end
%             
%             if invec.basis == 'p'
%                 invec = q2p(invec);
%             end            
%                                     
%             qvec = exp( -1j .* funcV(obj.q) / obj.hbar ) .* invec.y; 
%             pvec = fft(qvec);
%             pvec = exp( -1j .* funcT(obj.p) / obj.hbar) .* pvec;
%             qvec = ifft(pvec);
%             outvec = vec2state(obj, qvec, 'q');
%             
% 
%             if invec.basis == 'p'
%                 outvec = outvec.q2p();
%             end            
%             
%         end
%         
%         %function outvec = EvolveVT(obj, invec, funcT, funcV)
%             arguments
%                 obj
%                 invec %(1,:) {mustBeNumeric}
%                 funcT % mustbefunction @funcT
%                 funcV % mustbefunction @funcV
%             end
%            
%             
%             if invec.basis == 'q'
%                 invec = q2p(invec);
%             end            
%             
%             pvec = exp( -1j .* funcT(obj.p) / obj.hbar ) .* invec.y; 
%             qvec = fft(pvec);
%             qvec = exp( -1j .* funcV(obj.q) / obj.hbar) .* qvec;
%             pvec = ifft(qvec);
%             outvec = vec2state(obj, pvec, 'p');
%             
%             if invec.basis == 'q'
%                 outvec = outvec.p2q()
%             end
% 
%         end
%         
%         
%         %function a = AAA(obj)
%             a = 1
%         end
%         
%         
%         %function y = test(obj, func)
%             arguments
%                 obj
%                 func
%             end
%             disp(isa(@func, 'function_handle'))
%             y = 1;
%         end
    end
end

function dtypecmp(obj, vec)
    if ~strcmp(obj.dtype, class(vec))
        msg = {'data type error:'
            sprintf('%s.dtype = %s, but', class(obj), obj.dtype)
            sprintf('class(invec) = %s', class(vec) )
            };
        error(sprintf('%s\n', msg{:}));
    end
end

%function vec = expT4vec(obj, funcT, invec)
%if class(invec) == 'FundamentalState'
%    invec = invec.prep();
%    invec = invec.y;
%end
%vec = @(invec) transpose( exp( -1j * funcT(obj.p) * tau /obj.hbar) )  .* invec;
%end

%function vec = expV4vec(obj, funcV, invec)

%if class(invec) == 'FundamentalState'
%    invec = invec.qrep();
%    invec = invec.y;
%end
%vec = @(invec) transpose( exp( -1j * funcV(obj.q) * tau /obj.hbar) )  .* invec;
%end


