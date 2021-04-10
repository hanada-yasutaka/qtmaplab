classdef SimplecticIntegrator < handle
    %Classical Symplectic integrator for the one-dimensional split Hamiltonian H(q, p) = T(p) + V(q)
    %   詳細説明をここに記述
    
    properties
        dT % function of dT/dp
        dV % function of dV/dq
        tau % 
        evolve
    end
    
    methods
        function obj = SimplecticIntegrator(dT, dV, varargin)
            %UNTITLED2 このクラスのインスタンスを作成
            %   詳細説明をここに記述
                        
            obj.dT = dT;
            obj.dV = dV;

            par = inputParser;
            
            addOptional(par, 'tau', 1, @isnumeric)
            addOptional(par, 'order', 1, @(x) x == floor(x) );
            %addOptional(par, 'OPorder', 'TV', @(x) ismember(x, ["TV" "VT" "TVT" "VTV"]) );
            
            parse(par, varargin{:} );            
            obj.tau     = par.Results.tau;            
            order = par.Results.order;
            
            %obj.order = order;
            %obj.tau = 1;
            %obj.evolve = @(q, p) evo
            if order == 1
                obj.evolve = @(x) obj.evolveQP(x);
            elseif order == -1
                obj.evolve = @(x) obj.evolvePQ(x);                               
            elseif order == 2
                obj.evolve = @(x) obj.evolvePQP(x);                
            elseif order == -2
                obj.evolve = @(x) obj.evolveQPQ(x);               

            end
            
        end
        
        %function x = evolve(obj, x)
        %    x = evolveTV(x);
        %    %op = @(x) obj.evolveTV(x);
        %    
        %end
        
        function x = evolveP(obj, x, c)
            arguments
                obj
                x
                c = 1
            end
            q = x(1, :);
            p = x(2, :);
            tau = obj.tau;
            x = [q; p - obj.dV(q) * tau * c];
        end
        function x = evolveQ(obj, x, c)
            arguments
                obj
                x
                c = 1                
            end
            q = x(1, :);
            p = x(2, :);
            tau = obj.tau;
            x = [q + obj.dT(p) * tau * c; p];

        end
                      
        function x = evolveQP(obj, x, z)
            % 1st order simplectic integrator for TV order
            arguments
                obj
                x
                z = 1
            end
            c = [1 1];            
            x = obj.evolveP(x, c(1));
            x = obj.evolveQ(x, c(2));
        end
        
        function x = evolvePQ(obj, x, z)
            % 1st order simplectic integrator for TV order
            arguments
                obj
                x
                z = 1
            end
            c = [1 1];            
            x = obj.evolveQ(x, c(1));
            x = obj.evolveP(x, c(2));
        end
        
        
        function x = evolvePQP(obj, x, z)
            % 1st order simplectic integrator for TV order
            arguments
                obj
                x
                z = 1
            end
            c = [1/2 1];            
            x = obj.evolveP(x, c(1)*z);
            x = obj.evolveQ(x, c(2)*z);
            x = obj.evolveP(x, c(1)*z);            
        end
        
        function x = evolveQPQ(obj, x, z)
            % 1st order simplectic integrator for TV order
            arguments
                obj
                x
                z = 1
            end
            c = [1/2 1];            
            x = obj.evolveQ(x, c(1)*z);
            x = obj.evolveP(x, c(2)*z);
            x = obj.evolveQ(x, c(1)*z);            
        end
        
        
        
        
        
        
        function traj = getTraj(obj, x, tmax)
            traj = [[], []];
            for i=1:tmax
                x = obj.evolveTV(x);
                traj = horzcat(traj, x);
            end                      
        end
        
        %function x = evoVT(obj, x, tau)
            % 1st order simplectic integrator for VT order            
        %    q = q + obj.dT(p)*tau;            
        %    p = p - obj.dV(q)*tau;
        %end
                
        %function [q, p] = evo_TVT(obj,q, p, tau)
            % 1st order simplectic integrator for TVT order            
        %    q = q + obj.dT(p)*tau/2;
        %    p = p - obj.dV(q)*tau;
        %    q = q - obj.dT(p)*tau/2;
        %end
        
        %function outputArg = evo_TVT(obj,q, p, tau)
            %METHOD1 このメソッドの概要をここに記述
            %   詳細説明をここに記述
        %    p = p + obj.dV(q)*tau/2;
        %    q = q - obj.dT(p)*tau;
        %    p = p - obj.dV(q)*tau/2;
        %end
        
        
        
    end
end

