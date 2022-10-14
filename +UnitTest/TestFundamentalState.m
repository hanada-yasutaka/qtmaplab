classdef TestState < matlab.unittest.TestCase

    properties
        fig
        state
        sysinfo
        basis
    end
 
    methods(TestMethodSetup)
        function exportpath(obj)
            %addpath("../");
            obj.fig = @figure;
            
            r = randi([0 1], 1);
            if mod(r, 2) == 0
                obj.basis = 'p';
            else
                obj.basis = 'q';
            end

            obj.sysinfo = UnitTest.set_test_system();
            obj.state = State(obj.sysinfo, obj.basis);
            
        end
    end
    
    methods(TestMethodTeardown)
        %function closeFigure(obj)
        %    close(obj.TestFigure)
        %end
    end
     
 
    methods(Test)
        
        function addminus(obj)
            dim = obj.sysinfo.dim;
            v1 = rand(dim, 1) + 1j * rand(dim, 1);
            v2 = rand(dim, 1) + 1j * rand(dim, 1);
            s1 = State(obj.sysinfo, obj.basis, v1);
            s2 = State(obj.sysinfo, obj.basis, v2);
            tol = 1e-15;
            
            s = s1 + s2;
            v = v1 + v2;
            assertEqual(obj, s.y, v, 'AbsTol', tol);

            s = s1 - s2;
            v = v1 - v2;
            assertEqual(obj, s.y, v, 'AbsTol', tol);
            
        end
            
        function inner(obj)
            dim = obj.sysinfo.dim;
            v1 = rand(dim, 1) + 1j * rand(dim, 1);
            s1 = State(obj.sysinfo, obj.basis, v1);
            s1 = s1.normalize();
            a = inner(s1, s1);
            tol = 1e-15;
            assertEqual(obj, a, 1, 'AbsTol', tol);
        end
        
        function scaler_multdiv(obj)
            dim = obj.sysinfo.dim;            
            v1 = rand(dim, 1) + 1j * rand(dim, 1);
            a = rand(1);
            s1 = State(obj.sysinfo, obj.basis, v1);
            s = s1/a;
            tol = 1e-15;
            assertEqual(obj, s.y, v1 / a, 'AbsTol', tol);
            s = s1 .* a;
            assertEqual(obj, s.y, v1 .* a, 'AbsTol', tol);            
        end
            
        
%        function delta(obj)
%            a = obj.state.delta(0.00001, 'verbose', true);
%        end
        
        
        function coherent(obj)
            domain = obj.state.domain;
            dim = obj.state.dim;
            qc = ( domain(1,1) + domain(1,2) )/2;
            pc = ( domain(2,1) + domain(2,2) )/2;
            
            fig = obj.fig();
            ax1 = subplot(2,2,1);
            ax2 = subplot(2,2,4);
            ax0 = subplot(2,2,3);
            hold(ax1, 'on')
            hold(ax2, 'on')
            hold(ax0, 'on')            

            qpc = [(rand(1) - 1/2) + qc, (rand(1) -1/2) + pc];
            cs = obj.state.coherent(qpc(1), qpc(2), 'verbose', true);
            
            plot(ax1, [qpc(1), qpc(1)], [-40, 1], '--k');
            plot(ax2, [-40, 1], [qpc(2), qpc(2)], '--k');                        

            qvec = cs.qrep();

            [~, ind] = max( abs2(qvec) );
            tol = obj.state.hbar*2*pi;
            plot(ax1, qvec.x, log10(abs2(qvec)), 'ok')                        
            assertEqual(obj, qpc(1), qvec.x(ind), 'AbsTol', tol);
            assertTrue(obj, abs(norm(qvec) - 1) < 1e-14);
            
            pvec = cs.prep();
            plot(ax2, log10(abs2(pvec)), pvec.x, 'ok')            
            [~, ind] = max( abs2(pvec) );
            assertEqual(obj, qpc(2), pvec.x(ind), 'AbsTol', tol);
            assertTrue(obj, abs(norm(pvec) - 1) < 1e-14);
            
            qvec2 = pvec.p2q();
            plot(ax1, qvec2.x, log10(abs2(qvec2)), '-r','LineWidth', 2)            
            abs(qvec2.y  - qvec.y);
            assertTrue(obj, all( abs(qvec2.y  - qvec.y) < 1e-14, 'all') )
            
            pvec2 = qvec.q2p();
            plot(ax2, log10(abs2(pvec2)), pvec2.x, '-r','LineWidth', 2)                        
            assertTrue(obj, all( abs(pvec2.y  - pvec.y) < 1e-14, 'all') )
            
            [x,y,z] = cs.hsmrep();
            contour(ax0, x, y, z, 20);
            scatter(ax0, [qpc(1)], [qpc(2)], 50, 'filled',"ok")
                            
            
            xlim(ax1, domain(1,:))
            xlim(ax2, [-40 1])   
            xlim(ax0, domain(1,:))
            
            ylim(ax1, [-40 1])
            ylim(ax2, domain(2,:))
            ylim(ax0, domain(2,:))
                        
            xlabel(ax1, 'q')
            ylabel(ax2, 'p')
            xlabel(ax0, 'q')
            
            ylabel(ax1, '|psi(q)|^2')
            xlabel(ax2, '|psi(p)|^2')               
            ylabel(ax0, 'p')            
            
            %legend(ax1, 'qc','qrep()','qrep().q2p().p2q()','Location','southeast')
            %legend(ax2, 'pc','prep()','prep().p2q().q2p()','Location','southeast')
            
            legend(ax1, 'qc','qrep()','qrep().q2p().p2q()','Position',[0.52 0.8 0.1 0.1])
            legend(ax2, 'pc','prep()','prep().p2q().q2p()','Position',[0.8 0.52 0.1 0.1])
            legend(ax0, sprintf("qc=%.3f pc=%.3f", qpc(1), qpc(2)) )            
            
            hold(ax1, 'off')
            hold(ax2, 'off')
            
            savepath = strcat(pwd, '/UnitTestRes');
            
            if exist(savepath, 'dir')
                fname = 'coherent.png';
                path = sprintf('%s/%s', savepath, fname);
                saveas(gcf, path);
                fprintf("save: %s/\n", path);
            else
                waitforbuttonpress
            end
            %clf(fig)
            close(fig);
            %clf,clc
        end
    end 
end

