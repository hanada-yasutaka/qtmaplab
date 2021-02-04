classdef StateTest < matlab.unittest.TestCase

    properties
        %TestFigure
        state
    end
 
    methods(TestMethodSetup)
        function createFigure(obj)
            % comment
            %obj.TestFigure = figure;
            dim = randi([10 100], 1);
            r = randi([0 1], 1);
            if r
                domain = [-5 5; -2 5];
            else
                domain = mp('[0 5; -2 5]');
            end
            basis = 'q';

            scl=SystemInfo(dim, domain, basis);
            obj.state = FundamentalState(scl);
        end
    end
 
    methods(TestMethodTeardown)
        function closeFigure(obj)
            %close(testCase.TestFigure)
        end
    end
 
    methods(Test)
 
        function coherent(obj)
            domain = obj.state.domain;
            disp(class(domain))

            qpc = domain(1,1) + ( domain(1,2) - domain(1,1) ) * rand(1, 2, class(domain));
            cs = obj.state.coherent(qpc(1), qpc(2), false);
            normtest(obj, cs)
            %cs = obj.state.coherent(qpc(1), qpc(2), true);            

        end
 
        %function defaultCurrentObject(obj)
        %    csp = obj.state.coherent(-4.9, -1.1, false);
        %end
 
    end 
end


function normtest(testobj, vec)
  if isa(vec.y, 'mp')
      tolerance = vec.eps * 2000
  else
      tolerance = vec.eps * 200;
  end
  y = vec.norm() - 1;
  disp(y);
  testobj.verifyLessThanOrEqual(y, tolerance);
end