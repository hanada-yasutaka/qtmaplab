classdef qtmaplab < handle
    % qtmaplab provides a solver of an eigenvalue problem and a time-evolution of wavepacket for the quantum map (symplectic integrator). 
    % It supports a double-precision (by Matlab) and a multiple-precision by using Advanpix toolbox.
    methods(Static)
        function Info()
            % Info displays the information of qtmaplab
            version = "0.3";
            disp(version);
        end
        
        function Test()
            % Test runs the UnitTests
            runtests("+UnitTest/TestFundamentalState.m")            
        end
    end
end
