classdef qtmaplab < handle
    methods(Static)
        function Info()
            version = "0.1";
            disp(version);
        end
        
        function Test()
            %clc;
            %clear;
            %close all;
            %set(gcf,'Visible','on')
            runtests("+UnitTest/TestFundamentalState.m")            
        end
    end
end

