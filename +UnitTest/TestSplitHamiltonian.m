classdef TestSplitHamiltonian  < matlab.unittest.TestCase
    %TESTSPLITHAMILTONIAN このクラスの概要をここに記述
    %   詳細説明をここに記述
    
    properties
        fig
        SH
    end
    
    methods
        function obj = TestSplitHamiltonian(inputArg1,inputArg2)
            %TESTSPLITHAMILTONIAN このクラスのインスタンスを作成
            %   詳細説明をここに記述
            obj.Property1 = inputArg1 + inputArg2;
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 このメソッドの概要をここに記述
            %   詳細説明をここに記述
            outputArg = obj.Property1 + inputArg;
        end
    end
end

