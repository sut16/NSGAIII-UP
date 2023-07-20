classdef F1 < PROBLEM
% <multi/many> <real> 

    methods
        %% Default settings of the problem
        function Setting(obj)
            if isempty(obj.M); obj.M = 2; end
            if isempty(obj.D); obj.D = obj.M+8; end
            obj.lower    = [zeros(1,1),-ones(1,obj.D-1)];
            obj.upper    = ones(1,obj.D);
            obj.encoding = 'real';
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            [cou,~]=size(PopDec);
            cal_sub=2:obj.D;
            g= PopDec(:,2:end)-repmat(0.9*sin(cal_sub.*(pi/obj.D)),cou,1);
            PopDec_odd= 3:2:obj.D;PopDec_even=2:2:obj.D;
            g_even=g(:,2:2:end);g_odd=g(:,1:2:end);
            norm_odd=sqrt(sum(repmat(PopDec_odd,cou,1).^2,2));
            SumG_odd=sum( abs(g_odd).^0.7 ,2);
            f1=1-cos(0.5*pi*PopDec(:,1))+2./norm_odd.*SumG_odd;
            f2=10-10*sin(0.5*pi*PopDec(:,1))+2./sqrt(sum(repmat(PopDec_even,cou,1).^2,2)).*sum( abs(g_even).^0.7 ,2);
            PopObj = [f1,f2];
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            R = UniformPoint(N,obj.M);
            R(:,2)=10*(1-sqrt((1-((R(:,1)-1).^2))));
           
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            if obj.M == 2
                R = obj.GetOptimum(10000);
            elseif obj.M == 3
                a = linspace(0,pi/2,10)';
                R = {sin(a)*cos(a'),sin(a)*sin(a'),cos(a)*ones(size(a'))};
            else
                R = [];
            end
        end
    end
end