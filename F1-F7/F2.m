classdef F2 < PROBLEM
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
            %SumGandSin_odd=sum(( sqrt(sum((g_odd+sin(pi*g_odd)/pi).^2,2)) ),2);
            SumGandSin_odd=sum(abs(g_odd+sin(pi*g_odd)/pi),2);
            norm_even=sqrt(sum(repmat(PopDec_even,cou,1).^2,2));
            %SumGandSin_evn=sum( sqrt(sum((g_even+sin(pi*g_even)/pi).^2,2)),2);
            SumGandSin_even=sum(abs(g_even+sin(pi*g_even)/pi),2);
            f1=PopDec(:,1)+2./norm_odd.*SumGandSin_odd;
            f2=zeros(cou,1);
            for index=1:cou
                if(PopDec(index,1)<=0.05)
                    f2(index)=1-19*PopDec(index,1)+2/norm_even(index)*SumGandSin_even(index);
                else
                    f2(index)=1/19-(1/19)*PopDec(index,1)+2/norm_even(index)*SumGandSin_even(index);
                end
            end
            PopObj = [f1,f2];
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            R = UniformPoint(N,obj.M);
            [cou,~]=size(R);
            k=-0.95/0.05;
            for index=1:cou
                if(R(index,1)<=0.05)
                    R(index,2)=k*R(index,1)+1;
                else
                    R(index,2)=(1/k)*(R(index,1)-1);
                end
            end
           
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