classdef F7 < PROBLEM
% <multi/many> <real> 

    methods
        %% Default settings of the problem
        function Setting(obj)
            if isempty(obj.M); obj.M = 2; end
            if isempty(obj.D); obj.D = obj.M+8; end
             obj.lower    = zeros(1,obj.D);
            obj.upper    = ones(1,obj.D);
            obj.encoding = 'real';
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            [cou,~]=size(PopDec);
            xi=PopDec(:,2:end);x1=PopDec(:,1);
            yi=xi-sin(0.5*pi*xi);
            g=2*sin(0.5*pi*x1).*((obj.D-1)+sum(yi.^2-cos(2*pi*yi),2));
            f1=(1+g).*x1;
            f2=0.5*(1+g).*( (1-x1.^0.1)+  (1-x1.^0.5).^2.*   (cos(3*pi*x1).^2)    );
            PopObj = [f1,f2];
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            R = UniformPoint(N,obj.M);
%             temp=R;
%             [cou,~]=size(R);
%              cout=0;
%             for index=1:cou
%                 if((temp(index,1)>0.18&&temp(index,1)<0.4)||(temp(index,1)>0.52&&temp(index,1)<0.68))
%                     R(index-cout,:)=[];
%                      cout=cout+1;
%                 end
%             end
            [cou,~]=size(R);
            for index=1:cou
                x=R(index,1);
                R(index,2)=0.5*(  (1-x.^0.1)+  (1-x.^0.5).^2.*   (cos(3*pi*x).^2)  );
            end
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            if obj.M == 2
%                 R = obj.GetOptimum(100);
            x      = linspace(0,1,10000)';
                y      = 0.5*(  (1-x.^0.1)+  (1-x.^0.5).^2.*   (cos(3*pi*x).^2)  );
                nd     = NDSort([x,y],1)==1;
                x(~nd) = nan;
                R      = [x,y];
            elseif obj.M == 3
                a = linspace(0,pi/2,10)';
                R = {sin(a)*cos(a'),sin(a)*sin(a'),cos(a)*ones(size(a'))};
            else
                R = [];
            end
        end
    end
end