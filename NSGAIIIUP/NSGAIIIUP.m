classdef NSGAIIIUP < ALGORITHM
% <multi/many> <real/binary/permutation> <constrained/none>
  methods
        function main(Algorithm,Problem)
            num=4; %% Utopian factor
            [Z,Problem.N] = UniformPoint(Problem.N,Problem.M);
            Population    = Problem.Initialization();
            Zmin          = min(Population(all(Population.cons<=0,2)).objs,[],1);
            Z_org=Z;
            Z=Z+num*ones(1,Problem.M);
            %B  = eye(Problem.D);
            %m  = 0.5*(Problem.upper - Problem.lower);
            %ps = 0.5;
            while Algorithm.NotTerminated(Population)
               MatingPool = TournamentSelection(2,Problem.N,sum(max(0,Population.cons),2));
                Offspring  = OperatorGA(Population(MatingPool));
                
                %RSBX operator
                % Offspring  = ARSBX(Population(MatingPool),{B,m,ps});
               
                Zmin       = min([Zmin;Offspring(all(Offspring.cons<=0,2)).objs],[],1);
                Population = EnvironmentalSelection([Population,Offspring],Problem.N,Z,Zmin,num);
            end
        end
    end
end
