function [xgb,fittXgb] = cPSO(fittFun, param)
%compact Particle Swarm Optimization algorithm, the function receives the
%fitness function and a variable param, with the follow values
%   param.lowLim:   Vector with the low limits of the decision variables
%   param.upLim:    Vector with the up limits of the decision variables
%   param.NP:       Value of the virtual population
%   param.D:        Number of decision variables
%   param.maxEval:  Max number of function evaluations
%   lowLim  = param.lowLim;

%All values are copied to loval variables
    upLim   = param.upLim;
    lowLim  = param.lowLim;
    Np      = param.NP;
    D       = param.D;
    maxEval = param.maxEval;

%Save in memory the vectors that will be used in the algorithm
    muV     = zeros(D,1);   %\,PV
    stdV    = zeros(D,1);   %/
    xgb     = zeros(D,1);   %Global best
    xlb     = zeros(D,1);   %Local best
    xt      = zeros(D,1);   %trial
    vt      = 2*rand(D,1) - 1;  %velocity vector
    
    c0=-0.2;
    c1=-0.07;
    c2=3.74;

% Init the vectors samplig solution by means of PV
    for i = 1:D
        muV(i)  = 0;
        stdV(i) = 10;
        xgb(i)  = sampleSolution(muV(i),stdV(i));
        xt(i)   = sampleSolution(muV(i),stdV(i));
    end
    fittXgb = fittFun(denorm(xgb, lowLim, upLim)); %Get global best value
    tot_evals = 1;

    while tot_evals < maxEval
        for i=1:D
            r1=rand;
            r2=rand;
            xlb(i)  = sampleSolution(muV(i),stdV(i));   %sample randmly local best
            vt(i)   = c0*vt(i) + c2*r1*(xgb(i)-xt(i)) + c1*r2*(xlb(i)-xt(i));   %Update velocity of particle
            vt(i)   = min([vt(i), 1]);                  %Set bounds to the velocity of particle
            vt(i)   = max([vt(i),-1]);
            xt(i)   = xt(i) + vt(i);                    %Update position of xt

            if(xt(i) < -1)          %The particle position is update in a toroidal way
                xt(i) = xt(i)+2;
            elseif (xt(i) > 1)
                xt(i) = xt(i)-2;
            end

        end

        %Eval fitness value of x trial and local best
        fittXt = fittFun(denorm(xt, lowLim, upLim));
        fittXlb= fittFun(denorm(xlb, lowLim, upLim));
        tot_evals = tot_evals + 2;
        
        %Compete and update the PV accordigly
        %the mean is also updated in a toroidal way
        if(fittXt < fittXlb)
            for i=1:D
                mu = muV(i);
                muV(i)  = muV(i) + (1/Np)*(xt(i) - xlb(i));
                if(muV(i) < -1)
                    muV(i) = muV(i)+2;
                elseif (muV(i) > 1)
                    muV(i) = muV(i)-2;
                end
                stdV(i) = sqrt(abs(stdV(i)^2 + mu^2 - muV(i)^2 + (1/Np)*(xt(i)^2 - xlb(i)^2))); %update of std
                stdV(i) = min([10, stdV(i)]);  %Tis is to avoid the grow of std, the max value is 10
            end
        else
            %same as If part wit xt and xlb positions changed
            for i=1:D
                mu = muV(i);
                muV(i)  = muV(i) + (1/Np)*(xlb(i)-xt(i));
                if(muV(i) < -1)
                    muV(i) = muV(i)+2;
                elseif (muV(i) > 1)
                    muV(i) = muV(i)-2;
                end
                stdV(i) = sqrt(abs(stdV(i)^2 + mu^2 - muV(i)^2 + (1/Np)*(xlb(i)^2 - xt(i)^2)));
                stdV(i) = min([10, stdV(i)]);
            end
        end
        %If xt is better solution than global best is updated (compete).
        if(fittXt < fittXgb)
            fittXgb = fittXt;
            xgb     = xt;
        end
    end
end