function [elite,fittEle] = rcGA(fittFun, param)
%Real value genetic algorithm, the function receives the fitness function
%and a param with the following values:
%   param.F:        The value of F for the differential evolution operator
%   param.lowLim:   Vector with the low limits of the decision variables
%   param.upLim:    Vector with the up limits of the decision variables
%   param.CR:       Probability for the binary crossover
%   param.NP:       Value of the virtual population
%   param.D:        Number of decision variables
%   param.maxEval:  Max number of function evaluations
%   param.pe:       Bool value that decide if the algorithm  is persistent
%   param.age:      If the algorithm isn't persistent defines the max value
%   Returns the best solution found elite, and the fitess value.
    lowLim  = param.lowLim;
    upLim   = param.upLim;
    CR      = param.CR; %by default the algorithm does not have crossover
    Np      = param.NP;
    D       = param.D;
    maxEval = param.maxEval;
    persis  = param.pe;
    ageDE   = param.age;

    %Save in memory the used vectors.
    muV     = zeros(D,1);
    stdV    = zeros(D,1);
    elite   = zeros(D,1);   %elite
    xoff    = zeros(D,1);   %offspring
    ageAct  = 0;            %used if non persistent is enabled

    %init the elite vector
    for i=1:D
        muV(i)  = 0;
        stdV(i) = 10;
        elite(i)  = sampleSolution(muV(i),stdV(i));
    end
    fittEle = fittFun(denorm(elite, lowLim, upLim));
    

    for nEval = 2 : maxEval
        for i=1:D
            xoff(i) = sampleSolution(muV(i),stdV(i));
        end
%         %By default the algorithm does not have crossover so descoment if
%         %used
        for i=1:D
            %binary crossover
            if(rand > CR)
                xoff(i) = elite(i);
            else
                xoff(i) = sampleSolution(muV(i),stdV(i));
            end
        end
        %Eval the fitness of the offspring
        fittOff = fittFun(denorm(xoff, lowLim, upLim));
        %compete between elite and xoff and the PV is updated accordigly 
        %the mean is updated in a toroidal way
        if(fittOff < fittEle || (ageDE <= ageAct && ~persis))
            [muV, stdV] = upd_PV(muV, stdV, xoff, elite, Np);
            elite   = xoff;
            fittEle = fittOff;
            ageAct  = 0; %Used if non persistent
        else
            [muV, stdV] = upd_PV(muV, stdV, elite, xoff, Np);
            ageAct  = ageAct + 1; %used if non persistent
        end
    end
end