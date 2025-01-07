function [elite,fittEle] = CScDE(fittFun, param)
%compact Diferential evolution with DE/rand/1/bin operator and Compound  
% Sinusoidal parameter adaptation. Implemented according to the papper "A 
% compact compound sinusoidal differential evolution algorithm for solving 
% optimisation problems in memory-constrained environments"
%   The function receives the fitness function and a parameter that must
%   have the next values:
%   param.nWaves:   Number of waves being composed
%   param.freq:     Frequency of the sinusoidal formula
%   param.lowLim:   Vector with the low limits of the decision variables
%   param.upLim:    Vector with the up limits of the decision variables
%   param.NP:       Value of the virtual population
%   param.D:        Number of decision variables
%   param.maxEval:  Max number of function evaluations



    %Load parameters to local variables
    lowLim  = param.lowLim;
    upLim   = param.upLim;
    Np      = param.NP;
    D       = param.D;
    maxEval = param.maxEval;
    nWaves  = param.nWaves;
    freq    = param.freq;

    muV     = zeros(D,1);   %means
    stdV    = zeros(D,1);   %Standard deviation
    elite   = zeros(D,1);   %Elite value


    for i=1:D
        muV(i)  = 0;    %\,Provability vector
        stdV(i) = 10;   %/
        elite(i)  = sampleSolution(muV(i),stdV(i)); 
    end

    fittEle = fittFun(denorm(elite, lowLim, upLim)); %Get fitness value elite, is necesary denormalize the solution.

    for nEval = 2 : maxEval
        % F and CR values update-------------------------------------------
        F_it    = 0;
        for w=1:nWaves
            F_it    = F_it + (0.5*sin(2*pi*freq*(nEval-2)/w) + 1);
        end
        F_it    = (1/nWaves)*F_it;
        CR_it   = 0.6 + 0.1*F_it;
        % F and CR values update-------------------------------------------
        % Exponential crossover--------------------------------------------
        xoff=elite;
        i = randi(D);
        pos_temp=i;
            
        xr  = sampleSolution(muV(i),stdV(i));
        xs  = sampleSolution(muV(i),stdV(i));
        xt  = sampleSolution(muV(i),stdV(i));

        xoff(i) = xt + F_it*(xr - xs);
        if(xoff(i) < -1)
            xoff(i) = xoff(i)+2;
        elseif (xoff(i) > 1)
            xoff(i) = xoff(i)-2;
        end
        
        while rand <= CR_it
            i=i+1;
            if i > D
                i=1;
            end

            if i== pos_temp %To prevent inf cicles
                break;
            end
            
            xr  = sampleSolution(muV(i),stdV(i));
            xs  = sampleSolution(muV(i),stdV(i));
            xt  = sampleSolution(muV(i),stdV(i));

            xoff(i) = xt + F_it*(xr - xs);

            if(xoff(i) < -1)
                xoff(i) = xoff(i)+2;
            elseif (xoff(i) > 1)
                xoff(i) = xoff(i)-2;
            end

        end
        % Exponential crossover--------------------------------------------
        % Competition------------------------------------------------------
        fittOff = fittFun(denorm(xoff, lowLim, upLim));
        %If the offspring is better than the original elite the PV is updated accordingly.
        if(fittOff < fittEle)
            [muV, stdV] = upd_PV_D(muV, stdV, xoff, elite, Np); %No toroidal way update
            elite   = xoff;
            fittEle = fittOff;
        else % this part is the same of the if part just invert the elite and offspring to update the PV
            [muV, stdV] = upd_PV(muV, stdV, elite, xoff, Np);   %No toroidal way update
        end
        % Competition------------------------------------------------------
    end
end