function [x,fbak] = ISPO(fittFun, param)
% function [x,fbak, plt] = ISPO(fittFun, param)
%Algorithm implemented according to the presented in the papper "Simplified
% Intelligence Single Particle Optimization Based Neural Network for Digit 
% Recognition".
%   The algorithm receives:
%   param.A:        Acceleration
%   param.P:        Acceleration power factor
%   param.B:        Learning coefficient
%   param.S:        Learning factor reduction ratio
%   param.eps:      Minimum threshold on learning factor
%   param.H:        Particle learning step
%   param.D:        Problem dimension
%   param.maxEval:  Max number of function evaluations
%   param.lowLim:   Lower bound of the decision variables(vector)
%   param.upLim     Upper bound of the decision variables(vector)
%   fittFun:        Objective function.
%   The algorithm returns:
%   x:              Vector containing the best solution
%   fback:          Fitness value of the best solution
%   plt:            If uncomented vector containing the history of fitness

    %Read all parameters to local variables, (just for readability).
    A   = param.A;
    P   = param.P;
    B   = param.B;
    S   = param.S;
    eps = param.eps;
    H   = param.H;
    D   = param.D;
    maxEval = param.maxEval;
    lowLim    = param.lowLim;
    upLim     = param.upLim;
    
    %Start total evals
    tot_Eval = 0;

%     %Start variables used to save plot each 10D function evaluations
%     smpl    = 10*D;
%     plt     = zeros(1001,1);
%     indPlt  = 1;
    
    %start the best solution
    x   = unifrnd(lowLim, upLim, [D,1]);
    fx  = fittFun(x);
    tot_Eval = tot_Eval + 1;
    
%     %Save the first point
%     plt(indPlt) = fx;
%     indPlt   = indPlt + 1;
    
    %Start the old fitt value in memory to future comparations
    fbak    = fx;

    while tot_Eval < maxEval
        for i = 1:D
            L   = 0;
            xi_bak  = x(i);
            %Update each decsion variable----------------------------------
            for t = 1:H
                v   = (A/(t^P))*(rand-0.5)+B*L;
                x(i)= x(i)+v;
                
                %Fix bound constrains--------------------------------------
                while x(i) > upLim
                    x(i)    = x(i)-(upLim - lowLim);
                end
                while x(i) < lowLim
                    x(i)    = x(i)+(upLim - lowLim);
                end
                %Fix bound constrains--------------------------------------
                
                %Eval new solution-----------------------------------------
                fx  = fittFun(x);
                tot_Eval = tot_Eval + 1;
                %Eval new solution-----------------------------------------
                
                %Compare fittness of the new solution----------------------
                if fx <= fbak   %update learning factor if improvement
                    L       = v;        
                    fbak    = fx;   %update new position
                    xi_bak  = x(i);
                else  %update learning factor if no improvement
                    if L ~= 0
                        L   = L/S;
                    end
                    
                    if L < eps
                        L   = 0;
                    end
                    x(i)  = xi_bak; %Restore old position if no improvement
                end
                %Compare fittness of the new solution----------------------

%                 %Save best solution in array-----------------------------
%                 if mod(tot_Eval,smpl)==0
%                     plt(indPlt) = fbak;
%                     indPlt  = indPlt + 1;
%                 end
%                 %Save best solution in array-----------------------------
                
                %Check if the Max num of FE is reached---------------------
                if(tot_Eval >= maxEval)
                    break;
                end
                %Check if the Max num of FE is reached---------------------
            end
            %Update each decsion variable----------------------------------
            
            %Check if the Max num of FE is reached-------------------------
            if(tot_Eval >= maxEval)
                break;
            end
            %Check if the Max num of FE is reached-------------------------
        end
    end
    disp(tot_Eval)
end