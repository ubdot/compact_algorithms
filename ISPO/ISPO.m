function [x,fbak] = ISPO(fittFun, param)
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
    tot_Eval = 0;

    x   = unifrnd(lowLim, upLim, [D,1]);
    fx  = fittFun(x);

    
    tot_Eval = tot_Eval + 1;
    while tot_Eval < maxEval
        for i = 1:D
            L   = 0;
            fbak    = fx;
            xi_bak  = x(i);
            %Update each decsion variable----------------------------------
            for t = 1:H
                v   = (A/(t^P))*(rand-0.5)+B*L;
                x(i)= x(i)+v;
                
                while x(i) > upLim
                    x(i)    = x(i)-(upLim - lowLim);
                end
                while x(i) < lowLim
                    x(i)    = x(i)+(upLim - lowLim);
                end

                fx  = fittFun(x);
                tot_Eval = tot_Eval + 1;
                %Compare fittness of the new solution----------------------
                if fx <= fbak   %update learning factor if improvement
                    L   = v;        
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
            end
            %Update each decsion variable----------------------------------
        end
    end
end