function [xgb,fittXgb] = cPSO(fittFun, param)
% function [xgb,fittXgb, plt] = cPSO(fittFun, param)
%compact Particle Swarm Optimization algorithm, the function receives the
%fitness function and a variable param, with the follow values
%   param.lowLim:   Vector with the low limits of the decision variables
%   param.upLim:    Vector with the up limits of the decision variables
%   param.NP:       Value of the virtual population
%   param.D:        Number of decision variables
%   param.maxEval:  Max number of function evaluations
%   param.c0:       Constant to scalate the velocity
%   param.c1:       Constant to scalate lb difference
%   param.c2:       Constant to scalate gb difference
%   lowLim  = param.lowLim;

%All values are copied to local variables
    upLim   = param.upLim;
    lowLim  = param.lowLim;
    Np      = param.NP;
    D       = param.D;
    maxEval = param.maxEval;
    c0      = param.c0;
    c1      = param.c1;
    c2      = param.c2;

%Save in memory the vectors that will be used in the algorithm
    muV     = zeros(D,1);   %\,PV
    stdV    = zeros(D,1);   %/
    xgb     = zeros(D,1);   %Global best
    xlb     = zeros(D,1);   %Local best

%Init xt and vt
    xt      = unifrnd(-1,1,[D,1]);  %[-1, 1] according to paper
    vt      = unifrnd(-1,1,[D,1]);  %[0, 1] according to paper
   
%     %Start variables used to save plot each 10D function evaluations
%     smpl    = 10*D;
%     plt     = zeros(1001,1);
%     indPlt  = 1;

% Init the xgb solution by means of PV-------------------------------------
    for i = 1:D
        muV(i)  = 0;
        stdV(i) = 10;
        xgb(i)  = sampleSolution(muV(i),stdV(i));
%         xt(i)   = sampleSolution(muV(i),stdV(i));
    end    
    fittXgb = fittFun(denorm(xgb, lowLim, upLim)); %Get global best value
    tot_evals = 1;
% Init the xgb solution by means of PV-------------------------------------


%     %Save the first point
%     plt(indPlt)    = fittXgb;
%     indPlt  = indPlt + 1;


    while tot_evals < maxEval
        %Update particle---------------------------------------------------
        for i=1:D
            %sample xlb by means of PV-------------------------------------
            xlb(i)  = sampleSolution(muV(i),stdV(i));
            %sample xlb by means of PV-------------------------------------

            %Velocity and position update----------------------------------
            r1=rand;
            r2=rand;
            vt(i)   = c0*vt(i) + c1*r2*(xlb(i)-xt(i)) + c2*r1*(xgb(i)-xt(i));   
            %Fix velocity bounds-------------------------------------------
            vt(i)   = min([vt(i), 1]);                  
            vt(i)   = max([vt(i),-1]);
            %Fix velocity bounds-------------------------------------------
            xt(i)   = xt(i) + vt(i);%Update position of xt
            %Fix bounds----------------------------------------------------
            if(xt(i) < -1)          
                xt(i) = xt(i)+2;
            elseif (xt(i) > 1)
                xt(i) = xt(i)-2;
            end
            %Fix bounds----------------------------------------------------
            %Velocity and position update----------------------------------
        end
        %Update particle---------------------------------------------------

        %competition-------------------------------------------------------
        fittXt = fittFun(denorm(xt, lowLim, upLim));
        tot_evals = tot_evals + 1;
        %Prevent exceed  the budget----------------------------------------
        if tot_evals < maxEval
            fittXlb= fittFun(denorm(xlb, lowLim, upLim));
            tot_evals = tot_evals + 1;
            
            %Compete xt and xlb--------------------------------------------
            if(fittXt < fittXlb)
                [muV,stdV]= upd_PV_D(muV, stdV, xt, xlb, Np);
                %update xlb if is outperformed-----------------------------
                if(fittXt < fittXgb) 
                    fittXgb = fittXt;
                    xgb     = xt;
                end
                %update xlb if is outperformed-----------------------------
            else
                [muV,stdV]= upd_PV_D(muV, stdV, xlb, xt, Np);
                %update xlb if is outperformed-----------------------------
                if(fittXlb < fittXgb)
                    fittXgb = fittXlb;
                    xgb     = xlb;
                end
                %update xlb if is outperformed-----------------------------
            end
            %Compete xt and xlb--------------------------------------------
        elseif fittXt < fittXgb
            fittXgb = fittXt;
            xgb     = xt;
        end
        %Prevent exceed  the budget----------------------------------------
        %competition-------------------------------------------------------
        %Save best solution in array---------------------------------------
%         if mod(tot_evals-1,smpl)==0
%             plt(indPlt) = fittXgb;
%             indPlt      = indPlt + 1;
%         end
        %Save best solution in array---------------------------------------
    end
%     disp([muV, stdV]);
    disp(tot_evals);
end