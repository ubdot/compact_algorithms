function [elite,fittele] = CSM(fitt,param)
% function [elite,fittele, plt] = CSM(fitt,param)
    %Based on the papper "", hibrid between one solution and compact
    %algorithm, some considerations has been done that the original article
    %does not mention.
    %The algorithm receives parameters:
    %   param.D:    Dimension of the problem
    %   param.NP:   Poblation size(virtual, to update the model)
    %   param.B:    Value used to apply the non uniform mutation 
    %   param.maxEval:  max value of objective function evaluations
    %   param.lowLim:   Vector containing the low bounds of the decision variables
    %   param.upLim:    Vector containing the up bounds of the decision variables  
    %   fitt:           Objective function

    %Parameter initialization----------------------------------------------
    D   = param.D;
    Np  = param.NP;
    B   = param.B; 
    maxFE   = param.maxEval;
    lowlim  = param.lowLim;
    uplim   = param.upLim;
    
    muV     = zeros(D,1);
    stdV    = 10*ones(D,1);
    elite   = zeros(D,1);
    lbest   = zeros(D,1);
    xoff    = zeros(D,1);
   
    budget1 = floor(0.2*maxFE); %Budget assignment values set according to 
    budget2 = ceil(0.8*maxFE);  %values mentioned in the original papper
    totFE = 0;    

    % Reserve space to save best solutions in array
%     smpl    = 10*D;
%     plt     = zeros(1001,1);
%     indPlt  = 1;
    %Parameter initialization----------------------------------------------
    
    %Solution initialization-----------------------------------------------
    for i=1:D
        elite(i) = sampleSolution(muV(i), stdV(i)); %Used for both parts
        lbest(i) = sampleSolution(muV(i), stdV(i)); %Used for the second part of the algorithm
    end
    fittlb  = fitt(denorm(elite, lowlim, uplim));
    fittele = fitt(denorm(lbest, lowlim, uplim));
    totFE   = totFE + 2;
    %save first solution
%     plt(indPlt)    = fittele;
%     indPlt  = indPlt + 1;
    %Solution initialization-----------------------------------------------
    while totFE < maxFE
        if totFE < budget1
        %Single non uniform mutation on elite -----------------------------
        %The article does not specify if is distintion between the
        %different parts to apply the mutation, for the tests it was
        %selected that the mutation is applied considering the budget(not
        %the max FE) to apply the mutation (with distintion)
            xoff    = elite;
            xoff    = SNUM(xoff, budget1, totFE, B);%--
%             xoff    = SNUM(xoff, maxFE, totFE, B);
            fittxoff    = fitt(denorm(xoff, lowlim, uplim));
            totFE = totFE + 1;
            %competition---------------------------------------------------
            if fittxoff <= fittele
                fittele = fittxoff;
                elite   = xoff;
            end
            %competition---------------------------------------------------

            %save best solution
%             if mod(totFE,smpl)==0
%                 plt(indPlt) = fittele;
%                 indPlt  = indPlt + 1;
%             end
        %Single non uniform mutation on elite -----------------------------
        else
        %Multiple non uniform mutation on compact scheme-------------------
            if totFE == budget1
%                 muV = elite;%--
                muV = -elite; %According to papper note
            end
            for i=1:D
                xoff(i) = sampleSolution(muV(i), stdV(i));
            end
            xoff = MNUM(xoff, budget2, totFE-budget1 + 1, B);%--
%             xoff = MNUM(xoff, maxFE, totFE, B);
            fittxoff    = fitt(denorm(xoff, lowlim, uplim));
            totFE       = totFE + 1;
            
            %competition---------------------------------------------------
            if(fittxoff < fittlb)
                [muV, stdV] = upd_PV_D(muV, stdV, xoff, lbest, Np);%Non toroidal way, and fixed max and min value
                lbest   = xoff;
                fittlb  = fittxoff;
            else
                [muV, stdV] = upd_PV_D(muV, stdV, lbest, xoff, Np);
            end
            if fittlb < fittele
                fittele = fittlb;
                elite   = lbest;
            end
            %competition---------------------------------------------------
            %save best solution
%             if mod(totFE,smpl)==0
%                 plt(indPlt) = fittele;
%                 indPlt  = indPlt + 1;
%             end
        %Multiple non uniform mutation on compact scheme-------------------
        end
    end
    
%     disp(fittele)
%     disp([muV, stdV])
end

function xoff = MNUM(xoff_p, MaxIt, it, B)
    %NUM, is applied to all variables on the solution vector, to attact the
    %problem as not separable
    xoff    = xoff_p;
    rows    = size(xoff_p, 1);
    for i = 1:rows
        if rand < 0.5
            xoff(i) = xoff(i) + del_num(MaxIt, it, 1 - xoff(i), B);
        else
            xoff(i) = xoff(i) - del_num(MaxIt, it, xoff(i) - (-1), B);
        end
        
        %fix bounds--------------------------------------------------------
        while xoff(i) < -1
            xoff(i) = xoff(i)+2;
        end
        while xoff(i) > 1
            xoff(i) = xoff(i)-2;
        end
        %fix bounds--------------------------------------------------------
    end
end

function xoff=SNUM(xoff_p, MaxIt, it, B)
    %NUM on a single variable on the solution vector, to attact the problem
    %as serparable
    xoff    = xoff_p;
    rows    = size(xoff_p,1);
    i   = randi(rows);
    if rand < 0.5
        xoff(i) = xoff(i) + del_num(MaxIt, it, 1 - xoff(i), B);
    else
        xoff(i) = xoff(i) - del_num(MaxIt, it, xoff(i) - (-1), B);
    end

    %fix bounds------------------------------------------------------------
    while xoff(i) < -1
        xoff(i) = xoff(i)+2;
    end 
    while xoff(i) > 1
        xoff(i) = xoff(i)-2;
    end
    %fix bounds------------------------------------------------------------
end

function delt=del_num(MaxIt, it, y, B)
    %Detlta used to apply the mutation
    delt    = y*(1-rand^((1-(it/MaxIt))^B));
end