function [elite,fittEle] = OEMDE(param,fitt)
% function [elite,fittEle, plt] = OEMDE(param,fitt)
% The algorithm receives :
% param.pop:    population size
% param.CR:     crossover probability
% param.D:      problem dimension
% param.upLim:  problem constrains up lim
% param.lowLim: problem constrains low lim
% param.maxEval:max number of objective functions evaluations


%   Detailed explanation goes here
    Np      = param.pop;
    D       = param.D;
    upLim   = param.upLim;
    lowLim  = param.lowLim;
    maxEval = param.maxEval;
    CR      = param.CR;
    ind1    = zeros(6,1); %used to sample pop to use DE scheme
%6 is the minimun quantity required to execute the different DE schemes

%used to save the best solutions
%     smpl    = 10*D;
%     plt     = zeros(1001,1);
%     indPlt  = 1;
    
    break_alg   = false; %used to stop algorithm if maxFE is reached
    
    %reserve memory to save solutions--------------------------------------
    x   = zeros(2*Np,D+1); %Array form [fx, x_0, x_1,..., x_D]
    totEval = 0;
    for i=1:Np
        x(i,2:D+1)=unifrnd(lowLim,upLim,[1,D]); %init randomly solution
        x(i,1)=fitt(x(i,2:D+1)');
        totEval = totEval + 1;
    end
    %reserve memory to save solutions--------------------------------------
    
    %save inital solution
%     plt(indPlt)    = min(x(1:Np,1));
%     indPlt  = indPlt + 1;

    while(totEval < maxEval)
        %Calculation opositon based solution-------------------------------
        for i=1:Np
            x(Np+i,2:D+1)=(lowLim + upLim)*ones(1,D) - x(i,2:D+1);
            x(Np+i,1)=fitt(x(Np+i,2:D+1)');
            totEval = totEval + 1;
            
            %save best solution
%             if mod(totEval,smpl)==0
%                 plt(indPlt) = min(x(:,1));
%                 indPlt  = indPlt + 1;
%             end

            %break the cycle if max solutions is reached-------------------
            if totEval >= maxEval
                break_alg =true;
                break;
            end
            %break the cycle if max solutions is reached-------------------
        end
        %Calculation opositon based solution-------------------------------

        %Set to work just the Np best solutions----------------------------
        [~,ind]     = mink(x(:,1), Np);
        x(1:Np,:)   = x(ind,:);
%         x = sortrows(x,1, "ascend");
        best    = 1; %indentify the best solution
        %Set to work just the Np best solutions----------------------------

        %if algorithm has reached maxFE break the loop---------------------
        if(break_alg)
            break;
        end
        %if algorithm has reached maxFE break the loop---------------------
        
        %DE over each population element-----------------------------------
        for i=1:Np
            %Sample the necessary elements---------------------------------
%             ind = 1:Np;
%             ind = randsample(ind, 5); %matlab function is slow
            ind1(1)  = i;
            ind1     = rand_sample(Np, 5, ind1);
            %Sample the necessary elements---------------------------------
            
            %Ensemble over variable----------------------------------------
            for j=1:D
                F = unifrnd(0.1, 1.5, 1);
                mut     = randi(5);
                if mut == 1 %rand/1
                        x(Np+1,j+1)   = x(ind1(2),1+j) + F*(x(ind1(3),1+j) - x(ind1(4),1+j));
                elseif mut == 2 %best/1
                        x(Np+1,j+1)   = x(best,1+j) + F*(x(ind1(2),1+j) - x(ind1(3),1+j));
                elseif mut == 3 %target-to-best/1
                        x(Np+1,j+1)   = x(i,1+j) + F*(x(best,1+j) - x(i,1+j)) + F*(x(ind1(2),1+j) - x(ind1(3),1+j));
                elseif mut == 4 %rand/2
                        x(Np+1,j+1)   = x(ind1(2),1+j) + F*(x(ind1(3),1+j) - x(ind1(4),1+j)) + F*(x(ind1(5),1+j) - x(ind1(6),1+j));
                else %best/2
                        x(Np+1,j+1)   = x(best,1+j) + F*(x(ind1(2),1+j) - x(ind1(3),1+j)) + F*(x(ind1(4),1+j) - x(ind1(5),1+j));
                end

                while(x(Np+1,j+1)>upLim)
                    x(Np+1,j+1) = x(Np+1,j+1) - (upLim-lowLim);
                end
                while(x(Np+1,j+1)<lowLim)
                    x(Np+1,j+1) = x(Np+1,j+1) + (upLim-lowLim);
                end
            end
            %Ensemble over variable----------------------------------------
            
            %binary crossover----------------------------------------------
            k   = randi(D);
            for j=1:D
                if(rand >= CR && j~=k)
                    x(Np+1,j+1)    = x(i,1+j);
                end
            end
            %binary crossover----------------------------------------------

            %competition---------------------------------------------------
            x(Np+1,1)=fitt(x(Np+1,2:D+1)');
            totEval = totEval + 1;

            if(x(i,1) >= x(Np+1,1))
                x(i,:)  = x(Np+1,:);
                %update best position--------------------------------------
                if i~=best && x(i,1) < x(best,1)
                    best = i;
                end
                %update best position--------------------------------------
            end
            %competition---------------------------------------------------

            %save best solution
%             if mod(totEval,smpl)==0
%                 plt(indPlt) = x(best,1);
%                 indPlt  = indPlt + 1;
%             end
            
            %break loop if max FE is reached-------------------------------
            if totEval >= maxEval
                break;
            end
            %break loop if max FE is reached-------------------------------
        end
        %DE over each population element-----------------------------------
    end
    elite   = x(best, 2:D+1)';
    fittEle = x(best, 1);
end






