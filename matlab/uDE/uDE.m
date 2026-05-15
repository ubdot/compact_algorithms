function [elite,fittEle] = uDE(param,fitt)
% function [elite,fittEle, plt] = uDE(param,fitt)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Load parameters to local variables
Npop    = param.pop;
D       = param.D;
upLim   = param.upLim;
lowLim  = param.lowLim;
maxEval = param.maxEval;
F       = param.F;
resLim  = param.resLim;
CR      = param.CR;

%initialized used variables
pop     = zeros(D, Npop);
Fpop    = zeros(Npop, 1);
trial   = zeros(D,1);
totEval = 0;
elite   = 1;

%used for the sampling function
ind=zeros(4,1);

%used to break the cycles in case that the max number of FE is reached
breakLoop   = false;

%array to save the best positions
% smpl    = 10*D;
% plt     = zeros(1001,1);
% indPlt  = 1;

%The population is initialized
for i=1:Npop
    pop(:,i)    = unifrnd(lowLim, upLim,  [D, 1]);
    Fpop(i)     = fitt(pop(:,i));
    totEval     = totEval + 1;
    %Identification of elite position
    if(i ~= elite && Fpop(i)<Fpop(elite))
        elite   = i;
    end
end
% plt(indPlt)    = Fpop(elite);
% indPlt  = indPlt + 1;

while(totEval < maxEval)
    %DE on the population, counting the generations number
    for gen = 1:resLim
        %For each member of the population
        for i=1:Npop
            %ind =1:Np;
            %ind(i)=[];
            %ind=randsample(ind,3); %matlab function is slow
            %Sample individuals--------------------------------------------
            ind(1)  = i; %This is to prevent sample the same element in the pop
            ind     = rand_sample(Npop, 3, ind);
            k   = randi(D);
            %Sample individuals--------------------------------------------

            %Binary crossover----------------------------------------------
            for j=1:D
                if(rand < CR || j==k)
                    %Mutant vector-----------------------------------------
                    trial(j)    = pop(j, ind(2)) + F*(pop(j, ind(3)) - pop(j, ind(4)));
%                     trial(j)    = pop(j, ind(1)) + F*(pop(j, ind(1)) - pop(j, ind(1)));% Error in the first code
                    %Mutant vector-----------------------------------------
                    %Bounds fix--------------------------------------------
                    while(trial(j)>upLim)
                        trial(j) = trial(j) - (upLim-lowLim);
                    end
                    while(trial(j)<lowLim)
                        trial(j) = trial(j) + (upLim-lowLim);
                    end
                    %Bounds fix--------------------------------------------
                else
                    trial(j)    = pop(j, i);
                end
            end
            %Binary crossover----------------------------------------------
            
            %competition---------------------------------------------------
            Ftrial  = fitt(trial);
            totEval = totEval + 1;
            if(Ftrial < Fpop(i))
                pop(:,i)    = trial;
                Fpop(i)     = Ftrial;
                %elite position update-------------------------------------
                if(i ~= elite && Fpop(i)<Fpop(elite))
                    elite   = i;
                end
                %elite position update-------------------------------------
            end
            %competition---------------------------------------------------
%             %Save the best solution in array
%             if mod(totEval,smpl)==0
%                 plt(indPlt) = Fpop(elite);
%                 indPlt  = indPlt + 1;
%             end
            %if max number of FE has been reached break the loop-----------
            if totEval >= maxEval
                breakLoop   = true;
                break;
            end
            %if max number of FE has been reached break the loop-----------
        end
        %if max number of FE has been reached break the loop---------------
        if breakLoop
            break;
        end
        %if max number of FE has been reached break the loop---------------
    end
    %restart all the population but the elite------------------------------
    if ~breakLoop
        for i=1:Npop
            if(i==elite)
                continue;
            end
            pop(:,i)    = unifrnd(lowLim, upLim,  [D, 1]);
            Fpop(i)     = fitt(pop(:,i));
            totEval     = totEval + 1;
            if(Fpop(i)<Fpop(elite))
                elite   = i;
            end

%             %Save the best solution in array
%             if mod(totEval,smpl)==0
%                 plt(indPlt) = Fpop(elite);
%                 indPlt  = indPlt + 1;
%             end
            if totEval >= maxEval
                break;
            end
        end
    end
    %restart all the population but the elite------------------------------
end
% disp([totEval,indPlt]);
fittEle = Fpop(elite);
elite   = pop(:,elite);
end