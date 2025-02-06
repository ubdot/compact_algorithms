%% Create pool to make the calculations in parallel
% parpool(5)
%% Number of executions according to CEC 2017
% num_exe  = 51;
% tot_prob = 30;
res104   = zeros(tot_prob,num_exe);
% rng(1,"twister");
%% Parameters used with uDE
param.pop       = 5;
param.D         = D;
param.upLim     = 100;
param.lowLim    =-100;
param.maxEval   = param.D*10000;
param.F         = 0.89;
param.resLim    = 90*param.D;
param.CR        = 0.46;

disp('running uDE with parameters:')
disp(param);

%% Algorithm run
% plt2=zeros(1001,58);
for prob=1:tot_prob
    fitt=@(x) cec17_func(x,prob);
    for j=1:num_exe
        [~,min_Fitt] = uDE(param,fitt);
        res104(prob,j) = min_Fitt;
%         [~,~,plt] = uDE(param,fitt);
%         plt2(:,j) = plt;
    end
    disp(prob)
end

%% Save file results
% filename = sprintf("CEC_2017_comparison_results/final/D%d_uDE.mat",param.D);
% save(filename, "res104");
