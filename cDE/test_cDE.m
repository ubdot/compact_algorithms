%% Create pool to make the calculations in parallel
% parpool(5)
%% Number of executions according to CEC 2017
% num_exe  = 51;
% tot_prob = 30;
res104   = zeros(num_exe,tot_prob);
% rng(1,"twister");
%% Parameters used with cDE
param.D         = 10;
param.maxEval   = param.D*10000;
param.lowLim    = -100;
param.upLim     = 100;
param.F         = 0.5;
param.NP        = 300;
param.CR        = 1/2^(1/(0.25*param.D));
param.cross_op  = 2; %used for the cDE algorithm, exp crossover

disp('running cDE-exp with parameters:')
disp(param);
%% Algorithm run
% plt2=zeros(1001,58);
for j54 = 1:tot_prob
    fitt = @(x) cec17_func(x,j54);
    for i54=1:num_exe %parfor i54=1:num_exe
        [~,min_Fitt]=cDE(fitt, param);
        res104(i54, j54) = min_Fitt;
%         [~,~, plt]=cDE(fitt, param);
%         plt2(:,i54) = plt;
    end
    disp(j54);
end

%% Save file results
% filename = sprintf("CEC_2017_comparison_results/final/D%d_cDE_final.mat",param.D);
% save(filename, "res104");
% % delete(gcp("nocreate"));