%% Create pool to make the calculations in parallel
% parpool(5)

%% Number of executions according to CEC 2017
% num_exe  = 51;
% tot_prob = 30;
% rng(1,"twister");
res104   = zeros(num_exe,tot_prob);
%% Used parameters with cPSO
param.c0        = -0.2;
param.c1        = -0.07;
param.c2        = 3.74;
param.NP        = 300;
param.D         = D; % This parameter is modified(manual) each run
param.maxEval   = param.D*10000;
param.lowLim    = -100;
param.upLim     = 100;
disp('running cPSO with parameters:')
disp(param);

%% Algorithm run
% plt2=zeros(1001,58);
for j54 = 1:tot_prob
    fitt = @(x) cec17_func(x,j54);
    for i54=1:num_exe %parfor i54=1:num_exe
        [~,min_Fitt]=cPSO(fitt, param);
%         [~,~,plt]=cPSO(fitt, param);
        res104(i54, j54) = min_Fitt;
    end
    disp(j54);
end

%% Save file
% filename = sprintf("cec_17/results/D%d_cPSO_final.mat",param.D);
% save(filename, "res104");

%%
% delete(gcp("nocreate"));