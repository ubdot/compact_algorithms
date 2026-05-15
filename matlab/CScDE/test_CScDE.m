%% Create pool to make the calculations in parallel
% parpool(5)

%% Number of executions according to CEC 2017
% num_exe  = 51;
% tot_prob = 30;
res104   = zeros(num_exe,tot_prob);
% rng(1,"twister");

%% Used parameters with CScDE
param.D         = D; % This parameter is modified(manual) each run
param.lowLim    = -100*ones(param.D,1);
param.upLim     = 100*ones(param.D,1);
param.NP        = 300;
param.maxEval   = param.D*10000;
param.nWaves    = ceil(log(param.D)/log(2));
param.freq      = 1/param.D;
disp('running CScDE with parameters:')
disp(param);


%% Algorithm run

% plt2=zeros(1001,58);
for j54 = 1:tot_prob
    fitt = @(x) cec17_func(x,j54);
    for i54=1:num_exe %parfor i54=1:num_exe
        [~,min_Fitt]  = CScDE(fitt, param);
        res104(i54, j54) = min_Fitt;
%         [~,~, plt]  = CScDE(fitt, param);
%         plt2(:,i54) = plt;
    end
    disp(j54);
end

%% Save file
% filename = sprintf("cec_17/results/D%d_CScDE_final.mat",param.D);
% save(filename, "res104");

%%
% delete(gcp("nocreate"));