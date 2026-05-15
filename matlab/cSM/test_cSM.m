%% Create pool to make the calculations in parallel
% parpool(5)

%% Number of executions according to CEC 2017
% num_exe  = 51;
% tot_prob = 30;
res104   = zeros(num_exe,tot_prob);
% rng(1,"twister");

%%  Used parameters with cSM
param.NP        = 2; %In our tests show better results than 1
param.D         = D;
param.B         = 6;
param.maxEval   = param.D*10000;
param.lowLim    = -100;
param.upLim     = 100;
disp('running cSM with parameters:')
disp(param);

%% Algorithm run

% plt2=zeros(1001,58);
for j54 = 1:tot_prob
    fitt = @(x) cec17_func(x,j54);
    for i54=1:num_exe %parfor i54=1:num_exe
        [~,min_Fitt]=CSM(fitt, param);
        res104(i54, j54) = min_Fitt;
%         [~,~,plt]=CSM(fitt, param);
%         plt2(:,i54) = plt;
    end
    disp(j54);
end

%% Save file
% filename = sprintf("cec_17/results/D%d_cSM_final.mat",param.D);
% save(filename, "res104");

%%
% delete(gcp("nocreate"));