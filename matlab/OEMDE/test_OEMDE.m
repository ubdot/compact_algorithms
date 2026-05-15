%% Create pool to make the calculations in parallel
% parpool(5)

%% Number of executions according to CEC 2017
% num_exe  = 51;
% tot_prob = 30;
res104   = zeros(num_exe,tot_prob);
% rng(1,"twister");
%% Used parameters with OEMDE
param.pop   = 6;
param.D     = D;
param.upLim = 100;
param.lowLim= -100;
param.CR    = 0.9;
param.maxEval   = 10000*param.D;
disp('running OEMDE with parameters:')
disp(param);

%% Algorithm run

% plt2=zeros(1001,58);
for prob=3:3
    fitt=@(x) cec17_func(x,prob);
    for j=1:num_exe
        [~,min_Fitt]    = OEMDE(param,fitt);
        res104(j,prob)  = min_Fitt;
%         [~,~,plt] = OEMDE(param,fitt);
%         plt2(:,j) = plt;
    end
    disp(prob)
end

%% Save file
% filename = sprintf("CEC_2017_comparison_results/final/D%d_OEMDE.mat",param.D);
% save(filename, "res104");

%% Delete parpool
% delete(gcp("nocreate"));