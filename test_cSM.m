param.NP        = 2;
param.D         = 30;
param.B         = 6;
param.maxEval   = param.D*10000;
param.lowLim    = -100;
param.upLim     = 100;
disp(param);
%%
num_exe  = 51;
tot_prob = 30;
res104   = zeros(num_exe,tot_prob);
rng(1,"twister");
% parpool(5)
for j54 = 6:6
    fitt = @(x) cec17_func(x,j54);
    parfor i54=1:num_exe
        [~,min_Fitt]=CSM(fitt, param);
        res104(i54, j54) = min_Fitt;
    end
    disp(j54);
end

%%
% filename = sprintf("CEC_2017_comparison_results/final/D%d_cSM_final.mat",param.D);
% save(filename, "res104");