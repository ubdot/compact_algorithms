clc;
param.D         = 100;
param.lowLim    = -100*ones(param.D,1);
param.upLim     = 100*ones(param.D,1);
param.NP        = 1000;
param.maxEval   = param.D*10000;
param.nWaves    = floor(log(param.D)/log(2));
param.freq      = 1/param.D;

fitt = @(x) cec17_func(x,1);
% [~,fx]  = CScDE(fitt, param);

disp(param);
%%
num_exe  = 51;
tot_prob = 30;
res104   = zeros(num_exe,tot_prob);
rng(1,"twister");
% parpool(5)
for j54 = 19:19
    fitt = @(x) cec17_func(x,j54);
    parfor i54=1:num_exe
        [~,min_Fitt]  = CScDE(fitt, param);
        res104(i54, j54) = min_Fitt;
    end
    disp(j54);
end

%%
% filename = sprintf("CEC_2017_comparison_results/final/D%d_CScDE_final.mat",param.D);
% save(filename, "res104");