param.A     = 1;
param.P     = 10;
param.B     = 2;
param.S     = 4;
param.eps   = 1e-5;
param.H     = 30;
param.D     = 100;
param.maxEval   = param.D*10000;
param.lowLim= -100;
param.upLim = 100;
%%
num_exe  = 51;
tot_prob = 30;
res104   = zeros(num_exe,tot_prob);
rng(1,"twister");
% parpool(5)
for j54 = 1:tot_prob
    fitt = @(x) cec17_func(x,j54);
    parfor i54=1:num_exe
        [~,min_Fitt]=ISPO(fitt, param);
        res104(i54, j54) = min_Fitt;
    end
    disp(j54);
end

%%
filename = sprintf("CEC_2017_comparison_results/D%d_ISPO_final.mat",param.D);
save(filename, "res104");