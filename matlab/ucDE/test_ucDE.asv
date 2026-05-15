%% Create pool to make the calculations in parallel
% parpool(5)
%% Number of executions according to CEC 2017
% num_exe  = 51;
% tot_prob = 30;
res104  =zeros(tot_prob,num_exe);
% rng(1,"twister");
%% Parameters used with ucDE
F           = 0.89;%0.8962;
CR_bin      = 0.46;%0.4621;
CR_exp      = 0.21;%0.2127;
upop_budget = 0.90;%0.90;
num_upp     = 90;%90;%97D executions
Np          = 1;%1;%10 for 10D, 30 for 30D, 50 for 50D, and 100 for 100D
% D           = 10;

fprintf("running ucDE with parameters:\nF:\t\t\t\t%g\nCR_bin:\t\t\t%g\nCR_exp:\t\t\t%g\nupop_budget:\t%g\nnum_upp:\t\t%d\nNp:\t\t\t\t%d\nD:\t\t\t\t%d\n", F,CR_bin,CR_exp,upop_budget,num_upp,Np,D);
%% Algorithm run
% plt2=zeros(1001,58);
% parpool(5)
for prob=1:tot_prob
    SEED= randi(2^32-1, [num_exe,1]);
    parfor j=1:num_exe
        command = sprintf('"ucDE_endCom2.exe" %d %f %f %f %f %d %d %d %d',prob, F, CR_bin, CR_exp, upop_budget, num_upp, Np, SEED(j), D);
        [~,temp]=system(command);
        min_Fitt=sscanf(temp,'%e');
        res104(prob,j) = min_Fitt;
    end
    disp(prob);
end

%% Save file results
% filename = sprintf("CEC_2017_comparison_results/final/D%d_ucDE_endCom.mat",D);
% save(filename, "res104");
% %%
% disp(datetime("now"))

