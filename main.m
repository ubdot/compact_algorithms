%%
%util functions
addpath(genpath('cec_17'));
%must compile file inside the folder
%use command
%"mex cec17_func.cpp -DWINDOWS"
addpath(genpath('util'));

%Algorithms
addpath(genpath('cSM'));
addpath(genpath('cPSO'));
addpath(genpath('CScDE'));
addpath(genpath('ISPO'));
addpath(genpath('cDE'));
addpath(genpath('rcGA'));
addpath(genpath('OEMDE'));
addpath(genpath('ucDE'));
addpath(genpath('uDE'));

%%
num_exe  = 1;
tot_prob = 30;
D        = 10;
rng(1,"twister");

%% CPSO
test_cPSO

%% CScDE
test_CScDE

%% cSM
test_cSM

%% ISPO
test_ispo

%% OEuDE
test_OEMDE

%% cDE
test_cDE

%% uDE
test_uDE

%% ucDE
%must be compiled and the exe must be in the root file
test_ucDE
