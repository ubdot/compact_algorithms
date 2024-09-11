%% Example of calling the functions of compact algorithms
clear
% clc
param.F         = 0.5;
param.D         = 20;
param.CR        = 1/2^(1/(0.25*param.D));
param.NP        = 300;
param.lowLim    = -100*ones(param.D,1);
param.upLim     = 100*ones(param.D,1);
param.maxEval   = 5000*param.D;
param.cross_op  = 2; %used for the cDE algorithm, exp crossover
param.pe        = true;
param.age       = 5000;
param.pa        = 0.25;
param.rho       = 0.01;
%%
rng(1,"twister"); %Random seed
n_max=10;
fx=zeros(n_max,1);
x=zeros(n_max, param.D);
minFx = Inf;
%% Calling function 

for i=1:n_max
%     [x(i,:),fx(i)]=cDE(@fitt, param);
%     [x(i,:),fx(i)]=cPSO(@fitt, param);
    [x(i,:),fx(i)]=rcGA(@fitt, param);
%     disp([i, fx(i)]);
end
%% Get stats
meanfx=mean(fx);
maxfx=max(fx);
minfx=min(fx);
stdfx=std(fx);

disp(["Mean", "Min", "Max", "STD"])
disp([meanfx, minfx, maxfx, stdfx])

%% Fittness function
function fx=fitt(x)
%Modify the bounds accordingly the function
    %sphere [-100, 100]
%     fx=x'*x; 

    %Schwefel 2.22 [-10, 10]
%     fx=sum(abs(x))+prod(abs(x));

    %Schwefel 1.2 [-100, 100]
%     fx=0;
%     for i=1:30
%     fx=fx+sum(x(1:i))^2;
%     end

    %Schwefel 2.21 [-100, 100]
%     fx=max(abs(x));

    %Rosenbrock [-30,30]
%     fx=sum(100*(x(2:30)-(x(1:30-1).^2)).^2+(x(1:30-1)-1).^2);

    %Schwefel 2.26 [-500, 500]
% fx=1/30*sum(-x.*sin(sqrt(abs(x))));

    %Rastrigin [-5.12, 5.12]
fx=sum(x.^2-10*cos(2*pi.*x))+10*30; 

    %Ackley1 [-35, 35]
%     fx=-20*exp(-.2*sqrt(sum(x.^2)/30))-exp(sum(cos(2*pi.*x))/30)+20+exp(1);

    %Griewank [-100, 100]
%     fx=sum(x.^2)/4000-prod(cos(x./sqrt([1:30]')))+1; 

end