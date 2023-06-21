param_true.kon1 = 3;
param_true.ron1 = 0.5;
param_true.kon2 = 6;
param_true.ron2 = 0.7;
param_true.koff = 5;
param_true.roff = 0.5;
param_true.mu = 20;
param_true.delta = 1;
param_true.x0 = [1,0,0];
param_true.tottime = 2000;
param_true.q = 0.4;
% param_true.kon1 = 1;
% param_true.ron1 = 0.8;
% param_true.kon2 = 3;
% param_true.ron2 = 0.7;
% param_true.koff = 5;
% param_true.roff = 3;
% param_true.mu = 30;
% param_true.delta = 1;
% param_true.x0 = [1,0,0];
% param_true.tottime = 2000;
% param_true.q = 0.6;

% param_true.kon1 = 5;
% param_true.ron1 = 0.4;
% param_true.kon2 = 6;
% param_true.ron2 = 0.6;
% param_true.koff = 3;
% param_true.roff = 0.5;
% param_true.mu = 10;
% param_true.delta = 1;
% param_true.x0 = [1,0,0];
% param_true.tottime = 2000; 
% param_true.q = 0.6;


% param_true.kon1 = 3;
% param_true.ron1 = 0.4;
% param_true.kon2 = 6;
% param_true.ron2 = 0.7;
% param_true.koff = 3;
% param_true.roff = 0.6;
% param_true.mu = 40;
% param_true.delta = 1;
% param_true.x0 = [1,0,0];
% param_true.tottime = 2000; 
% param_true.q = 0.7;


% Simulation algorithm for GTM.
[x,t] = simulQM(param_true);
tq = 1500:0.1:param_true.tottime;
xq = interp1(t,x(:,3),tq,'previous');%插值  函数关系由t,x(:,3)给出，插值tq对应的xq
data = xq;
statis_data = statisData(data);%%%%用参数产生的模拟数据得到的汇总统计量，在这里相当于真实数据
statis_ther = statisQM(param_true,4);%%用提前设定好的真实参数计算的汇总统计量
%eps0 = sqrt(sum(statis_data./statis_ther).^2);
eps0 = sqrt(sum(log(statis_data./statis_ther).^2));
N = 1000;
prior = @() [5*rand(),5*rand(),logunif(-1,1),logunif(-1,1),5*rand(),logunif(-1,1),50*rand,1,rand];
f = @(k) statisQM(k,4);%k是param   用推断出的参数计算汇总统计量
%rho = @(s) sqrt(sum(statis_data./s).^2);%s是由f计算得到的汇总统计量
rho = @(s) sqrt(sum(log(statis_data./s).^2));
T = 6;
epsilon = 1;
%[result,flag] = ABCRejectionSampler(N,prior,f,rho,epsilon,T,1);
proposal_sigma = 0.2;
proposal = @(x) lognrnd(x,proposal_sigma);
proposal_pdf = @(kon1_post,kon1_prior,kon2_post,kon2_prior,ron1_post,ron1_prior,ron2_post,ron2_prior,koff_post,koff_prior,roff_post,roff_prior,mu_post,mu_prior,q_post,q_prior)...
    lognpdf(mu_post,log(mu_prior),proposal_sigma) * lognpdf(kon1_post,log(kon1_prior),proposal_sigma) *...
    lognpdf(kon2_post,log(kon2_prior),proposal_sigma) * lognpdf(ron1_post,log(ron1_prior),proposal_sigma) *...
    lognpdf(ron2_post,log(ron2_prior),proposal_sigma) * lognpdf(koff_post,log(koff_prior),proposal_sigma) *...
    lognpdf(roff_post,log(roff_prior),proposal_sigma) * lognpdf(q_post,log(q_prior),proposal_sigma);

[result,~] = ABCPRCSampler(N,prior,f,rho,epsilon,T,proposal,proposal_pdf,1);
figureResult(data,result(:,end),param_true);%取出了result的最后一列
%    N - the number of ABC posterior samples ABC后验样本的数量 result(:,end)取result元素的最后一列
%    prior - function that generates iid samples from the parameter joint
%        prior distribution 从参数联合先验分布生成独立同分布样本的函数
%    f - function that computes statictis given a parameters set
%    给定参数集计算汇总统计量的函数
%    rho - discrepancy metric, treated as a function of simulated data only 差异指标，仅视为模拟数据的函数
%    epsilon - the discrepancy acceptance threshold 容差接受阈值
%    T - number of rounds
%    gene - sequence number of the gene.
%save(sprintf('E:/code/5-BayesGTM-main/BayesGTM-main/results1/example/%d_%.1f_%d_%.1f_%d_%d.mat',param_true.kon,param_true.ron,param_true.koff,param_true.roff,param_true.mu))