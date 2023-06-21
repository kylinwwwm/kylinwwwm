%% master program for MEFs scRNA-seq data
clear;clc;
addpath(genpath(pwd))

tic;
if isempty(gcp('nocreate'))
    parpool(4);
end
toc;
%% load preprocessed data加载预处理的数据
data_all = csvread('data/MEF_QC_all.csv');
save_folder  = fullfile(pwd,'results1/results_MEF');

% check the uninferred gene检查未推断的基因
gene_number_wait = [];
for gene_number = 1:size(data_all,1)%data_all的行数  即基因数
    filename = sprintf('//result_gene_%d.mat',gene_number);%文件名 字符串格式
    if exist([save_folder,filename])==0
        gene_number_wait = [gene_number_wait,gene_number];
    end
end
fprintf('开始运行了\n');

%% main program
parfor infer_index = 1:length(gene_number_wait)
    % Delete 5% of the tail data 删除百分之五的尾部数据
    gene = gene_number_wait(infer_index);
    data = data_all(gene,:);
    data = data(data>=0);
    [~,inter] = mink(data,floor(0.95*length(data)));%[B,I] = mink(A,k) 计算 A 的 k 个最小值的索引，并在 I 中返回这些索引
    data = data(inter);
    data_mean = mean(data);
    data_var = var(data);
    data_noise = data_var/data_mean^2;
    
    % inference
    statis_data = statisData(data);
    rho = @(s) sqrt(sum(log(statis_data./s).^2));%s表示模拟数据，statis_data
    f = @(k) statisQM(k,4);
    N = 1000;
    T = 5;
    epsilon = 1;
    % [kon ron koff roff mu delta]
    % OFF dwell time distribution f_off(t) = ron^(kon) * t^(kon-1) * e^(-ron * t) / gamma(kon)
    % ON dwell time distribution f_on(t) = roff^(koff) * t^(koff-1) * e^(-roff * t) / gamma(koff)
    prior = @() [5*rand(),5*rand(),logunif(-1,1),logunif(-1,1),5*rand(),logunif(-1,1),50*rand,1,rand];%koff,roff,kon,ron,rsyn的先验分布，分别是均匀分布U(0,5)，基于10的区间[-1，1]的对数均匀分布...
    proposal_sigma = 0.2;
    proposal = @(x) lognrnd(x,proposal_sigma);%lognrnd(mu,sigma)生成均值为mu，标准差为sigma的随机数 转移核选取对数正态分布  扰动
    proposal_pdf = @(kon1_post,kon1_prior,kon2_post,kon2_prior,ron1_post,ron1_prior,ron2_post,ron2_prior,koff_post,koff_prior,roff_post,roff_prior,mu_post,mu_prior,q_post,q_prior)...
        lognpdf(mu_post,log(mu_prior),proposal_sigma) * lognpdf(kon1_post,log(kon1_prior),proposal_sigma) *...
        lognpdf(kon2_post,log(kon2_prior),proposal_sigma) * lognpdf(ron1_post,log(ron1_prior),proposal_sigma) *...
        lognpdf(ron2_post,log(ron2_prior),proposal_sigma) * lognpdf(koff_post,log(koff_prior),proposal_sigma) *...
        lognpdf(roff_post,log(roff_prior),proposal_sigma) * lognpdf(q_post,log(q_prior),proposal_sigma);
    [result,~] = ABCPRCSampler(N,prior,f,rho,epsilon,T,proposal,proposal_pdf,1);
    
    % save result
    filename = sprintf('//result_gene_%d',gene);
    parsave([save_folder,filename],gene,data,result);
    fprintf('模拟已接受基因%d的结果\n',gene);
end
load('E:\code\1_test\results1\results_MEF\result_gene_1')