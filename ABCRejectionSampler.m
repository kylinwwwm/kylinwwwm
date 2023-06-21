function [result,flag] = ABCRejectionSampler(N,prior,f,rho,epsilon,T,gene)
%% Rejection Sampler for approximate Bayesian computaion
% Inputs:
%    N - the number of ABC posterior samples ABC后验样本的数量 1000
%    prior - function that generates iid samples from the parameter joint
%        prior distribution 从参数联合先验分布生成独立同分布样本的函数
%    f - function that computes statictis given a parameters set
%    给定参数集计算汇总统计量的函数
%    rho - discrepancy metric, treated as a function of simulated data only 差异指标，仅视为模拟数据的函数
%    epsilon - the discrepancy acceptance threshold 容差接受阈值
%    T - number of rounds轮数   5
%    gene - sequence number of the gene.
%
% Outputs:
%    result - a matrix of ABC prior samples 一个ABC先验样本的矩阵

result = cell(N,T);%生成维度为N*T的元胞数组，cell是matlab中的一种数据类型：元胞
result0 = [];
total_time = 0;
flag = true;
while total_time < N && size(result0,2) < 10*N % A && B首先判断A,若A为假则不必再判断B,size(result0,2)返回result0的列数
    tic;
    % generate trial from the prior
    theta_trial = prior(); % 从给定的先验分布中抽取初始参数，初始参数向量为theta_trial
    param.kon1 = theta_trial(1);
    param.kon2 = theta_trial(2);
    param.ron1 = theta_trial(3);
    param.ron2 = theta_trial(4);
    param.koff = theta_trial(5);
    param.roff = theta_trial(6);
    param.mu = theta_trial(7);
    param.delta = theta_trial(8);
    param.q = theta_trial(9);
    % compute theorical statictis of parameters
    static_theo = f(param);%f为在给定参数集的情况下计算汇总统计量
    dist = rho(static_theo);%rho，差异指标，仅仅被视为模拟数据的函数（距离度量）
    param.dist = dist;
    % accept or reject
    if dist <= epsilon
        ind = size(result0,2)+1;%size(result0,2)返回result0第2维的尺寸
        result0(ind).kon1 = param.kon1;
        result0(ind).kon2 = param.kon2;
        result0(ind).ron1 = param.ron1;
        result0(ind).ron2 = param.ron2;
        result0(ind).koff = param.koff;
        result0(ind).roff = param.roff;
        result0(ind).mu = param.mu;
        result0(ind).delta = param.delta;
        result0(ind).q = param.q;
        result0(ind).dist = param.dist;
      
%         disp(ind)
    end
    elapsedTime = toc;%计算程序完成时间
    total_time = total_time + elapsedTime;
end

if total_time > 300
    flag = false;
    fprintf('Gene %d :wrong!!\n',gene);
    pause(5) %pause(5)表示matlab暂停运行5s
else
    [~,index0] = sort([result0.dist]);%sort函数默认矩阵按列升序，~表示排序过后的矩阵，此处省略不用输出。index0表示排序后元素在原始矩阵中的位置,所以~中元素越往后，代表抽到的参数产生的数据与真实数据差距越大，越不可取，因此result0只取前1到N个
    result0 = result0(index0(1:N));%因为后面需要选取上一次偏差的中位数作为下一次容差的阈值，因此需要对容差进行排序
    for index2 = 1:N
        result(index2,1) = {result0(index2)};%result的第一列为result0中的元素，并且第一个元素参数产生的数据与真实数据差别最小
    end
end
end  %这一部分执行了从分布当中采得第一批差异小于初始阈值得样本，并对第一批样本产生得数据与真实数据之间的差值进行排序，取前N个差值对应的样本粒子作为第1次采样的结果
