function [result,flag] = ABCPRCSampler(N,prior,f,rho,epsilon,T,proposal,proposal_pdf,gene)
%% Sequential Monte Carlo Sampler for approximate Bayesian computaion
% Inputs:
%    N - the number of particles 粒子数
%    prior - prior distribution sampler 先验分布采样器
%    f - function that generates simulated data give a parameters set
%    给定参数集生成模拟数据的函数
%    rho - discrepancy metric, treated as a function of simulated data
%    only距离度量，仅是模拟数据的函数
%    epsilon - a sequence of discrepancy acceptance thresholds 容差接受阈值序列
%    T - number of rounds
%    proposal - proposal kernel process, generates new trials 建议核过程，生成新的试验
%    proposal_pdf - the PDF of the proposal kernel建议核概率密度函数
%    gene - sequence number of the gene.
%
% Outputs:
%    result - the sequence of particles一系列粒子

% initialise
[result,flag] = ABCRejectionSampler(N,prior,f,rho,epsilon,T,gene);
W = (1/N)*ones(N,T);

% sequential sampling
if flag == true
    for t = 2:T
        epsilon = prctile(cellfun(@(c)c.dist,result(:,t-1)),50);%容差阈值的选择策略，对于第t次迭代，将接受阈值ε_t设置为从第t−1迭代获得的结果的偏差的中位数
%         fprintf('The eps of %d rounds is %f\n',t,epsilon);
        % generate new generation of particles生成新一代粒子
        for i = 1:N
            % rejections steps
            r = inf;% 无穷大
            total_time = 0;
            q = inf;
           %while total_time < 20 && r > epsilon % 只要满足该条件，就会一直循环，只有不满足该条件时也就是说只有r <= epsilon时才会跳出循环
           while total_time < 20 && r > epsilon || (q < 0 ||  q > 1)
                tic; % tic计时单位是s
                % sample from weighted particles
                j = randsample(N,1,true,W(:,t-1));
                %y =randsample(n,k,true,w)从1：n中以权重w抽k个值，相当于抽的是t-1次抽样群体中某一个样本的索引
                param_temp = result{j,t-1} ;%抽到的就是第t-1次的抽样结果中第j个样本，实值是result{j,t-1}
                
                % generate new particle based on proposal kernel基于建议核生成新的粒子
                %  根据建议分布进行一个扰动
                %param_proposal = proposal(log([param_temp.kon1,param_temp.kon2,param_temp.ron1,param_temp.ron2,param_temp.koff,param_temp.roff,param_temp.mu]));%param_proposal = proposal(log([param_temp.kon1,param_temp.kon2,param_temp.ron1,param_temp.ron2,param_temp.koff,param_temp.roff,param_temp.mu,param_q]));
                % q_proposal = proposal(log(param_temp.q))
                param_proposal = proposal(log([param_temp.kon1,param_temp.kon2,param_temp.ron1,param_temp.ron2,param_temp.koff,param_temp.roff,param_temp.mu,param_temp.q]));
                param.kon1 = param_proposal(1);
                param.kon2 = param_proposal(2);
                param.ron1 = param_proposal(3);
                param.ron2 = param_proposal(4);
                param.koff = param_proposal(5);
                param.roff = param_proposal(6);
                param.mu = param_proposal(7);
                param.q = param_proposal(8);
                q = param.q;
                param.delta = 1;
                result{i,t} = param;
                static_temp = f(param);%用推断出的参数直接计算汇总统计量
                r = rho(static_temp);% rho是距离度量
                result{i,t}.dist = r;
                elapsedTime = toc;
                total_time = total_time + elapsedTime;
           end

            if total_time > 20
                flag = false;
                break;
            end
            % recompute particle weight using optimal backward kernel
            % 使用最佳后向核重新计算粒子权重
            K_t = proposal_pdf(result{i,t}.kon1,result{j,t-1}.kon1,...
                    result{i,t}.kon2,result{j,t-1}.kon2,...
                    result{i,t}.ron1,result{j,t-1}.ron1,...
                    result{i,t}.ron2,result{j,t-1}.ron2,...
                    result{i,t}.koff,result{j,t-1}.koff,...
                    result{i,t}.roff,result{j,t-1}.roff,...
                    result{i,t}.mu,result{j,t-1}.mu,...
                    result{i,t}.q,result{j,t-1}.q);%result{j,t-1}.q
            L_t_1 = proposal_pdf(result{j,t-1}.kon1,result{i,t}.kon1,...
                    result{j,t-1}.kon2,result{i,t}.kon2,...
                    result{j,t-1}.ron1,result{i,t}.ron1,...
                    result{j,t-1}.ron2,result{i,t}.ron2,...
                    result{j,t-1}.koff,result{i,t}.koff,...
                    result{j,t-1}.roff,result{i,t}.roff,...
                    result{j,t-1}.mu,result{i,t}.mu,...
                    result{j,t-1}.q,result{i,t}.q);%result{j,t-1}.q,result{i,t}.q
            pi_t = 1/5*1/5*1/5*1/50*(result{i,t}.ron1 * result{i,t}.ron2 * result{i,t}.roff);
            pi_t_1 = 1/5*1/5*1/5*1/50*(result{j,t-1}.ron1 * result{j,t-1}.ron2 * result{j,t-1}.roff);
                
            
            W(i,t) = (pi_t * L_t_1)/ (K_t * pi_t_1);
        end
%        W
        if flag == false
            break;
        end
        %采样得到一组N*1的粒子之后，通过转移核给每个粒子分配权重，再根据权重对这些粒子进行重采样，重采样得到的新的N*1的粒子作为第t次抽样的结果。下一次采样的初始权重还是1/N，采样之后分配权重并进行重采样
        % resample
        if t < T
            result_rs  = result(:,t);%表示result的第t列
            W(:,t) = W(:,t)./sum(W(:,t));%好像是按列归一化
            J = randsample(N,N,true,W(:,t));%重采样得到的样本已经经过判断满足容差阈值，且不经过扰动
            result(:,t) = result_rs(J);
        end
        
        %   re-set weights
        W(:,t) = 1/N;
    end
end
