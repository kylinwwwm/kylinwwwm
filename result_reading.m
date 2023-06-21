all_q = [];
all_tau_off = [];
all_tau_on = [];
all_bf = [];
all_bs = [];
all_koff = [];
memory_index = [];
gene_number_inf = [];
gene_number_non = [];%收集未推断的基因
for gene_number = 1:1971
    filename = sprintf('E:\\code\\1_test\\results1\\results_mESC\\result_gene_%d',gene_number);
    result_data = load(filename);
    result = getfield(result_data,'result');
    last_iteration_result = result(:,end);%结果中的最后一列

    a = last_iteration_result(1);
    if isempty(a{1,1}) == 1
        gene_number_non(end + 1) = gene_number;
        continue
    end
    dist = cellfun(@(c) c.dist,last_iteration_result);
    [min_dist,min_ind] = min(dist);%min_ind 对应某基因推断出的最优参数索引
    if isfloat(min_dist) == 0
        continue
    end
    
    gene_number_inf(end + 1) = gene_number;
    param_vector = last_iteration_result(min_ind);%某基因推断出的最优参数值
    kon1 = cellfun(@(c) c.kon1,param_vector);
    kon2 = cellfun(@(c) c.kon2,param_vector);
    ron1 = cellfun(@(c) c.ron1,param_vector);
    ron2 = cellfun(@(c) c.ron2,param_vector);
    koff = cellfun(@(c) c.koff,param_vector);
    roff = cellfun(@(c) c.roff,param_vector);
    mu = cellfun(@(c) c.mu,param_vector);
    q = cellfun(@(c) c.q,param_vector);
    all_koff(end + 1) = koff;
    all_q(end + 1) = q;
    tau_off1 = kon1./ron1;% 伽马分布的均值
    tau_off2 = kon2./ron2;
    tau_off = q .* tau_off1 + (1-q) .* tau_off2; %平均off态驻留时间
    all_tau_off(end + 1) = tau_off;
    tau_on = koff./roff; %平均on态驻留时间
    all_tau_on(end + 1) = tau_on;
    bf = 1./(tau_off +tau_on); % 爆发频率
    all_bf(end + 1) = bf;
    bs = mu .* tau_on; % 爆发大小
    all_bs(end + 1) = bs;
    memory_index_1 = q .* kon1 + (1-q) .* kon2;
    memory_index(end + 1) = memory_index_1;
end
save all_result_mESC.mat gene_number_inf all_q all_tau_off all_tau_on all_bf all_bs memory_index all_koff
save mESC_gene_number_non.mat gene_number_non

%subplot(1,2,1)
% scatplot(all_q,log10(all_tau_off),'circles')
% xlabel('q')
% ylabel('log(<\tau_{off}>)')