function  [] = figureResult(data,result,param_true)

kon1 = cellfun(@(c) c.kon1,result);
kon2 = cellfun(@(c) c.kon2,result);%result中最后一列参数的kon  分别取出最后一次迭代产生的1000个元胞数组中的参数
ron1 = cellfun(@(c) c.ron1,result);
ron2 = cellfun(@(c) c.ron2,result);
koff = cellfun(@(c) c.koff,result);
roff = cellfun(@(c) c.roff,result);
mu = cellfun(@(c) c.mu,result);
q = cellfun(@(c) c.q,result);
dist = cellfun(@(c) c.dist,result);

%推断出的结果
tau_off1 = kon1./ron1;% 伽马分布的均值
tau_off2 = kon2./ron2;
tau_off = q .* tau_off1 + (1-q) .* tau_off2;
tau_on = koff./roff;
bf = 1./(tau_off +tau_on); % 爆发频率
bs = mu .* tau_on; % 爆发大小

%预先设置好的真实参数
param_true.tau_off1 = param_true.kon1./param_true.ron1;
param_true.tau_off2 = param_true.kon2./param_true.ron2;
param_true.tau_off = param_true.q * param_true.tau_off1 + (1-param_true.q) * param_true.tau_off2;
param_true.tau_on = param_true.koff./param_true.roff;
param_true.bf = 1./(param_true.tau_on + param_true.tau_off);
param_true.bs = param_true.mu .* param_true.tau_on;


[f, tau] = ksdensity([log10(tau_off), log10(tau_on)]);%[f,xi] = ksdensity(x)计算样本x的密度估计，返回在xi点的密度f，此时我们使用plot(xi,f)就可以绘制出谱密度曲线。该函数，首先统计样本x在各个区间的概率(与hist有些相似)，再自动选择xi，计算对应的xi点的谱密度

tau_offCenter = tau(find(f == max(f)), 1);%find返回max(f)的第一个索引，
tau_onCenter = tau(find(f == max(f)), 2);
%Idx = knnsearch(log10([tau_off,tau_on]),[tau_offCenter,tau_onCenter]);%knnsearch在第一个参数矩阵中搜索离第二个参数矩阵中每一行最近的行的索引（欧氏距离）

[~,min_ind] = min(dist);%选择与真实数据差值最小的一组参数作为参数估计结果绘图
param_est = result{min_ind,1};
param_est.x0 = [1,0,0];
param_est.tottime = 2000;
[x_est,t_est] = simulQM(param_est);
tq_est = 1000:0.1:param_est.tottime;
xq_est = interp1(t_est,x_est(:,3),tq_est,'previous');
data_est = xq_est;
pts = linspace(1,max(data)+10,100);
[f_true,xi_true] = ksdensity(data + 1,pts,'Support','positive','Bandwidth',0.3); 

%%%%%%
Color = {'#EE2201'};
FaceColor = {'#00837E','#4DBBD4'};

f1 = figure;
histogram(data_est,'normalization','pdf','BinEdges',0:max(data_est),'EdgeColor','none','FaceColor',FaceColor{1})%画出由估计的参数产生的模拟数据
hold on
plot(xi_true-0.5,f_true,'r','LineWidth',1,'Color',Color{1},'LineWidth',1.5)%应该是画出提前设置好的真实参数产生数据的线图
xlim([0 max(data_est)+5])
xlabel('mRNA')
ylabel('PDF')
set(gca,'TickLength',[0.02,0.025]);
set(gca,'FontName','Arial','FontSize',6);
set(f1,'position',[400 400 200 120]);

f2 = figure;
h = histogram(mu,'normalization','pdf','EdgeColor','none','FaceColor',FaceColor{2});%画mRNA合成速率mu的分布图
[~, max_ind] = max(h.Values);%估计出的mu的众数
muCenter = mean(h.BinEdges(max_ind:max_ind + 1));
hold on;
plot([param_true.mu param_true.mu], get(gca, 'YLim'), '-r', 'LineWidth', 2) 
plot([muCenter muCenter], get(gca, 'YLim'), '--k', 'LineWidth', 2) 
xlabel('\mu')
ylabel('PDF')
set(gca,'TickLength',[0.02,0.025]);
set(gca,'FontName','Arial','FontSize',6);
set(f2,'position',[400 400 120 120]);
set(gca,'linewidth',1.5)
set(gca,'FontWeight','bold')
set(gca,'FontSize',10);

f3 = figure;
subplot(1,2,1)
% log10(tau_off)
% log10(tau_on)
scatplot(log10(tau_off),log10(tau_on),'circles')
hold on
plot([log10(param_true.tau_off) log10(param_true.tau_off)], get(gca, 'YLim'), '-r', 'LineWidth', 1) %红色实线对应真是参数
plot([tau_offCenter tau_offCenter], get(gca, 'YLim'), '--k', 'LineWidth', 1) %黑色虚线代表推断出参数
plot(get(gca, 'XLim'), [log10(param_true.tau_on) log10(param_true.tau_on)], '-r', 'LineWidth', 1) %红色实线
plot(get(gca, 'XLim'), [tau_onCenter tau_onCenter], '--k', 'LineWidth', 1) %黑色虚线
xlim(get(gca, 'XLim'))
ylim(get(gca, 'YLim'))
xlabel('log(<\tau_{off}>)')
ylabel('log(<\tau_{on}>)')
set(gca,'TickLength',[0.02,0.025]);
set(gca,'FontName','Arial','FontSize',6);
axis square
axis on
box on

subplot(1,2,2)
scatplot(log10(bf),log10(bs),'circles')
hold on
plot([log10(param_true.bf) log10(param_true.bf)], get(gca, 'YLim'), '-r', 'LineWidth', 1) 
[f, log_burst] = ksdensity([log10(bf), log10(bs)]);
log_bfCenter = log_burst(find(f == max(f)), 1);
log_bsCenter = log_burst(find(f == max(f)), 2);
plot([log_bfCenter log_bfCenter], get(gca, 'YLim'), '--k', 'LineWidth', 1) 
plot(get(gca, 'XLim'), [log10(param_true.bs) log10(param_true.bs)], '-r', 'LineWidth', 1) 
plot(get(gca, 'XLim'), [log_bsCenter log_bsCenter], '--k', 'LineWidth', 1)

xlim(get(gca, 'XLim'))
ylim(get(gca, 'YLim'))
xlabel('log(bf)')
ylabel('log(bs)')
set(gca,'TickLength',[0.02,0.025]);
set(gca,'FontName','Arial','FontSize',6);
axis square
box on
set(f3,'position',[400 400 400 200]);
end


