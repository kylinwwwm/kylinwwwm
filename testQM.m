% Parameters setting
% OFF dwell time distribution f_off(t) = ron^(kon) * t^(kon-1) * e^(-ron * t) / gamma(kon)
% ON dwell time distribution f_on(t) = roff^(koff) * t^(koff-1) * e^(-roff * t) / gamma(koff)
% param_true.kon1 = 3;
% param_true.ron1 = 0.5;
% param_true.kon2 = 6;
% param_true.ron2 = 0.7;
% param_true.koff = 5;
% param_true.roff = 0.5;
% param_true.mu = 20;
% param_true.delta = 1;
% param_true.x0 = [1,0,0];
% param_true.tottime = 2000;
% param_true.q = 0.4;
all_statis_data = [];
for kon1 = 1:0.5:10.5
    param_true.kon1 = kon1;
    param_true.ron1 = 0.5;
    param_true.kon2 = 6;
    param_true.ron2 = 0.7;
    param_true.koff = 1;
    param_true.roff = 0.5;
    param_true.mu = 20;
    param_true.delta = 1;
    param_true.x0 = [1,0,0];
    param_true.tottime = 2000;
    param_true.q = 1;





% Simulation algorithm for GTM.
    [x,t] = simulQM(param_true);
    tq = 1500:0.1:param_true.tottime;
    xq = interp1(t,x(:,3),tq,'previous');%插值  函数关系由t,x(:,3)给出，插值tq对应的xq  x(:,3)表示x中的第三个元素，也就是mRNA分子数                                           
    data = xq;
    time_series = cat(2,t,x);
    time_series2 = cat(1,tq,xq);
    time_series2 = time_series2';
    statis_data = statisData(data);
    true_statisQM = statisQM(param_true,4);
    all_statis_data = [all_statis_data;statis_data];

    save statisData.mat all_statis_data
end
% all_statis_data = load('statisData');
color=[[53 94 103]/255;[246 178 107]/255];
x=1:0.5:10.5;%x轴上的数据，第一个值代表数据开始，第二个值代表间隔，第三个值代表终止
mean = all_statis_data(:,1); %a数据y值
cv2 = all_statis_data(:,2);
fano = all_statis_data(:,3);
sk = all_statis_data(:,4);
kt = all_statis_data(:,5);
bc = all_statis_data(:,6);
subplot(2,3,1)
plot(x,mean,'-','Color',color(1,:),'Marker','o','Markersize',3,'LineWidth',2); %线性，颜色，标记
axis([0,11,0,12])  %确定x轴与y轴框图大小
set(gca,'XTick',[0:5:10]) %x轴范围1-6，间隔1
%set(gca,'XTickLabel',{' ','Mean','CV^2','Fano','Skewness','Kurtosis','Bimodal coefficient'})
set(gca,'YTick',[0:3:12]) %y轴范围0-700，间隔100
legend('Mean');   %右上角标注
xlabel('k_{on}')  %x轴坐标描述
ylabel('Mean') %y轴坐标描述
set(gca,'linewidth',1.5)
set(gca,'FontWeight','bold')
set(gca,'FontName','Times New Rome','FontSize',10);

subplot(2,3,2)
plot(x,cv2,'-','Color',color(2,:),'Marker','o','Markersize',3,'LineWidth',2);
axis([0,11,0,9])
set(gca,'XTick',[0:5:10]) %x轴范围1-6，间隔1
%set(gca,'XTickLabel',{' ','Mean','CV^2','Fano','Skewness','Kurtosis','Bimodal coefficient'})
set(gca,'YTick',[0:3:9]) %y轴范围0-700，间隔100
legend('Noise strength');   %右上角标注
xlabel('k_{on}')  %x轴坐标描述
ylabel('Noise strength') %y轴坐标描述
set(gca,'linewidth',1.5)
set(gca,'FontWeight','bold')
set(gca,'FontName','Times New Rome','FontSize',10);

subplot(2,3,3)
plot(x,fano,'-','Color',color(1,:),'Marker','o','Markersize',3,'LineWidth',2);
axis([0,11,3,15])
set(gca,'XTick',[0:5:10]) %x轴范围1-6，间隔1
%set(gca,'XTickLabel',{' ','Mean','CV^2','Fano','Skewness','Kurtosis','Bimodal coefficient'})
set(gca,'YTick',[3:3:15]) %y轴范围0-700，间隔100
legend('Fano factor');   %右上角标注
xlabel('k_{on}')  %x轴坐标描述
ylabel('Fano factor') %y轴坐标描述
set(gca,'linewidth',1.5)
set(gca,'FontWeight','bold')
set(gca,'FontName','Times New Rome','FontSize',10);

subplot(2,3,4)
plot(x,sk,'-','Color',color(2,:),'Marker','o','Markersize',3,'LineWidth',2);
axis([0,11,0,6])
set(gca,'XTick',[0:5:10]) %x轴范围1-6，间隔1
%set(gca,'XTickLabel',{' ','Mean','CV^2','Fano','Skewness','Kurtosis','Bimodal coefficient'})
set(gca,'YTick',[0:3:6]) %y轴范围0-700，间隔100
legend('Skewness');   %右上角标注
xlabel('k_{on}')  %x轴坐标描述
ylabel('Skewness') %y轴坐标描述
set(gca,'linewidth',1.5)
set(gca,'FontWeight','bold')
set(gca,'FontName','Times New Rome','FontSize',10);

subplot(2,3,5)
plot(x,kt,'-','Color',color(1,:),'Marker','o','Markersize',3,'LineWidth',2);
axis([0,11,0,20])
set(gca,'XTick',[0:5:10]) %x轴范围1-6，间隔1
%set(gca,'XTickLabel',{' ','Mean','CV^2','Fano','Skewness','Kurtosis','Bimodal coefficient'})
set(gca,'YTick',[0:4:20]) %y轴范围0-700，间隔100
legend('Kurtosis');   %右上角标注
xlabel('k_{on}')  %x轴坐标描述
ylabel('Kurtosis') %y轴坐标描述
set(gca,'linewidth',1.5)
set(gca,'FontWeight','bold')
set(gca,'FontName','Times New Rome','FontSize',10);

subplot(2,3,6)
plot(x,bc,'-','Color',color(2,:),'Marker','o','Markersize',3,'LineWidth',2);
axis([0,11,0.4,1])
set(gca,'XTick',[0:5:10]) %x轴范围1-6，间隔1
%set(gca,'XTickLabel',{' ','Mean','CV^2','Fano','Skewness','Kurtosis','Bimodal coefficient'})
set(gca,'YTick',[0.4:0.2:1]) %y轴范围0-700，间隔100
legend('Bimodal coefficient');   %右上角标注
xlabel('k_{on}')  %x轴坐标描述
ylabel('Bimodal coefficient') %y轴坐标描述
set(gca,'linewidth',1.5)
set(gca,'FontWeight','bold')
set(gca,'FontName','Times New Rome','FontSize',10);
% figure
% histogram(data,'normalization','pdf')%横坐标表示数据的取值，纵坐标表示该数据出现的概率
% set(gca,'ytick',0:0.1:1)