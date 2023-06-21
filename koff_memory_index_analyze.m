result = load('all_result');
MEF_feature = load('MEF_feature_result');
memory_index = result.all_koff;
bimodality_coefficient = MEF_feature.data_bc;
[Y,I] = sort(memory_index);
x1 = [];
x2 = [];
x3 = [];
x4 = [];
x5 = [];
for i = I
    if memory_index(i) <= 2
        x1(end + 1) = bimodality_coefficient(i);
    end
end

for i = I
    if memory_index(i) <= 4 && memory_index(i) > 2
        x2(end + 1) = bimodality_coefficient(i);
    end
end


for i = I
    if memory_index(i) <= 6 && memory_index(i) > 4
        x3(end + 1) = bimodality_coefficient(i);
    end
end

for i = I
    if memory_index(i) <= 8 && memory_index(i) > 6
        x4(end + 1) = bimodality_coefficient(i);
    end
end

for i = I
    if memory_index(i) > 8
        x5(end + 1) = bimodality_coefficient(i);
    end
end


figure;%memory index---bc
bc = [x1,x2,x3,x4,x5];
mean_bc = [mean(x1),mean(x2),mean(x3),mean(x4),mean(x5)];
bc_tmp = [ones(size(x1)),2*ones(size(x2)),3*ones(size(x3)),4*ones(size(x4)),5*ones(size(x5))];
h = boxplot(bc,bc_tmp,'Symbol','o','Color','k','Whisker',3,'Labels',["[0,2]","(2,4]","(4,6]","(6,8]","(8,inf]"]); %whisker可以决定异常值显示的多少
ylabel("Bimodality coefficient","Fontname","Times New Roman");
xlabel("The Range of koff Memory Index","Fontname","Times New Roman");
% mycolor = ['#00837E','#4DBBD4','#00837E','#4DBBD4','#00837E'];
% mycolor = [255,240,245;123,104,238;152,251,152;211,160,152;255,222,173] ./255;
set(h,'LineWidth',1.5)
mycolor = [251,202,197;185,230,249;251,202,197;185,230,249;251,202,197] ./255;%RGB颜色表
boxobj = findobj(gca,'Tag','Box');
% outliersobj = findobj(gca,'Tag','Outliers');
% lineobj = findobj('gca','Tag','Median');

for i = 1:5
    patch(get(boxobj(i),'XData'),get(boxobj(i),'YData'),mycolor(i,:),'FaceAlpha',0.5,'LineWidth',1.1);
end
hold on;
set(gca,'FontSize',12,'Fontname', 'Times News Roman');

hold on;
plot(mean_bc(1),'bo','Marker','v','Color','k','LineWidth',1.5)
hold on;
plot(2*ones(size(x2)),mean_bc(2),'bo','Marker','v','Color','k','LineWidth',1.5)
hold on;
plot(3*ones(size(x3)),mean_bc(3),'bo','Marker','v','Color','k','LineWidth',1.5)
hold on;
plot(4*ones(size(x4)),mean_bc(4),'bo','Marker','v','Color','k','LineWidth',1.5)
hold on;
plot(5*ones(size(x5)),mean_bc(5),'bo','Marker','v','Color','k','LineWidth',1.5)
hold on;