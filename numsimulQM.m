function [bs,bf] = numsimulQM(param)
%% 参数设置
koff = param.koff; %状态切换反应步数参数
roff = param.roff; %状态切换速率参数
kon1 = param.kon1;
ron1 = param.ron1;
kon2 = param.kon2;
ron2 = param.ron2;
mu = param.mu; %转录速率
delta = param.delta; %降解速率
x = param.x0; %初始模拟点[OFF,ON,mRNA]=[1,0,0],代表初始状态处于OFF态，mRNA分子数为0
tottime = param.tottime;



%% 反应矩阵
r_mu = [-1 1 0;
    1 -1 0;
    0 0 1;
    0 0 -1]; %前两行代表状态切换，后两行代表mRNA的合成和降解

t = 0;
n = 0; %设置爆发次数计数器
j = 0; %设置mRNA总量计数器
%% 模拟反应
while (t(end) < tottime)
    if x(end,1) == 1 %OFF态
        q = rand;
        r = rand;
        if r < q
            tau_off = gamrnd(kon1,1/ron1);
        else
            tau_off = gamrnd(kon2,1/ron2); %OFF态驻留时间
        end
        t_temp = t(end) + tau_off;
        while t(end) < t_temp
            a_0 = delta * x(end,3);
            r1 = rand;
            tau_1 = (1/a_0) * log(1/r1);
            t = [t;t(end) + tau_1];
            x = [x;x(end,:) + r_mu(4,:)];
        end
        t = [t(1:end-1);t_temp];
        x = [x(1:end-1,:);x(end-1,:) + r_mu(1,:)]; %从OFF切换到ON
        
    elseif x(end,1) == 0 % ON状态
        tau_on = gamrnd(koff,1/roff);
        t_temp = t(end) + tau_on; %ON态驻留时间
        i = 0;%设置一次爆发产生的mRNA数量计数器
        while t(end) < t_temp
            a_mu(1) = mu * x(end,2);
            a_mu(2) = delta * x(end,3);
            a_0 = sum(a_mu);
            r1 = rand;
            tau_2 = (1/a_0) * log(1/r1);
            
            r2 = rand;
            for iter = 1:2
                if (sum(a_mu(1:iter)) >= r2 * a_0)
                    next_mu = iter;
                    if next_mu == 1
                        i = i + 1;
                    end
                    break;
                end
            end
            t = [t;t(end) + tau_2];
            x = [x;x(end,:) + r_mu(next_mu + 2,:)];
        end
        j = i + j;
        t = [t(1:end-1);t_temp];
        x = [x(1:end-1,:);x(end-1,:) + r_mu(2,:)];%状态切换
        n = n + 1;
    end
end

bs = j / n; 
bf = n / 2000;