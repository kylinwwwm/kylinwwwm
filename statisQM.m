function statis = statisQM(param,k)
%% This code computes some statistics of the GTM.计算统计量
%% 参数设置

koff = param.koff; %状态切换反应步数参数
roff = param.roff; %状态切换速率参数
kon1 = param.kon1;
ron1 = param.ron1;
kon2 = param.kon2;
ron2 = param.ron2;
mu = param.mu; %转录速率
delta = param.delta; %降解速率
q = param.q; %选择第一条通路的概率


Laplace_s = 0:k;
Lfon1 = (ron1./(Laplace_s + ron1)).^kon1; % Laplace fon1(x)索引比取值大1
Lfon2 = (ron2./(Laplace_s + ron2)).^kon2; % Laplace fon2(x)
Lfon = q .* Lfon1 + (1 - q) .* Lfon2;
Lfoff = (roff./(Laplace_s + roff)).^koff;% Laplace foff(x)
mean_tau_off1 = kon1/ron1; % 伽马分布的均值
mean_tau_off2 = kon2/ron2;
mean_tau_off = q * mean_tau_off1 + (1-q) * mean_tau_off2;
mean_tau_on = koff/roff;
c = 1/(mean_tau_off + mean_tau_on);%BF
LFon1 = (1 - Lfon1(2:end))./Laplace_s(2:end); % Laplace Fon1(x)
LFon2 = (1 - Lfon2(2:end))./Laplace_s(2:end);% 索引等于取值
LFon = q .* LFon1 + (1 - q) .* LFon2;
LFoff = (1 - Lfoff(2:end))./Laplace_s(2:end); % Laplace Foff(x)
vk = (1 - Lfon(1:end))./(1 - Lfon(1:end) .* Lfoff(1:end));% 取值比索引小一
b1 = c * mu * mean_tau_on  ;
b2 = (mu^2 * mean_tau_on / (2 * (mean_tau_on + mean_tau_off))) - (vk(2) * LFoff(1) * mu^2 / (2 * (mean_tau_on + mean_tau_off)));%vk的索引比取值大1：v1=v(2)
b3 = (c * mu^3 * mean_tau_on / 6) - (c * vk(2) * mu^3 * LFoff(1) / 6) - (c * vk(3) * mu^3 * LFoff(2) / 6) + (c * vk(2) * vk(3) * mu ^3 * (Lfoff(2) - Lfoff(3)) / 6);
b4 = (c * mu^4 * mean_tau_on / 24) - (c * vk(2) * vk(4) * mu^4 * (1 - Lfon(3) * Lfoff(4)) * (Lfoff(2) - Lfoff(3))) / (12 * (1 - Lfoff(3) * Lfon(3))) - ...
    (c * mu^4 * (vk(2) * LFoff(1) + vk(3) * LFoff(2) + vk(4) * LFoff(3)) / 24) + (c * mu^4 * (vk(2) * vk(3) * (Lfoff(2) - Lfoff(3)) + vk(3) * vk(4) * (Lfoff(3) - Lfoff(4)) + ...
    vk(2) * vk(4) * (Lfoff(4) - Lfoff(2))) / 24);
% LFon = [mean_tau_off (1-Lfon(2:end))./Laplace_s(2:end)]; % Laplace Fon(x)
% LFoff = [mean_tau_on (1-Lfoff(2:end))./Laplace_s(2:end)]; % Laplace Foff(x)
bf = 1/(mean_tau_off + mean_tau_on); % 爆发频率
bs = mu * mean_tau_on; % 爆发大小
          
m = b1;
v = 2*b2+b1-b1^2;
cv2 = v/(m^2);
fano = v/m;
sk = (6*b3+6*b2+b1-3*b1*(2*b2+b1)+2*b1^3)/(2*b2+b1-b1^2)^1.5 + 1;%不知道为啥加了个1
kt = ((24*b4+36*b3+14*b2+b1-4*b1*(6*b3+6*b2+b1)+...
   6*b1^2*(2*b2+b1)-3*b1^4)/(2*b2+b1-b1^2)^2 )+ 1;%不知道为啥加了个1
% kt = (b1+14*b2+36*b3+24*b4-7*b1^2-12*b2^2-36*b1*b2-24*b1*b3+12*b1^3+24*b1^2*b2-6*b1^4)/(2*b2+b1-b1^2)^2 ;
bc = ((sk - 1)^2 + 1)/kt;
%statis = [m,cv2,fano,sk,kt,bc];   
statis = [m,cv2,fano,sk];  
end
