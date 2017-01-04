% 阵列方向图（使用MVDR方法计算加权向量）
% 阵元数为16，间距为半波长，信号频率为1.3G，指向角为0度，干扰三个在－40，20，45
% 考虑波束指向方向无信号情况
% 2007.5.22无信号

clear all
clc
f0 = 1.3*10^9;                    % 信号频率
fs = 3*10^9;                   % 采样频率
snr = [40;40];                  % 干扰的信噪比
w0 = 20/180*pi;                  % 目标指向角
M = 9;                         % 阵元数为M
K = 2;                          % 干扰信号数目
seta =[20/180*pi,45/180*pi];   % 干扰方向          
N = 100;                        % 采样点数
c = 3*10^8;                     % 光速
d = 0.5*c/f0;                   % 阵元间距
x = zeros(M,1);                 % 数据矢量
R = zeros(M,M);                 % 接收数据协方差矩阵
tic
% 信号源
nn = 1:N;%下面设计的是干扰信号,nn为列N的行向量
for k=1:K
    for n=1:N
        fai(1,n) = rand;
    end
    s(k,:) = exp(i*2*pi*(f0*nn/fs+fai*k));%构造信号
end

% 构造方向矢量
theta = -90*pi/180:pi/180:90*pi/180;
WW = length(theta);
P = zeros(1,WW);
for m=1:M
        a(m,:) = exp(-i*pi*(m-1)*sin(theta)); % 方向向量，用于方向搜索
end
 
% 接收数据协方差矩阵R
X=[];
for n=1:N                  
    for m=1:M
        for k=1:K
             A(m,k) = exp(-i*pi*(m-1)*sin(seta(k)));
        end
    end
    noise(:,n)=randn(M,1);
    y = A*(10.^(snr/20).*s(:,n))+noise(:,n);
    R = R+y*y';
    
    X=[X,y];
end
R = R/N;
Ri = inv(R);

% 构造加权向量
for k=1:M
    a0(k,1) =exp(-i*pi*(k-1)*sin(w0));  % 指向向量 
end 
W = Ri*a0/(a0'*Ri*a0);
f = abs(W'*a).^2;
figure;
f = f./max(f);
%load QRD_LS_f
plot(theta*180/pi,10*log10(f),':r');
%theta*180/pi,10*log10(f_LS),'b-');%
% title('power pattern of linear phased array');
xlabel('\theta/deg'); 
ylabel('阵列增益（dB）');
grid on
toc
