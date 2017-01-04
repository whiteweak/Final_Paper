%MVDR谱估计
clear,clc,close all;
N=256;
noise=(randn(1,N)+1j*randn(1,N))/sqrt(2);  %高斯白噪声
f1=0.1;
f2=0.25;
f3=0.27;
SNR1=30;
SNR2=30;
SNR3=27;
A1=sqrt(2)*10^(SNR1/20);
A2=sqrt(2)*10^(SNR2/20);
A3=sqrt(2)*10^(SNR3/20);
n=0:N-1;
sig1=A1*exp(1j*2*pi*f1*n+1j*2*pi*rand(1,1));
sig2=A2*exp(1j*2*pi*f2*n+1j*2*pi*rand(1,1));
sig3=A2*exp(1j*2*pi*f3*n+1j*2*pi*rand(1,1));
u=sig1+sig2+sig3+noise;

M=16;      %自相关矩阵的阶数,阶数越高，谱分辨率越高
for k=1:N-M
    xs(:,k)=u(k+M-1:-1:k)';   %构造样本矩阵
end
R=xs*xs'/(N-M);

%计算MVDR谱
NF=2048;     %谱峰搜索点数
f=linspace(-0.5,0.5,NF);
for n=1:NF
    a=exp(1j*2*pi*f(n)*(0:M-1)');   %定义向量a(w),为什么没负号？？？
    Pmvdr(n)=1/(a'*inv(R)*a);
end
Pmvdr=10*log10(abs(Pmvdr)/max(abs(Pmvdr))); %归一化MVDR谱
plot(f,Pmvdr),grid
set(gca,'xtick',[f1,f2,f3]);
title('MVDR谱估计'),ylabel('归一化MVDR谱/dB'),xlabel('w/2\pi')