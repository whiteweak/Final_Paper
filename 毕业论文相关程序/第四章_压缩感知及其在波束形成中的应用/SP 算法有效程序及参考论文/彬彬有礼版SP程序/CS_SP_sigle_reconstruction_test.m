%压缩感知重构算法测试  
clear all;close all;clc;  
M = 64;%观测值个数  
N = 256;%信号x的长度  
K = 12;%信号x的稀疏度  
Index_K = randperm(N);  
x = zeros(N,1);  
x(Index_K(1:K)) = 5*randn(K,1);%x为K稀疏的，且位置是随机的  
Psi = eye(N);%x本身是稀疏的，定义稀疏矩阵为单位阵x=Psi*theta  
Phi = randn(M,N);%测量矩阵为高斯矩阵  
A = Phi * Psi;%传感矩阵  
y = Phi * x;%得到观测向量y  
C=1000;
Errors=zeros(C,1);
%% 恢复重构信号x  
tic  
for count=1:C
        theta = CS_SP( y,A,K );  
        x_r = Psi * theta;% x=Psi * theta
        Errors(count)=norm(x_r-x);
end 
toc  
e=mean(Errors)
%% 绘图  
figure;  
plot(x_r,'ko-');%绘出x的恢复信号  
hold on;  
plot(x,'r*--');%绘出原信号x  
hold off;  

legend('Recovery','Original');
xlabel('Length Of Singal');  
ylabel('Discrete Singal Amplitude');  
title('Recovered And Original Singal ');  
% error=x_r-x;
% fprintf('\n恢复残差：');  
% norm(x_r-x)%恢复残差  





