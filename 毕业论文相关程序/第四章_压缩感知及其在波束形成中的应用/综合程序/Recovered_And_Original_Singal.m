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

            %%   SP方法的恢复信号与原始信号对比图以及恢复误差计算
Error_SP=zeros(C,1);
%% 恢复重构信号x  
tic  
for count=1:C
        theta = CS_SP( y,A,K );  
        x_r_SP = Psi * theta;% x=Psi * theta
        Error_SP(count)=norm(x_r_SP-x);
end 
toc  
e_SP=mean(Error_SP)
%% 绘图  
figure(1);  
plot(x_r_SP,'ko-');%绘出x的恢复信号  
hold on;  
plot(x,'r*--');%绘出原信号x  
hold off;  
legend('Recovery','Original');
xlabel('Length Of Singal');  
ylabel('Discrete Singal Amplitude');  
title('Recovered And Original Singal(SP) ');  


%%   CoSaMP方法的恢复信号与原始信号对比图以及恢复误差计算
Error_CoSaMP=zeros(C,1);
%% 恢复重构信号x  
tic  
for count=1:C
        theta = CS_CoSaMP( y,A,K );  
        x_r_CoSaMP = Psi * theta;% x=Psi * theta
        Error_CoSaMP(count)=norm(x_r_CoSaMP-x);
end 
toc  
e_CoSaMP=mean(Error_CoSaMP)
%% 绘图  
figure(2);  
plot(x_r_CoSaMP,'ko-');%绘出x的恢复信号  
hold on;  
plot(x,'r*--');%绘出原信号x  
hold off;  
legend('Recovery','Original');
xlabel('Length Of Singal');  
ylabel('Discrete Singal Amplitude');  
title('Recovered And Original Singal(CoSaMP) '); 



%%   OMP方法的恢复信号与原始信号对比图以及恢复误差计算
Error_OMP=zeros(C,1);
%% 恢复重构信号x  
tic  
for count=1:C
        theta = CS_OMP( y,A,K );  
        x_r_OMP = Psi * theta;% x=Psi * theta
        Error_OMP(count)=norm(x_r_OMP-x);
end 
toc  
e_OMP=mean(Error_OMP)
%% 绘图  
figure(3);  
plot(x_r_OMP,'ko-');%绘出x的恢复信号  
hold on;  
plot(x,'r*--');%绘出原信号x  
hold off;  
legend('Recovery','Original');
xlabel('Length Of Singal');  
ylabel('Discrete Singal Amplitude');  
title('Recovered And Original Singal(OMP) ');  




