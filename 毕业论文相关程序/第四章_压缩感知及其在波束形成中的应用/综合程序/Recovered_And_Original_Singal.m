%ѹ����֪�ع��㷨����  
clear all;close all;clc;  
M = 64;%�۲�ֵ����  
N = 256;%�ź�x�ĳ���  
K = 12;%�ź�x��ϡ���  
Index_K = randperm(N);  
x = zeros(N,1);  
x(Index_K(1:K)) = 5*randn(K,1);%xΪKϡ��ģ���λ���������  
Psi = eye(N);%x������ϡ��ģ�����ϡ�����Ϊ��λ��x=Psi*theta  
Phi = randn(M,N);%��������Ϊ��˹����  
A = Phi * Psi;%���о���  
y = Phi * x;%�õ��۲�����y  
C=1000;

            %%   SP�����Ļָ��ź���ԭʼ�źŶԱ�ͼ�Լ��ָ�������
Error_SP=zeros(C,1);
%% �ָ��ع��ź�x  
tic  
for count=1:C
        theta = CS_SP( y,A,K );  
        x_r_SP = Psi * theta;% x=Psi * theta
        Error_SP(count)=norm(x_r_SP-x);
end 
toc  
e_SP=mean(Error_SP)
%% ��ͼ  
figure(1);  
plot(x_r_SP,'ko-');%���x�Ļָ��ź�  
hold on;  
plot(x,'r*--');%���ԭ�ź�x  
hold off;  
legend('Recovery','Original');
xlabel('Length Of Singal');  
ylabel('Discrete Singal Amplitude');  
title('Recovered And Original Singal(SP) ');  


%%   CoSaMP�����Ļָ��ź���ԭʼ�źŶԱ�ͼ�Լ��ָ�������
Error_CoSaMP=zeros(C,1);
%% �ָ��ع��ź�x  
tic  
for count=1:C
        theta = CS_CoSaMP( y,A,K );  
        x_r_CoSaMP = Psi * theta;% x=Psi * theta
        Error_CoSaMP(count)=norm(x_r_CoSaMP-x);
end 
toc  
e_CoSaMP=mean(Error_CoSaMP)
%% ��ͼ  
figure(2);  
plot(x_r_CoSaMP,'ko-');%���x�Ļָ��ź�  
hold on;  
plot(x,'r*--');%���ԭ�ź�x  
hold off;  
legend('Recovery','Original');
xlabel('Length Of Singal');  
ylabel('Discrete Singal Amplitude');  
title('Recovered And Original Singal(CoSaMP) '); 



%%   OMP�����Ļָ��ź���ԭʼ�źŶԱ�ͼ�Լ��ָ�������
Error_OMP=zeros(C,1);
%% �ָ��ع��ź�x  
tic  
for count=1:C
        theta = CS_OMP( y,A,K );  
        x_r_OMP = Psi * theta;% x=Psi * theta
        Error_OMP(count)=norm(x_r_OMP-x);
end 
toc  
e_OMP=mean(Error_OMP)
%% ��ͼ  
figure(3);  
plot(x_r_OMP,'ko-');%���x�Ļָ��ź�  
hold on;  
plot(x,'r*--');%���ԭ�ź�x  
hold off;  
legend('Recovery','Original');
xlabel('Length Of Singal');  
ylabel('Discrete Singal Amplitude');  
title('Recovered And Original Singal(OMP) ');  




