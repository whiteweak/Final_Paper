%ѹ����֪�ع��㷨����CS_Reconstuction_KtoPercentage.m  
%   ���Ʋο������е�Fig.2  
%   �ο����ף�Joel A. Tropp and Anna C. Gilbert   
%   Signal Recovery From Random Measurements Via Orthogonal Matching  
%   Pursuit��IEEE TRANSACTIONS ON INFORMATION THEORY, VOL. 53, NO. 12,  
%   DECEMBER 2007.  
%   Elapsed time is 1448.966882 seconds.(@20150418night)  
clear all;close all;clc;  
%% �������ó�ʼ��  
CNT = 1000;%����ÿ��(K,M,N)���ظ���������  
N = 256;%�ź�x�ĳ���  
Psi = eye(N);%x������ϡ��ģ�����ϡ�����Ϊ��λ��x=Psi*theta  
M_set = [52,100,148,196,244];%����ֵ����  
Percentage = zeros(length(M_set),N);%�洢�ָ��ɹ�����  
%% ��ѭ��������ÿ��(K,M,N)  
tic  
for mm = 1:length(M_set)  
    M = M_set(mm);%���β���ֵ����  
    K_set = 1:5:ceil(M/2);%�ź�x��ϡ���Kû��Ҫȫ��������ÿ��5����һ���Ϳ�����  
    PercentageM = zeros(1,length(K_set));%�洢�˲���ֵM�²�ͬK�Ļָ��ɹ�����  
    for kk = 1:length(K_set)  
       K = K_set(kk);%�����ź�x��ϡ���K  
       fprintf('K=%d,M=%d\n',K,M);  
       P = 0;  
       for cnt = 1:CNT %ÿ���۲�ֵ����������CNT��  
            Index_K = randperm(N);  
            x = zeros(N,1);  
            x(Index_K(1:K)) = 5*randn(K,1);%xΪKϡ��ģ���λ���������                  
            Phi = randn(M,N);%��������Ϊ��˹����  
            A = Phi * Psi;%���о���  
            y = Phi * x;%�õ��۲�����y  
            theta = CS_OMP(y,A,K);%�ָ��ع��ź�theta  
            x_r = Psi * theta;% x=Psi * theta  
            if norm(x_r-x)<1e-6%����в�С��1e-6����Ϊ�ָ��ɹ�  
                P = P + 1;  
            end  
       end  
       PercentageM(kk) = P/CNT*100;%����ָ�����  
    end  
    Percentage(mm,1:length(K_set)) = PercentageM;  
end  
toc  
save KtoPercentage1000test %����һ�β����ף��ѱ���ȫ���洢����  
%% ��ͼ  
S = ['-ks';'-ko';'-kd';'-kv';'-k*'];  
figure;  
for mm = 1:length(M_set)  
    M = M_set(mm);  
    K_set = 1:5:ceil(M/2);  
    L_Kset = length(K_set);  
    plot(K_set,Percentage(mm,1:L_Kset),S(mm,:));%���x�Ļָ��ź�  
    hold on;  
end  
hold off;  
xlim([0 125]);  
legend('M=52','M=100','M=148','M=196','M=244');  
xlabel('Sparsity level(K)');  
ylabel('Percentage recovered');  
title('Percentage of input signals recovered correctly(N=256)(Gaussian)');  