%ѹ����֪�ع��㷨����CS_Reconstuction_KtoPercentage.m  
%��ͬ�㷨��ϡ�����ָ��ٷֱ�֮��Ĺ�ϵ
clear all;close all;clc; 

S1 = ['-ks';'-ko';'-kd';'-kv';'-k*'];  
S2 = ['-rs';'-ro';'-rd';'-rv';'-r*'];
S3 = ['-bs';'-bo';'-bd';'-bv';'-b*'];  

%% �������ó�ʼ��  
CNT = 1000;%����ÿ��(K,M,N)���ظ���������  
N = 256;%�ź�x�ĳ���  
Psi = eye(N);%x������ϡ��ģ�����ϡ�����Ϊ��λ��x=Psi*theta  
M_set = [52,100,148,196,244];%����ֵ����

Percentage_OMP = zeros(length(M_set),N);%�洢�ָ��ɹ�����
Percentage_CoSaMP = zeros(length(M_set),N);
Percentage_SP = zeros(length(M_set),N);

%% ��ѭ��������ÿ��(K,M,N)  
tic  
for mm = 1:length(M_set)  
    M = M_set(mm);%���β���ֵ����  
    K_set = 1:5:ceil(M/2);%�ź�x��ϡ���Kû��Ҫȫ��������ÿ��5����һ���Ϳ����� 
    
    PercentageM_OMP = zeros(1,length(K_set));%�洢�˲���ֵM�²�ͬK�Ļָ��ɹ�����
    PercentageM_CoSaMP = zeros(1,length(K_set));
    PercentageM_SP = zeros(1,length(K_set));
    
    for kk = 1:length(K_set)  
       K = K_set(kk);%�����ź�x��ϡ���K
       fprintf('K=%d,M=%d\n',K,M); 
       
       %ͳ�ƻָ��ɹ����źŵ���
       P_OMP = 0; 
       P_CoSaMP = 0; 
       P_SP = 0; 
       
       for cnt = 1:CNT %ÿ���۲�ֵ����������CNT��  
            Index_K = randperm(N);  
            x = zeros(N,1);  
            x(Index_K(1:K)) = 5*randn(K,1);%xΪKϡ��ģ���λ���������
            
            %% OMP�����źŻָ����ָ����������
            Phi = randn(M,N);%��������Ϊ��˹����  
            A = Phi * Psi;%���о���  
            y = Phi * x;%�õ��۲�����y  
            theta = CS_OMP(y,A,K);%�ָ��ع��ź�theta  
            x_r = Psi * theta;% x=Psi * theta  
            if norm(x_r-x)<1e-6%����в�С��1e-6����Ϊ�ָ��ɹ�  
                P_OMP = P_OMP + 1;  
            end 
            
            %% CoSaMP�����źŻָ����ָ����������
            Phi = randn(M,N)/sqrt(M);%��������Ϊ��˹����  
            A = Phi * Psi;%���о���  
            y = Phi * x;%�õ��۲�����y  
            theta = CS_CoSaMP(y,A,K);%�ָ��ع��ź�theta  
            x_r = Psi * theta;% x=Psi * theta  
            if norm(x_r-x)<1e-6%����в�С��1e-6����Ϊ�ָ��ɹ�  
                P_CoSaMP = P_CoSaMP + 1;  
            end 
            
            %% SP�����źŻָ����ָ����������
            Phi = randn(M,N)/sqrt(M);%��������Ϊ��˹����  
            A = Phi * Psi;%���о���  
            y = Phi * x;%�õ��۲�����y  
            theta = CS_SP(y,A,K);%�ָ��ع��ź�theta  
            x_r = Psi * theta;% x=Psi * theta  
            if norm(x_r-x)<1e-6%����в�С��1e-6����Ϊ�ָ��ɹ�  
                P_SP = P_SP + 1;  
            end 
            
            
       end 
       
       PercentageM_OMP(kk) = P_OMP/CNT*100;%����ָ����� 
       PercentageM_CoSaMP(kk) = P_CoSaMP/CNT*100;
       PercentageM_SP(kk) = P_SP/CNT*100;  
    end  
    Percentage_OMP(mm,1:length(K_set)) = PercentageM_OMP;
    Percentage_CoSaMP(mm,1:length(K_set)) = PercentageM_CoSaMP;
    Percentage_SP(mm,1:length(K_set)) = PercentageM_SP;
    
end  
toc  
save KtoPercentage1000test %����һ�β����ף��ѱ���ȫ���洢����  
%% ��ͼ  

figure(1);  
for mm = 1:length(M_set)  
    M = M_set(mm);  
    K_set = 1:5:ceil(M/2);  
    L_Kset = length(K_set);  
    plot(K_set,Percentage_OMP(mm,1:L_Kset),S(mm,:));%���x�Ļָ��ź�  
    hold on;  
end  
hold off;  
xlim([0 125]);  
legend('M=52','M=100','M=148','M=196','M=244');  
xlabel('Sparsity level(K)');  
ylabel('Percentage recovered');  
title('Percentage of input signals recovered correctly(OMP)');  


figure(2);  
for mm = 1:length(M_set)  
    M = M_set(mm);  
    K_set = 1:5:ceil(M/2);  
    L_Kset = length(K_set);  
    plot(K_set,Percentage_CoSaMP(mm,1:L_Kset),S(mm,:));%���x�Ļָ��ź�  
    hold on;  
end  
hold off;  
xlim([0 125]);  
legend('M=52','M=100','M=148','M=196','M=244');  
xlabel('Sparsity level(K)');  
ylabel('Percentage recovered');  
title('Percentage of input signals recovered correctly(CoSaMP)');



figure(3);  
for mm = 1:length(M_set)  
    M = M_set(mm);  
    K_set = 1:5:ceil(M/2);  
    L_Kset = length(K_set);  
    plot(K_set,Percentage_SP(mm,1:L_Kset),S(mm,:));%���x�Ļָ��ź�  
    hold on;  
end  
hold off;  
xlim([0 125]);  
legend('M=52','M=100','M=148','M=196','M=244');  
xlabel('Sparsity level(K)');  
ylabel('Percentage recovered');  
title('Percentage of input signals recovered correctly(SP)');