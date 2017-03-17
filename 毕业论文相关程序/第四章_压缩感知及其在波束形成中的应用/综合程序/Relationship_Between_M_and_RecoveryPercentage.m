

%ѹ����֪�ع��㷨����CS_Reconstuction_MtoPercentage.m
%��ͬ�㷨���źŲ�������ָ��ٷֱ�֮��Ĺ�ϵ

clear all;close all;clc;  

S1 = ['-ks';'-ko';'-kd';'-kv';'-k*'];  
S2 = ['-rs';'-ro';'-rd';'-rv';'-r*'];
S3 = ['-bs';'-bo';'-bd';'-bv';'-b*']; 
%% �������ó�ʼ��  
CNT = 1000;%����ÿ��(K,M,N)���ظ���������  
N = 256;%�ź�x�ĳ���  
Psi = eye(N);%x������ϡ��ģ�����ϡ�����Ϊ��λ��x=Psi*theta  
K_set = [4,12,20,28,36];%�ź�x��ϡ��ȼ���  
Percentage_OMP = zeros(length(K_set),N);%�洢�ָ��ɹ�����
Percentage_CoSaMP = zeros(length(K_set),N);
Percentage_SP = zeros(length(K_set),N);
%% ��ѭ��������ÿ��(K,M,N)  
tic  
for kk = 1:length(K_set)  
    K = K_set(kk);%����ϡ���  
    M_set = K:5:N;%Mû��Ҫȫ��������ÿ��5����һ���Ϳ����� 
    
    PercentageK_OMP = zeros(1,length(M_set));%�洢��ϡ���K�²�ͬM�Ļָ��ɹ�����
    PercentageK_CoSaMP = zeros(1,length(M_set));
    PercentageK_SP = zeros(1,length(M_set));
    
    for mm = 1:length(M_set)  
       M = M_set(mm);%���ι۲�ֵ����  
       fprintf('K=%d,M=%d\n',K,M);
       
       %ͳ�ƻָ��ɹ����źŵ���
       P_OMP=0;
       P_CoSaMP = 0;
       P_SP=0;
       
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
       
       PercentageK_OMP(mm) = P_OMP/CNT*100;%����ָ����� 
       PercentageK_CoSaMP(mm) = P_CoSaMP/CNT*100;
       PercentageK_SP(mm) = P_SP/CNT*100;
       
    end  
    Percentage_OMP(kk,1:length(M_set)) = PercentageK_OMP;
    Percentage_CoSaMP(kk,1:length(M_set)) = PercentageK_CoSaMP;
    Percentage_SP(kk,1:length(M_set)) = PercentageK_SP;
end  
toc  
save MtoPercentage1000 %����һ�β����ף��ѱ���ȫ���洢����  
%% ��ͼ  
 
figure(1);  
for kk = 1:length(K_set)  
    K = K_set(kk);  
    M_set = K:5:N;  
    L_Mset = length(M_set);  
    plot(M_set,Percentage_OMP(kk,1:L_Mset),S1(kk,:));%���x�Ļָ��ź�  
    hold on;  
end  
hold off;  
xlim([0 256]);  
legend('K=4','K=12','K=20','K=28','K=36');  
xlabel('Number of measurements(M)');  
ylabel('Percentage recovered');  
title('Percentage of input signals recovered correctly(OMP)');  


figure(2);  
for kk = 1:length(K_set)  
    K = K_set(kk);  
    M_set = K:5:N;  
    L_Mset = length(M_set);  
    plot(M_set,Percentage_CoSaMP(kk,1:L_Mset),S(kk,:));%���x�Ļָ��ź�  
    hold on;  
end  
hold off;  
xlim([0 256]);  
legend('K=4','K=12','K=20','K=28','K=36');  
xlabel('Number of measurements(M)');  
ylabel('Percentage recovered');  
title('Percentage of input signals recovered correctly(CoSaMP)');


figure(3);  
for kk = 1:length(K_set)  
    K = K_set(kk);  
    M_set = K:5:N;  
    L_Mset = length(M_set);  
    plot(M_set,Percentage_SP(kk,1:L_Mset),S(kk,:));%���x�Ļָ��ź�  
    hold on;  
end  
hold off;  
xlim([0 256]);  
legend('K=4','K=12','K=20','K=28','K=36');  
xlabel('Number of measurements(M)');  
ylabel('Percentage recovered');  
title('Percentage of input signals recovered correctly(SP)');