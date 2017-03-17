

%压缩感知重构算法测试CS_Reconstuction_MtoPercentage.m
%不同算法的信号测量数与恢复百分比之间的关系

clear all;close all;clc;  

S1 = ['-ks';'-ko';'-kd';'-kv';'-k*'];  
S2 = ['-rs';'-ro';'-rd';'-rv';'-r*'];
S3 = ['-bs';'-bo';'-bd';'-bv';'-b*']; 
%% 参数配置初始化  
CNT = 1000;%对于每组(K,M,N)，重复迭代次数  
N = 256;%信号x的长度  
Psi = eye(N);%x本身是稀疏的，定义稀疏矩阵为单位阵x=Psi*theta  
K_set = [4,12,20,28,36];%信号x的稀疏度集合  
Percentage_OMP = zeros(length(K_set),N);%存储恢复成功概率
Percentage_CoSaMP = zeros(length(K_set),N);
Percentage_SP = zeros(length(K_set),N);
%% 主循环，遍历每组(K,M,N)  
tic  
for kk = 1:length(K_set)  
    K = K_set(kk);%本次稀疏度  
    M_set = K:5:N;%M没必要全部遍历，每隔5测试一个就可以了 
    
    PercentageK_OMP = zeros(1,length(M_set));%存储此稀疏度K下不同M的恢复成功概率
    PercentageK_CoSaMP = zeros(1,length(M_set));
    PercentageK_SP = zeros(1,length(M_set));
    
    for mm = 1:length(M_set)  
       M = M_set(mm);%本次观测值个数  
       fprintf('K=%d,M=%d\n',K,M);
       
       %统计恢复成功的信号点数
       P_OMP=0;
       P_CoSaMP = 0;
       P_SP=0;
       
       for cnt = 1:CNT %每个观测值个数均运行CNT次  
            Index_K = randperm(N);  
            x = zeros(N,1);  
            x(Index_K(1:K)) = 5*randn(K,1);%x为K稀疏的，且位置是随机的
            
            %% OMP方法信号恢复及恢复误差计算过程 
            Phi = randn(M,N);%测量矩阵为高斯矩阵  
            A = Phi * Psi;%传感矩阵  
            y = Phi * x;%得到观测向量y  
            theta = CS_OMP(y,A,K);%恢复重构信号theta  
            x_r = Psi * theta;% x=Psi * theta  
            if norm(x_r-x)<1e-6%如果残差小于1e-6则认为恢复成功  
                P_OMP = P_OMP + 1;
            end
                
            %% CoSaMP方法信号恢复及恢复误差计算过程 
            Phi = randn(M,N)/sqrt(M);%测量矩阵为高斯矩阵  
            A = Phi * Psi;%传感矩阵  
            y = Phi * x;%得到观测向量y  
            theta = CS_CoSaMP(y,A,K);%恢复重构信号theta  
            x_r = Psi * theta;% x=Psi * theta  
            if norm(x_r-x)<1e-6%如果残差小于1e-6则认为恢复成功  
                P_CoSaMP = P_CoSaMP + 1;      
                
            end  
            
            %% SP方法信号恢复及恢复误差计算过程 
            Phi = randn(M,N)/sqrt(M);%测量矩阵为高斯矩阵  
            A = Phi * Psi;%传感矩阵  
            y = Phi * x;%得到观测向量y  
            theta = CS_SP(y,A,K);%恢复重构信号theta  
            x_r = Psi * theta;% x=Psi * theta  
            if norm(x_r-x)<1e-6%如果残差小于1e-6则认为恢复成功  
                P_SP = P_SP + 1;      
                
            end  
       end  
       
       PercentageK_OMP(mm) = P_OMP/CNT*100;%计算恢复概率 
       PercentageK_CoSaMP(mm) = P_CoSaMP/CNT*100;
       PercentageK_SP(mm) = P_SP/CNT*100;
       
    end  
    Percentage_OMP(kk,1:length(M_set)) = PercentageK_OMP;
    Percentage_CoSaMP(kk,1:length(M_set)) = PercentageK_CoSaMP;
    Percentage_SP(kk,1:length(M_set)) = PercentageK_SP;
end  
toc  
save MtoPercentage1000 %运行一次不容易，把变量全部存储下来  
%% 绘图  
 
figure(1);  
for kk = 1:length(K_set)  
    K = K_set(kk);  
    M_set = K:5:N;  
    L_Mset = length(M_set);  
    plot(M_set,Percentage_OMP(kk,1:L_Mset),S1(kk,:));%绘出x的恢复信号  
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
    plot(M_set,Percentage_CoSaMP(kk,1:L_Mset),S(kk,:));%绘出x的恢复信号  
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
    plot(M_set,Percentage_SP(kk,1:L_Mset),S(kk,:));%绘出x的恢复信号  
    hold on;  
end  
hold off;  
xlim([0 256]);  
legend('K=4','K=12','K=20','K=28','K=36');  
xlabel('Number of measurements(M)');  
ylabel('Percentage recovered');  
title('Percentage of input signals recovered correctly(SP)');