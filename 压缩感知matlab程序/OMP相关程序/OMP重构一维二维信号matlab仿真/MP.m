clc;clear
%%  1. 时域测试信号生成
%产生长度为N=256的稀疏信号，其稀疏度K=23且这23个非零值随机分布于信号256个位置
%观测向量y的长度M=80，即采样率M/N=0.3
N=256;
K=23;
M=80;
x = zeros(N,1);
q = randperm(N);
x(q(1:K)) =randn(K,1);    %原始信号
%% 2. 测量矩阵 及观测值获得
Phi=randn(M,N); %测量矩阵 %  感知矩阵（高斯分布白噪声）M*N
matrixNorm = Phi.'*Phi;
matrixNorm = sqrt(diag(matrixNorm)).';
Phi = Phi./repmat(matrixNorm, [M,1]);  %注意，观测矩阵是要归一化的，因为原子范数要是1！
y=Phi*x ;       %获得线性测量
%% 3.用MP算法重构信号
iterations=K;                                      %  算法迭代次数(m>=K)
%signal_reconstruct=zeros(1,1);                    %  近似解矩阵(初始值为空矩阵)
r_n=y;     %  残差值M*1
x_rec=zeros(N,1);
for times=1:iterations
    for col=1:N                             %感知矩阵的所有列向量
        innerpro(col)=Phi(:,col)'*r_n;      %计算余量和感知矩阵每一列的内积
    end
    [val,pos]=max(abs(innerpro) );           %找出内积中绝对值最大的元素和它的对应的感知矩阵的列pos  
    x_rec(pos)=x_rec(pos)+innerpro(pos);    %计算新的近似x_rec
    r_n=r_n-innerpro(pos)*Phi(:,pos);       %更新残差       
end
norm(x_rec-x)/norm(x)                           %  重构误差
subplot(3,1,1);plot(x);title('origin');
subplot(3,1,2);plot(x_rec);title('reconstruct');
subplot(3,1,3);plot(r_n);title('残差');
