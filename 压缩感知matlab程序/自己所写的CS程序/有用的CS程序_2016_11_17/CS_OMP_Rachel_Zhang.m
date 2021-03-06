

%  1-D信号压缩传感的实现(正交匹配追踪法Orthogonal Matching Pursuit)
%  测量数M>=m*log(d/m),K是稀疏度,N是信号长度,可以近乎完全重构原始信号

clc;clear

%%  ****************** 时域测试信号生成 ****************** 
m=7;      %  稀疏度设为K(“稀疏度”可理解为信号中非零元素或系数较大的元素的个数，做FFT可以看出来)，相当于文章中的m
d=256;    %  信号长度，取的是离散信号,相当于文章中的d
N=64;     %  测量数(N>=m*log(d/m),至少40,但考虑到有出现错误的概率，故取大一些，作为安全余量)
f1=50;    %  信号频率，单位Hz,这里采用了四个信号叠加起来的方式构成了输入信号
f2=100;   
f3=200;   
f4=400;   
fs=100;   %  采样频率设为800，保证满足奈奎斯特频率(?),这里也可采用一些其他的数值
ts=1/fs;  %  采样间隔
Ts=1:d;   %  采样序列长度
s=0.3*cos(2*pi*f1*Ts*ts)+0.6*cos(2*pi*f2*Ts*ts)+0.1*cos(2*pi*f3*Ts*ts)-0.9*cos(2*pi*f4*Ts*ts);  
          %  完整信号s,由4个余弦信号叠加而来

%%  ****************** 时域信号压缩传感 ****************** 
Phi=randn(N,d);                 %  测量矩阵Phi(高斯分布白噪声)64*256的扁矩阵,这里是先设好了测量矩阵然后再进行相关计算 
v=Phi*s.';                      %  获得线性测量 ，s相当于文中的y矩阵

%%  ****************** 正交匹配追踪(OMP)法重构信号(本质上是L_1范数最优化问题) ****************** 
%匹配追踪：找到一个其标记看上去与收集到的数据相关（內积最大）的小波；在信号数据中去除这个标记的所有印迹；
%不断重复直到我们能用小波标记“解释”收集到的所有数据。
%本程序的做法是：构造一个傅里叶空间T，其列向量的维数和原始信号的维数相同为Nx1，这个空间由一个同样大小为
%64x256维的字典矩阵构成，匹配追踪的思路是在该字典矩阵中寻找一个与信号（列向量）x最匹配的列向量，来构建一个稀疏逼近，
%并计算信号向量x与该逼近列向量之间的差值记为残差（初始值设置为测量到的信号y，意为误差最大的时候），然后继续选择与信号
%残差最匹配（接近）的列向量，如此循环下去（迭代次数由稀疏度来决定，一般为其两倍），则原始信号x可以由这些作为稀疏逼近的
%列向量来线性和再加上最后的残差值来决定。如果残差值在可以忽略的范围内，则信号x就是这些原子的线性组合。

t=2*m;                                            %  算法迭代次数(t>=m)，设x是K-sparse（K稀疏）的
Psi=fft(eye(d,d));                        %  傅里叶正变换矩阵，eye(d,d)函数返回的是N*N的单位矩阵
Psi=Psi/sqrt(d);
T=Phi*Psi';                                       %  恢复矩阵(测量矩阵*正交反变换矩阵)，为64x256维

hat_y=zeros(1,d);                                 %  待重构的谱域(变换域)向量，为1x256维                     
Aug_t=[];                                         %  增量矩阵(初始值为空矩阵)
r_m=v;                                            %  初始化残差值矩阵，r_m为64x1维

for times=1:t;                                    %  迭代次数(有噪声的情况下,该迭代次数为K)
    for col=1:d;                                  %  恢复矩阵的所有列向量
        product(col)=abs(T(:,col)'*r_m);          %  求恢复矩阵的列向量和残差的投影系数即他们之间的内积值 
    end
    [val,pos]=max(product);       %  最大投影系数对应的位置，即找到一个其标记看上去与收集到的数据相关的小波
                                  %  val得到的是product这个1x256的內积值矩阵中最大值的-value，而pos得到的是
                                  %  最大值的位置-position。此时，我们找到了第1次迭代时信号中与r_n最相关的
                                  %  （关联最大的）小波成分即找到了恢复矩阵T的列向量中与r_n內积最大的那一列
    Aug_t=[Aug_t,T(:,pos)];       %  将上面得到的关联最大的小波（T的某一列向量）加入Aug_t增量矩阵，进行矩阵扩充    
    
    T(:,pos)=zeros(N,1);          %  将刚才得到的与残差矩阵r_n关联最大的那一列置零（实质上应该去掉，为了简单我把
                                  %  它置零），在数据中去除这个标记（小波）的所有印迹
    aug_y=(Aug_t'*Aug_t)^(-1)*Aug_t'*v;           %  最小二乘,使残差值最小，实际上相当于：aug_y=(1/(A'A))A'v
                                                  %  =A^(-1)Iv=A^(-1)v
    r_m=v-Aug_t*aug_y;                            %  计算本次迭代后的新的残差值
    pos_array(times)=pos;                         %  纪录最大投影系数(內积值)的位置
end
hat_y(pos_array)=aug_y;                           %  重构的谱域向量
hat_x=real(Psi'*hat_y.');                         %  做逆傅里叶变换重构得到时域信号

%%  ****************** 恢复信号和原始信号对比 ****************** 
figure(1);
hold on;
plot(hat_x,'k:*')                                 %  重建信号
plot(s,'r-o')                                       %  原始信号
legend('Recovery','Original')
error=norm(hat_x.'-s)/norm(s);                           %  重构误差