%  压缩传感的非线性共轭梯度方法实现（该程序可应用到核磁共振成像）
%  采用全变分技术和快速傅立叶变换技术
%  随机测量方式为傅立叶域的采样
%  参考文献：1.  Michael Lustig, David Donoho and John M. Pauly，
%  Sparse MRI: The Application of Compressed Sensing for Rapid MR Imaging.
%  2.  Dai Qi and Wei E.I. Sha, The Physics of Compressive Sensing and the 
%  Gradient-Based Recovery Algorithms.
%  编程人：沙威 香港大学
%  编程时间：2011年5月21日
%  电子邮件: wsha@eee.hku.hk
%  转载时请保留上面的注释

function TV_Norm
clc;clear

%  图像库
A={'phantom256.bmp' 'fruits256.bmp' 'cameraman256.bmp' 'lena256.bmp'...
    'peppers256.bmp' 'boat256.bmp' 'baboon256.bmp'};
I=imread(A{1}); %  1表示最经典的phantom图像，可以尝试其它的图像

%  图像归一化
[a,b]=size(I);
Scale=max(max(double(I)));
M_image=double(I)/Scale;

%  傅立叶域采样模式（程序mask_radial.m可以构造）
%  随机测量数约为图像像素数的30%
load mask_radial
%  生成傅立叶域的随机测量（measurement）
M_measure=FT_for(mask_matrix,M_image); 

%  图像迭代初值
M_0=zeros(a,b);

%  2范数和总变差的权重阈值
lambda=0.01;
%  梯度生成
grad=grad_2norm(mask_matrix,M_0,M_measure)+lambda*grad_1norm_tv(M_0);
gradx=grad;
dire=-grad;

%  线搜索（Line Search）变量
alpha=0.1;
beta=0.6;
tau0=1;
index_1=0;

%  梯度收敛标准
epsi=1e-3;
%  最大迭代次数
max_iter=100;
%  当前迭代次数
k=0;
%  最终恢复的图像
M_recover=M_0;

%  迭代收敛
while(norm(grad,'fro')>epsi && k<max_iter)
    
    %  初值
    tau=tau0;
    num=0;
    
    %  线搜索（Line Search）
    while ((f_2norm(mask_matrix,M_recover+tau*dire,M_measure)+...
            lambda*f_1norm_tv(M_recover+tau*dire))>...
           (f_2norm(mask_matrix,M_recover,M_measure)+...
            lambda*f_1norm_tv(M_recover)+alpha*tau*real(conj(grad).*dire)))
        tau=beta*tau;
        num=num+1;
    end
    
    %  自适应阈值
	if num>2
		tau0 = tau0*beta;
	end 
	if num<1
		tau0 = tau0/beta;
    end
    
    %  恢复图像修正
    M_recover=M_recover+tau*dire;
    grad_0=grad;  
    
    %  梯度显示
    grad_show=norm(grad,'fro');
    disp('梯度误差：')
    disp(grad_show)
    
    %  和原始图像相比的峰值信噪比PSNR
    errorx=sum(sum(abs(M_image-M_recover).^2));  %  MSE误差
    psnr=10*log10(1*1/(errorx/256/256));         %  PSNR
    disp('峰值信噪比：')
    disp(psnr)
    
    %  梯度修正
    grad=grad_2norm(mask_matrix,M_recover,M_measure)+lambda*grad_1norm_tv(M_recover);
    gamma=norm(grad,'fro')^2/norm(grad_0,'fro')^2;
    dire=-grad+gamma*dire;  
    
    %  迭代次数更新
    k=k+1;
    disp('迭代次数：')
    disp(k)
    
    %  2范数和总变差的权重阈值收缩（对有噪声的情况lambda可固定）
    lambda=lambda*0.90;

end

%  结果显示
figure(1)
image(abs(M_measure)); 
title('傅立叶域的随机测量')  

figure(2);
subplot(2,2,1)
imshow(uint8(Scale*M_image));
title('原始图像')
subplot(2,2,2)
imshow(uint8(Scale*FT_back(mask_matrix,M_measure)));
title('2范数恢复图像')
subplot(2,2,3)
imshow(uint8(Scale*M_recover));
title('总变差恢复图像')
subplot(2,2,4)
imshow(uint8(20*abs(Scale*M_image-Scale*M_recover)));  %  误差幅度放大20倍
title('恢复误差图像')

%  傅立叶正变换算子(随机测量)
function MM=FT_for(mask,M)
M=real(M);
MM=mask.*(fftshift(fft2(M)));

%  傅立叶反变换算子(随机测量的共轭算子)
function MM=FT_back(mask,M)
MM=real(ifft2(ifftshift(mask.*M)));

%  2范数
function TT=f_2norm(mask_matrix,T,S)  
TT=norm(FT_for(mask_matrix,T)-S,'fro')^2;

%  总变差
function TT=f_1norm_tv(Solution)
Solution=[Solution(:,1) Solution Solution(:,end)];
Solution=[Solution(1,:);Solution;Solution(end,:)];
df_x=(Solution(2:end-1,3:end)-Solution(2:end-1,1:end-2))/2;
df_y=(Solution(3:end,2:end-1)-Solution(1:end-2,2:end-1))/2;
TT=sum(sum(sqrt(df_x.^2+df_y.^2)));

%  2范数的梯度
function TT=grad_2norm(mask_matrix,T,S)   
TT=2*FT_back(mask_matrix,FT_for(mask_matrix,T)-S);

%  总变差的梯度
function TT=grad_1norm_tv(Solution)

epsx=1e-14;  %  防止梯度无限大

Solution=[Solution(:,1) Solution Solution(:,end)];
Solution=[Solution(1,:);Solution;Solution(end,:)];

%  自身导数
xx_1=Solution(2:end-1,2:end-1)-Solution(2:end-1,3:end);
yy_1=Solution(2:end-1,2:end-1)-Solution(3:end,2:end-1);
%  左边导数
xx_2=Solution(2:end-1,1:end-2)-Solution(2:end-1,2:end-1);
yy_2=Solution(2:end-1,1:end-2)-Solution(3:end,1:end-2);
%  上边导数
xx_3=Solution(1:end-2,2:end-1)-Solution(1:end-2,3:end);
yy_3=Solution(1:end-2,2:end-1)-Solution(2:end-1,2:end-1);

%  梯度
grad_1=sqrt(xx_1.^2+yy_1.^2+epsx);
grad_2=sqrt(xx_2.^2+yy_2.^2+epsx);
grad_3=sqrt(xx_3.^2+yy_3.^2+epsx);
            
%  总变差的梯度
TT=(xx_1./grad_1+yy_1./grad_1-xx_2./grad_2-yy_3./grad_3);



