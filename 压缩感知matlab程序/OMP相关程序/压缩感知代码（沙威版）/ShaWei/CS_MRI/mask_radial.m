%  该程序生成傅里叶变换域的射线采样模式用作随机测量
%  编程人：沙威 香港大学
%  编程时间：2011年5月21日
%  电子邮件: wsha@eee.hku.hk
%  转载时请保留上面的注释

clc;clear

N=256;   %  图像大小
M=76;    %  随机角度数 76
T=linspace(0,pi,M);     %  测量角度
mask_matrix=zeros(N,N); %  掩模矩阵（0，1分布）

%  角度循环
for i=1:M;  
    kk=tan(T(i));    %  斜率
    c1=sqrt(kk^2+1); %  常数用于计算点到直线距离
    %  像素循环
    for m=-N/2+0.5:N/2-0.5;
        for n=-N/2+0.5:N/2-0.5
            mm=m+N/2+0.5;
            nn=n+N/2+0.5;
            d=abs(kk*m-n)/c1;          %  距离公式
            if (d<1/2)
                mask_matrix(mm,nn)=1;  %  掩模设置1
            end
        end
    end
end

%  测量占图像尺寸比例
disp('测量比例')
disp(sum(sum(mask_matrix))/(N*N))
%  测量模式
imshow(255*uint8(mask_matrix))
%  保存
save mask_radial mask_matrix