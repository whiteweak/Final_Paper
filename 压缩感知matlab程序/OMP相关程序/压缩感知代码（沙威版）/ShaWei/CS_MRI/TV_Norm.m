%  ѹ�����еķ����Թ����ݶȷ���ʵ�֣��ó����Ӧ�õ��˴Ź������
%  ����ȫ��ּ����Ϳ��ٸ���Ҷ�任����
%  ���������ʽΪ����Ҷ��Ĳ���
%  �ο����ף�1.  Michael Lustig, David Donoho and John M. Pauly��
%  Sparse MRI: The Application of Compressed Sensing for Rapid MR Imaging.
%  2.  Dai Qi and Wei E.I. Sha, The Physics of Compressive Sensing and the 
%  Gradient-Based Recovery Algorithms.
%  ����ˣ�ɳ�� ��۴�ѧ
%  ���ʱ�䣺2011��5��21��
%  �����ʼ�: wsha@eee.hku.hk
%  ת��ʱ�뱣�������ע��

function TV_Norm
clc;clear

%  ͼ���
A={'phantom256.bmp' 'fruits256.bmp' 'cameraman256.bmp' 'lena256.bmp'...
    'peppers256.bmp' 'boat256.bmp' 'baboon256.bmp'};
I=imread(A{1}); %  1��ʾ����phantomͼ�񣬿��Գ���������ͼ��

%  ͼ���һ��
[a,b]=size(I);
Scale=max(max(double(I)));
M_image=double(I)/Scale;

%  ����Ҷ�����ģʽ������mask_radial.m���Թ��죩
%  ���������ԼΪͼ����������30%
load mask_radial
%  ���ɸ���Ҷ������������measurement��
M_measure=FT_for(mask_matrix,M_image); 

%  ͼ�������ֵ
M_0=zeros(a,b);

%  2�������ܱ���Ȩ����ֵ
lambda=0.01;
%  �ݶ�����
grad=grad_2norm(mask_matrix,M_0,M_measure)+lambda*grad_1norm_tv(M_0);
gradx=grad;
dire=-grad;

%  ��������Line Search������
alpha=0.1;
beta=0.6;
tau0=1;
index_1=0;

%  �ݶ�������׼
epsi=1e-3;
%  ����������
max_iter=100;
%  ��ǰ��������
k=0;
%  ���ջָ���ͼ��
M_recover=M_0;

%  ��������
while(norm(grad,'fro')>epsi && k<max_iter)
    
    %  ��ֵ
    tau=tau0;
    num=0;
    
    %  ��������Line Search��
    while ((f_2norm(mask_matrix,M_recover+tau*dire,M_measure)+...
            lambda*f_1norm_tv(M_recover+tau*dire))>...
           (f_2norm(mask_matrix,M_recover,M_measure)+...
            lambda*f_1norm_tv(M_recover)+alpha*tau*real(conj(grad).*dire)))
        tau=beta*tau;
        num=num+1;
    end
    
    %  ����Ӧ��ֵ
	if num>2
		tau0 = tau0*beta;
	end 
	if num<1
		tau0 = tau0/beta;
    end
    
    %  �ָ�ͼ������
    M_recover=M_recover+tau*dire;
    grad_0=grad;  
    
    %  �ݶ���ʾ
    grad_show=norm(grad,'fro');
    disp('�ݶ���')
    disp(grad_show)
    
    %  ��ԭʼͼ����ȵķ�ֵ�����PSNR
    errorx=sum(sum(abs(M_image-M_recover).^2));  %  MSE���
    psnr=10*log10(1*1/(errorx/256/256));         %  PSNR
    disp('��ֵ����ȣ�')
    disp(psnr)
    
    %  �ݶ�����
    grad=grad_2norm(mask_matrix,M_recover,M_measure)+lambda*grad_1norm_tv(M_recover);
    gamma=norm(grad,'fro')^2/norm(grad_0,'fro')^2;
    dire=-grad+gamma*dire;  
    
    %  ������������
    k=k+1;
    disp('����������')
    disp(k)
    
    %  2�������ܱ���Ȩ����ֵ�������������������lambda�ɹ̶���
    lambda=lambda*0.90;

end

%  �����ʾ
figure(1)
image(abs(M_measure)); 
title('����Ҷ����������')  

figure(2);
subplot(2,2,1)
imshow(uint8(Scale*M_image));
title('ԭʼͼ��')
subplot(2,2,2)
imshow(uint8(Scale*FT_back(mask_matrix,M_measure)));
title('2�����ָ�ͼ��')
subplot(2,2,3)
imshow(uint8(Scale*M_recover));
title('�ܱ��ָ�ͼ��')
subplot(2,2,4)
imshow(uint8(20*abs(Scale*M_image-Scale*M_recover)));  %  �����ȷŴ�20��
title('�ָ����ͼ��')

%  ����Ҷ���任����(�������)
function MM=FT_for(mask,M)
M=real(M);
MM=mask.*(fftshift(fft2(M)));

%  ����Ҷ���任����(��������Ĺ�������)
function MM=FT_back(mask,M)
MM=real(ifft2(ifftshift(mask.*M)));

%  2����
function TT=f_2norm(mask_matrix,T,S)  
TT=norm(FT_for(mask_matrix,T)-S,'fro')^2;

%  �ܱ��
function TT=f_1norm_tv(Solution)
Solution=[Solution(:,1) Solution Solution(:,end)];
Solution=[Solution(1,:);Solution;Solution(end,:)];
df_x=(Solution(2:end-1,3:end)-Solution(2:end-1,1:end-2))/2;
df_y=(Solution(3:end,2:end-1)-Solution(1:end-2,2:end-1))/2;
TT=sum(sum(sqrt(df_x.^2+df_y.^2)));

%  2�������ݶ�
function TT=grad_2norm(mask_matrix,T,S)   
TT=2*FT_back(mask_matrix,FT_for(mask_matrix,T)-S);

%  �ܱ����ݶ�
function TT=grad_1norm_tv(Solution)

epsx=1e-14;  %  ��ֹ�ݶ����޴�

Solution=[Solution(:,1) Solution Solution(:,end)];
Solution=[Solution(1,:);Solution;Solution(end,:)];

%  ������
xx_1=Solution(2:end-1,2:end-1)-Solution(2:end-1,3:end);
yy_1=Solution(2:end-1,2:end-1)-Solution(3:end,2:end-1);
%  ��ߵ���
xx_2=Solution(2:end-1,1:end-2)-Solution(2:end-1,2:end-1);
yy_2=Solution(2:end-1,1:end-2)-Solution(3:end,1:end-2);
%  �ϱߵ���
xx_3=Solution(1:end-2,2:end-1)-Solution(1:end-2,3:end);
yy_3=Solution(1:end-2,2:end-1)-Solution(2:end-1,2:end-1);

%  �ݶ�
grad_1=sqrt(xx_1.^2+yy_1.^2+epsx);
grad_2=sqrt(xx_2.^2+yy_2.^2+epsx);
grad_3=sqrt(xx_3.^2+yy_3.^2+epsx);
            
%  �ܱ����ݶ�
TT=(xx_1./grad_1+yy_1./grad_1-xx_2./grad_2-yy_3./grad_3);



