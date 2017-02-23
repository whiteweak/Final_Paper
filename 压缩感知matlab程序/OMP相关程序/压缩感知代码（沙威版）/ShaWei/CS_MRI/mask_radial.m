%  �ó������ɸ���Ҷ�任������߲���ģʽ�����������
%  ����ˣ�ɳ�� ��۴�ѧ
%  ���ʱ�䣺2011��5��21��
%  �����ʼ�: wsha@eee.hku.hk
%  ת��ʱ�뱣�������ע��

clc;clear

N=256;   %  ͼ���С
M=76;    %  ����Ƕ��� 76
T=linspace(0,pi,M);     %  �����Ƕ�
mask_matrix=zeros(N,N); %  ��ģ����0��1�ֲ���

%  �Ƕ�ѭ��
for i=1:M;  
    kk=tan(T(i));    %  б��
    c1=sqrt(kk^2+1); %  �������ڼ���㵽ֱ�߾���
    %  ����ѭ��
    for m=-N/2+0.5:N/2-0.5;
        for n=-N/2+0.5:N/2-0.5
            mm=m+N/2+0.5;
            nn=n+N/2+0.5;
            d=abs(kk*m-n)/c1;          %  ���빫ʽ
            if (d<1/2)
                mask_matrix(mm,nn)=1;  %  ��ģ����1
            end
        end
    end
end

%  ����ռͼ��ߴ����
disp('��������')
disp(sum(sum(mask_matrix))/(N*N))
%  ����ģʽ
imshow(255*uint8(mask_matrix))
%  ����
save mask_radial mask_matrix