function res = RESIDUE(H)
%RESIDUE ͨ�������������������ܺ�����Χ�߻���
%   ����Ȩֵ[b1,b2],���Ϊ��Ȩֵ��Ե�Χ�߻��֣�
%a0=1;b1=1;b2=2;
% syms z
% a0=1;
% H=(a0/(1-b1*z-b2*z^2)-2/(1-1.2*z+0.6*z^2))*a0*z/(z^2-b1*z-b2);
[N,D]=numden(H) ;
p1=sym2poly(N);%%%%��ȡ���Ӷ���ʽϵ��
p2=sym2poly(D);%%%%��ȡ��ĸ����ʽϵ��
[R,P,K]=residue(p1,p2);%%%%%%�õ������㼰�����㴦������
[b,index1]=unique(P,'first');%%%%%�ҵ���ͬ�ļ��㼰��Ӧ����
index2=find(abs(P(index1))<1);%%%ȥ��ģ����1�ļ���
index=index1(index2);        %%%�õ�ģС��1�ļ����λ��
res=sum(R(index));
% res=real(res);
end


% function res = RESIDUE(b1,b2)
% %RESIDUE ͨ�������������������ܺ�����Χ�߻���
% %   ����Ȩֵ[b1,b2],���Ϊ��Ȩֵ��Ե�Χ�߻��֣�
% %a0=1;b1=1;b2=2;
% syms z
% a0=1;
% H=(a0/(1-b1*z-b2*z^2)-2/(1-1.2*z+0.6*z^2))*a0*z/(z^2-b1*z-b2);
% [N,D]=numden(H) ;
% p1=sym2poly(N);%%%%��ȡ���Ӷ���ʽϵ��
% p2=sym2poly(D);%%%%��ȡ��ĸ����ʽϵ��
% [R,P,K]=residue(p1,p2);%%%%%%�õ������㼰�����㴦������
% [b,index1]=unique(P,'first');%%%%%�ҵ���ͬ�ļ��㼰��Ӧ����
% index2=find(abs(P(index1))<1);%%%ȥ��ģ����1�ļ���
% index=index1(index2);        %%%�õ�ģС��1�ļ����λ��
% res=sum(R(index));
% end
% %
