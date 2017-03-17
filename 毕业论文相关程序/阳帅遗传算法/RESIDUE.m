function res = RESIDUE(H)
%RESIDUE 通过计算留数来计算性能函数的围线积分
%   输入权值[b1,b2],输出为该权值点对的围线积分；
%a0=1;b1=1;b2=2;
% syms z
% a0=1;
% H=(a0/(1-b1*z-b2*z^2)-2/(1-1.2*z+0.6*z^2))*a0*z/(z^2-b1*z-b2);
[N,D]=numden(H) ;
p1=sym2poly(N);%%%%提取分子多项式系数
p2=sym2poly(D);%%%%提取分母多项式系数
[R,P,K]=residue(p1,p2);%%%%%%得到各极点及各极点处的留数
[b,index1]=unique(P,'first');%%%%%找到不同的极点及对应索引
index2=find(abs(P(index1))<1);%%%去除模大于1的极点
index=index1(index2);        %%%得到模小于1的极点的位置
res=sum(R(index));
% res=real(res);
end


% function res = RESIDUE(b1,b2)
% %RESIDUE 通过计算留数来计算性能函数的围线积分
% %   输入权值[b1,b2],输出为该权值点对的围线积分；
% %a0=1;b1=1;b2=2;
% syms z
% a0=1;
% H=(a0/(1-b1*z-b2*z^2)-2/(1-1.2*z+0.6*z^2))*a0*z/(z^2-b1*z-b2);
% [N,D]=numden(H) ;
% p1=sym2poly(N);%%%%提取分子多项式系数
% p2=sym2poly(D);%%%%提取分母多项式系数
% [R,P,K]=residue(p1,p2);%%%%%%得到各极点及各极点处的留数
% [b,index1]=unique(P,'first');%%%%%找到不同的极点及对应索引
% index2=find(abs(P(index1))<1);%%%去除模大于1的极点
% index=index1(index2);        %%%得到模小于1的极点的位置
% res=sum(R(index));
% end
% %
