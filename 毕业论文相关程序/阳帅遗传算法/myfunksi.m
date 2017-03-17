function ksi = myfunksi(W)
%UNTITLED7 计算性能函数
%   W为输入权值，ksi为均方误差
a0=W(:,1);
% a0=1;
b1=W(:,2);b2=W(:,3);
rt1=0.5*(b1+sqrt(b1.^2+4*b2));
rt2=0.5*(b1-sqrt(b1.^2+4*b2));

% res1=rt1./((1-rt1.*rt1).*(1-rt1.*rt2).*(rt1-rt2));%%计算rt1处的留数
% res2=rt2./((1-rt2.*rt1).*(1-rt2.*rt2).*(rt2-rt1));%%计算rt2处的留数
res1=a0.^2.*rt1./((1-rt1.*rt1).*(1-rt1.*rt2).*(rt1-rt2));%%计算rt1处的留数
res2=a0.^2.*rt2./((1-rt2.*rt1).*(1-rt2.*rt2).*(rt2-rt1));%%计算rt2处的留数

%PART2---该式为-2*phidx(z)*H(z)/z
% res3=-2*rt1./((1-1.2*rt1+0.6*rt1.^2).*(rt1-rt2));
% res4=-2*rt2./((1-1.2*rt2+0.6*rt2.^2).*(rt2-rt1));
res3=-2*a0.*rt1./((1-1.2*rt1+0.6*rt1.^2).*(rt1-rt2));
res4=-2*a0.*rt2./((1-1.2*rt2+0.6*rt2.^2).*(rt2-rt1));


ksi=3.5714+res1+res2+res3+res4;     %%%%性能表面函数

end

