clear all;
clc;
close all;

N=1000;
f1=0.1;
f2=0.25;
f3=0.27;
p=2*pi*rand(1,3);
n=1:N;
s1(n)=sqrt(10^3)*exp(j*2*pi*f1*n+j*p(1));
s2(n)=sqrt(10^3)*exp(j*2*pi*f2*n+j*p(2));
s3(n)=sqrt(10^2.7)*exp(j*2*pi*f3*n+j*p(3));
v=sqrt(1/2)*randn(1,N)+j*sqrt(1/2)*randn(1,N);
x=s1+s2+s3+v;


for i=1:N-3
    A(i,:)=[x(i+3),x(i+2),x(i+1),x(i)];
end
B=A.'*A'.';
[V D]=eig(B);
D=eig(D);



ifai=zeros(4,4);
 for i=1:4
    ifai=ifai+V(:,i)*V(:,i)'/D(i);
 end
 
 f=-0.5:0.0005:0.5;
 L=length(f);
 for i=1:L
    a=[1;exp(-j*2*pi*f(i));exp(-2*j*2*pi*f(i));exp(-3*j*2*pi*f(i))];
    p_MVDR(i)=1/(a'*ifai*a);
 end

 
 p_MVDR=p_MVDR/max(abs(p_MVDR));
 p_MVDR=10*log10(abs(p_MVDR));
 plot(f,p_MVDR);
 xlabel('w/2/pi');ylabel('归一化MVDR谱/dB');
 title('MVDR方法');
   



