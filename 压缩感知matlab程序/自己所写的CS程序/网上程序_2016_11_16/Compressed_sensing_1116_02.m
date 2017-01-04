% sparse_in_frequency.m

%
%This code demonstrate compressive sensing example. In this
%example the signal is sparse in frequency domain and random samples
%are taken in time domain.

close all;
clear all;

%setup path for the subdirectories of l1magic
% path(path, 'C:\MATLAB7\l1magic-1.1\Optimization');
% path(path, 'C:\MATLAB7\l1magic-1.1\Data');


%length of the signal
N=1024;

%Number of random observations to take
K=256;

%Discrete frequency of two sinusoids in the input signal
k1=29;
k2=100;

n=0:N-1;

%Sparse signal in frequency domain.
x=sin(2*pi*(k1/N)*n)+sin(2*pi*(k2/N)*n);

% This code demonstrates the compressive sensing using a sparse signal in frequency domain. The signal consists of summation of % two sinusoids of different frequencies in time domain. The signal is sparse in Frequency domain and therefore K random

% measurements are taken in time domain.


figure;
subplot(2,1,1);
plot(x)
grid on;
xlabel('Samples');
ylabel('Amplitude');
title('Original Signal,1024 samples with two different frequency sinsuoids');

xf=fft(x);

subplot(2,1,2);
plot(abs(xf))
grid on;
xlabel('Samples');
ylabel('Amplitude');
title('Frequency domain, 1024 coefficients with 4-non zero coefficients');

%creating dft matrix
B=dftmtx(N);
Binv=inv(B);  % The inverse discrete Fourier transform matrix, Binv, equals CONJ(dftmtx(N))/N.

%Taking DFT of the signal
xf = B*x.';

%Selecting random rows of the DFT matrix
q=randperm(N);

%creating measurement matrix
A=Binv(q(1:K),:);    % ��IDFT��������ѡK=256��

%taking random time measurements
y=(A*xf);   % ��x��fft���xf(1024-by-1)��������IDFT�õ�256��ʱ��ϡ�����ֵ��ͨ��plot(real(y))��ԭ����x�Աȣ�ע�������ʱ����ȡK=256������ֵ

%Calculating Initial guess
x0=A.'*y;  % ע�⣺���ָ�ʱ���ź�xprec��DFTֵxp�Ĺ��Ƴ�ֵx0��θ��� y ��ʱ��ϡ�����ֵ

%Running the recovery Algorithm
tic
xp=l1eq_pd(x0,A,[],y,1e-5); %�ָ���xp��Ƶ���ź�
toc

%recovered signal in time domain
xprec=real(Binv*xp);  % ��IDFTת����ʱ��

figure;
subplot(2,1,1)
plot(abs(xf))   % ԭ�źŵ�Ƶ��
grid on;
xlabel('Samples');
ylabel('Amplitude');
title('Original Signal, Discrete Fourier Transform');

subplot(2,1,2)
plot(abs(xp),'r') %ѹ�������ָ�����źŵ�Ƶ��
grid on;
xlabel('Samples');
ylabel('Amplitude');
title(sprintf('Recovered Signal, Discrete Fourier Transform sampled with %d samples',K));

figure;
subplot(2,1,1);
plot(x)
grid on;
xlabel('Samples');
ylabel('Amplitude');
title('Original Signal,1024 samples with two different frequency sinsuoids');

subplot(2,1,2)
plot(xprec,'r')
grid on;
xlabel('Samples');
ylabel('Amplitude');
title(sprintf('Recovered Signal in Time Domain'));