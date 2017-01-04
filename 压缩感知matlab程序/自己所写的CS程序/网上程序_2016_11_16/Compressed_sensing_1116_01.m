% sparse_in_time.m
%
%This code demonstrate compressive sensing example. In this
%example the signal is sparse in time domain and random samples
%are taken in frequency domain.

close all;
clear all;

%setup path for the subdirectories of l1magic
% path(path, 'C:\MATLAB7\l1magic-1.1\Optimization');
% path(path, 'C:\MATLAB7\l1magic-1.1\Data');

%number of samples per period
s=4;

%RF frequency
f=4e9;

%pulse repetition frequency
prf=1/30e-9;

%sampling frequency
fs=s*f;

%Total Simulation time
T=30e-9;

t=0:1/fs:T;

%generating pulse train
x=pulstran(t,15e-9,'gauspuls',f,0.5);

%length of the signal
N=length(x);

%Number of random observations to take
K=90;

figure;
subplot(2,1,1);
plot(t,x)
grid on;
xlabel('Time');
ylabel('Amplitude');
title(sprintf('Original Signal, UWB Pulse RF freq=%g GHz',f/1e9));

%taking Discrete time Fourier Transform of the signal
xf=fft(x);

subplot(2,1,2);
plot(abs(xf))
grid on;
xlabel('Samples');
ylabel('Amplitude');
title('Discrete Fourier Transform of UWB pulse');


%creating dft matrix
B=dftmtx(N);
Binv=inv(B);

%Selecting random rows of the DFT matrix
q=randperm(N);

%creating measurement matrix
A=B(q(1:K),:); % ��DFT������ȡǰK=90�У�B������481-by-481��

%taking random frequency measurements
y=(A*x.');  % y ��90-by-1��Ƶ�����ֵ

% Calculating Initial guess
x0=A.'*y;  % y ��90-by-1��Ƶ�����ֵ��ע��ָ�ʱ���ź�ʱ��ֵ��θ���

%Running the recovery Algorithm
tic
xp=l1eq_pd(x0,A,[],y,1e-5);
toc

figure;
subplot(2,1,1)
plot(t,x)
grid on;
xlabel('Time');
ylabel('Amplitude');
title(sprintf('Original Signal, UWB Pulse RF freq=%g GHz',f/1e9));


subplot(2,1,2)
plot(t,real(xp),'r')
grid on;
xlabel('Time');
ylabel('Amplitude');
title(sprintf('Recovered UWB Pulse Signal with %d random samples',K));