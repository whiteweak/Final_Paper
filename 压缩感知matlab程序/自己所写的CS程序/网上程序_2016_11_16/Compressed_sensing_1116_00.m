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
A=Binv(q(1:K),:);    % 在IDFT矩阵中任选K=256行

%taking random time measurements
y=(A*xf);   % 对x的fft后的xf(1024-by-1)的数据做IDFT得到256个时域稀疏采样值，通过plot(real(y))和原来的x对比，注意如何在时域中取K=256个采样值

%Calculating Initial guess
x0=A.'*y;  % 注意：待恢复时域信号xprec的DFT值xp的估计初值x0如何给？ y 是时域稀疏采样值

%Running the recovery Algorithm
tic
xp=l1eq_pd(x0,A,[],y,1e-5); %恢复的xp是频域信号
toc

%recovered signal in time domain
xprec=real(Binv*xp);  % 做IDFT转换到时域

figure;
subplot(2,1,1)
plot(abs(xf))   % 原信号的频谱
grid on;
xlabel('Samples');
ylabel('Amplitude');
title('Original Signal, Discrete Fourier Transform');

subplot(2,1,2)
plot(abs(xp),'r') %压缩采样恢复后的信号的频谱
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

%%%%%%%%%%%%%%%%%%漫长的分割线%%%%%%%%%%%%%%%%%%%%%%%%%%


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
A=B(q(1:K),:); % 在DFT矩阵中取前K=90行，B矩阵是481-by-481的

%taking random frequency measurements
y=(A*x.');  % y 是90-by-1的频域采样值

% Calculating Initial guess
x0=A.'*y;  % y 是90-by-1的频域采样值，注意恢复时域信号时初值如何给？

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

%%%%%%%%%%%%%%%%%%飘逸的分割线%%%%%%%%%%%%%%%%%%%%%%%%%%

% 用到的函数
% l1eq_pd.m
%
% Solve
% min_x ||x||_1  s.t.  Ax = b
%
% Recast as linear program
% min_{x,u} sum(u)  s.t.  -u <= x <= u,  Ax=b
% and use primal-dual interior point method
%
% Usage: xp = l1eq_pd(x0, A, At, b, pdtol, pdmaxiter, cgtol, cgmaxiter)
%
% x0 - Nx1 vector, initial point.
%
% A - Either a handle to a function that takes a N vector and returns a K 
%     vector , or a KxN matrix.  If A is a function handle, the algorithm
%     operates in "largescale" mode, solving the Newton systems via the
%     Conjugate Gradients algorithm.
%
% At - Handle to a function that takes a K vector and returns an N vector.
%      If A is a KxN matrix, At is ignored.
%
% b - Kx1 vector of observations.
%
% pdtol - Tolerance for primal-dual algorithm (algorithm terminates if
%     the duality gap is less than pdtol).  
%     Default = 1e-3.
%
% pdmaxiter - Maximum number of primal-dual iterations.  
%     Default = 50.
%
% cgtol - Tolerance for Conjugate Gradients; ignored if A is a matrix.
%     Default = 1e-8.
%
% cgmaxiter - Maximum number of iterations for Conjugate Gradients; ignored
%     if A is a matrix.
%     Default = 200.
%
% Written by: Justin Romberg, Caltech
% Email: jrom@acm.caltech.edu
% Created: October 2005
%

function xp = l1eq_pd(x0, A, At, b, pdtol, pdmaxiter, cgtol, cgmaxiter)

largescale = isa(A,'function_handle');

if (nargin < 5), pdtol = 1e-3;  end
if (nargin < 6), pdmaxiter = 50;  end
if (nargin < 7), cgtol = 1e-8;  end
if (nargin < 8), cgmaxiter = 200;  end

N = length(x0);

alpha = 0.01;
beta = 0.5;
mu = 10;

gradf0 = [zeros(N,1); ones(N,1)];

% starting point --- make sure that it is feasible
if (largescale)
  if (norm(A(x0)-b)/norm(b) > cgtol)
    disp('Starting point infeasible; using x0 = At*inv(AAt)*y.');
    AAt = @(z) A(At(z));
    [w, cgres, cgiter] = cgsolve(AAt, b, cgtol, cgmaxiter, 0);
    if (cgres > 1/2)
      disp('A*At is ill-conditioned: cannot find starting point');
      xp = x0;
      return;
    end
    x0 = At(w);
  end
else
  if (norm(A*x0-b)/norm(b) > cgtol)
    disp('Starting point infeasible; using x0 = At*inv(AAt)*y.');
    opts.POSDEF = true; opts.SYM = true;
    [w, hcond] = linsolve(A*A', b, opts);
    if (hcond < 1e-14)
      disp('A*At is ill-conditioned: cannot find starting point');
      xp = x0;
      return;
    end
    x0 = A'*w;
  end  
end
x = x0;
u = (0.95)*abs(x0) + (0.10)*max(abs(x0));

% set up for the first iteration
fu1 = x - u;
fu2 = -x - u;
lamu1 = -1./fu1;
lamu2 = -1./fu2;
if (largescale)
  v = -A(lamu1-lamu2);
  Atv = At(v);
  rpri = A(x) - b;
else
  v = -A*(lamu1-lamu2);
  Atv = A'*v;
  rpri = A*x - b;
end

sdg = -(fu1'*lamu1 + fu2'*lamu2);
tau = mu*2*N/sdg;

rcent = [-lamu1.*fu1; -lamu2.*fu2] - (1/tau);
rdual = gradf0 + [lamu1-lamu2; -lamu1-lamu2] + [Atv; zeros(N,1)];
resnorm = norm([rdual; rcent; rpri]);

pditer = 0;
done = (sdg < pdtol) | (pditer >= pdmaxiter);
while (~done)
  
  pditer = pditer + 1;
  
  w1 = -1/tau*(-1./fu1 + 1./fu2) - Atv;
  w2 = -1 - 1/tau*(1./fu1 + 1./fu2);
  w3 = -rpri;
  
  sig1 = -lamu1./fu1 - lamu2./fu2;
  sig2 = lamu1./fu1 - lamu2./fu2;
  sigx = sig1 - sig2.^2./sig1;
  
  if (largescale)
    w1p = w3 - A(w1./sigx - w2.*sig2./(sigx.*sig1));
    h11pfun = @(z) -A(1./sigx.*At(z));
    [dv, cgres, cgiter] = cgsolve(h11pfun, w1p, cgtol, cgmaxiter, 0);
    if (cgres > 1/2)
      disp('Cannot solve system.  Returning previous iterate.  (See Section 4 of notes for more information.)');
      xp = x;
      return
    end
    dx = (w1 - w2.*sig2./sig1 - At(dv))./sigx;
    Adx = A(dx);
    Atdv = At(dv);
  else
    w1p = -(w3 - A*(w1./sigx - w2.*sig2./(sigx.*sig1)));
    H11p = A*(sparse(diag(1./sigx))*A');
    opts.POSDEF = true; opts.SYM = true;
    [dv,hcond] = linsolve(H11p, w1p);
    if (hcond < 1e-14)
      disp('Matrix ill-conditioned.  Returning previous iterate.  (See Section 4 of notes for more information.)');
      xp = x;
      return
    end
    dx = (w1 - w2.*sig2./sig1 - A'*dv)./sigx;
    Adx = A*dx;
    Atdv = A'*dv;
  end
  
  du = (w2 - sig2.*dx)./sig1;
  
  dlamu1 = (lamu1./fu1).*(-dx+du) - lamu1 - (1/tau)*1./fu1;
  dlamu2 = (lamu2./fu2).*(dx+du) - lamu2 - 1/tau*1./fu2;
  
  % make sure that the step is feasible: keeps lamu1,lamu2 > 0, fu1,fu2 < 0
  indp = find(dlamu1 < 0);  indn = find(dlamu2 < 0);
  s = min([1; -lamu1(indp)./dlamu1(indp); -lamu2(indn)./dlamu2(indn)]);
  indp = find((dx-du) > 0);  indn = find((-dx-du) > 0);
  s = (0.99)*min([s; -fu1(indp)./(dx(indp)-du(indp)); -fu2(indn)./(-dx(indn)-du(indn))]);
  
  % backtracking line search
  suffdec = 0;
  backiter = 0;
  while (~suffdec)
    xp = x + s*dx;  up = u + s*du; 
    vp = v + s*dv;  Atvp = Atv + s*Atdv; 
    lamu1p = lamu1 + s*dlamu1;  lamu2p = lamu2 + s*dlamu2;
    fu1p = xp - up;  fu2p = -xp - up;  
    rdp = gradf0 + [lamu1p-lamu2p; -lamu1p-lamu2p] + [Atvp; zeros(N,1)];
    rcp = [-lamu1p.*fu1p; -lamu2p.*fu2p] - (1/tau);
    rpp = rpri + s*Adx;
    suffdec = (norm([rdp; rcp; rpp]) <= (1-alpha*s)*resnorm);
    s = beta*s;
    backiter = backiter + 1;
    if (backiter > 32)
      disp('Stuck backtracking, returning last iterate.  (See Section 4 of notes for more information.)')
      xp = x;
      return
    end
  end
  
  
  % next iteration
  x = xp;  u = up;
  v = vp;  Atv = Atvp; 
  lamu1 = lamu1p;  lamu2 = lamu2p;
  fu1 = fu1p;  fu2 = fu2p;
  
  % surrogate duality gap
  sdg = -(fu1'*lamu1 + fu2'*lamu2);
  tau = mu*2*N/sdg;
  rpri = rpp;
  rcent = [-lamu1.*fu1; -lamu2.*fu2] - (1/tau);
  rdual = gradf0 + [lamu1-lamu2; -lamu1-lamu2] + [Atv; zeros(N,1)];
  resnorm = norm([rdual; rcent; rpri]);
  
  done = (sdg < pdtol) | (pditer >= pdmaxiter);
  
  disp(sprintf('Iteration = %d, tau = %8.3e, Primal = %8.3e, PDGap = %8.3e, Dual res = %8.3e, Primal res = %8.3e',...
    pditer, tau, sum(u), sdg, norm(rdual), norm(rpri)));
  if (largescale)
    disp(sprintf('                  CG Res = %8.3e, CG Iter = %d', cgres, cgiter));
  else
    disp(sprintf('                  H11p condition number = %8.3e', hcond));
  end
  
end