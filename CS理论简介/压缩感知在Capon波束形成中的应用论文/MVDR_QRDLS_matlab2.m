% ���з���ͼ��ʹ��MVDR���������Ȩ������
% ��Ԫ��Ϊ16�����Ϊ�벨�����ź�Ƶ��Ϊ1.3G��ָ���Ϊ0�ȣ����������ڣ�40��20��45
% ���ǲ���ָ�������ź����
% 2007.5.22���ź�

clear all
clc
f0 = 1.3*10^9;                    % �ź�Ƶ��
fs = 3*10^9;                   % ����Ƶ��
snr = [40;40];                  % ���ŵ������
w0 = 20/180*pi;                  % Ŀ��ָ���
M = 9;                         % ��Ԫ��ΪM
K = 2;                          % �����ź���Ŀ
seta =[20/180*pi,45/180*pi];   % ���ŷ���          
N = 100;                        % ��������
c = 3*10^8;                     % ����
d = 0.5*c/f0;                   % ��Ԫ���
x = zeros(M,1);                 % ����ʸ��
R = zeros(M,M);                 % ��������Э�������
tic
% �ź�Դ
nn = 1:N;%������Ƶ��Ǹ����ź�,nnΪ��N��������
for k=1:K
    for n=1:N
        fai(1,n) = rand;
    end
    s(k,:) = exp(i*2*pi*(f0*nn/fs+fai*k));%�����ź�
end

% ���췽��ʸ��
theta = -90*pi/180:pi/180:90*pi/180;
WW = length(theta);
P = zeros(1,WW);
for m=1:M
        a(m,:) = exp(-i*pi*(m-1)*sin(theta)); % �������������ڷ�������
end
 
% ��������Э�������R
X=[];
for n=1:N                  
    for m=1:M
        for k=1:K
             A(m,k) = exp(-i*pi*(m-1)*sin(seta(k)));
        end
    end
    noise(:,n)=randn(M,1);
    y = A*(10.^(snr/20).*s(:,n))+noise(:,n);
    R = R+y*y';
    
    X=[X,y];
end
R = R/N;
Ri = inv(R);

% �����Ȩ����
for k=1:M
    a0(k,1) =exp(-i*pi*(k-1)*sin(w0));  % ָ������ 
end 
W = Ri*a0/(a0'*Ri*a0);
f = abs(W'*a).^2;
figure;
f = f./max(f);
%load QRD_LS_f
plot(theta*180/pi,10*log10(f),':r');
%theta*180/pi,10*log10(f_LS),'b-');%
% title('power pattern of linear phased array');
xlabel('\theta/deg'); 
ylabel('�������棨dB��');
grid on
toc
