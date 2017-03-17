%%*******************IIR-GA-Main***********************%%
%   ��Ϊ��GA�㷨ʵ��P113ϵͳ��ʶ������
%   ������Χa->(0,2)  b1->(-2,2)  b2->(-1,1)
%   select->ѡ����  cross->���溯��   mutation->���캯��
%   trans->����ת������(����Ҫ)
%   ����---����޸�2015-11-29

%% 1.�����ܱ���ͼ��
clear all;close all;clc
tic
%%%%---------------------�����ܱ��漰���Ȩ----------------------%%%
syms z;
H0=1/((1-1.2*z^(-1)+0.6*z^(-2))*(1-1.2*z+0.6*z^2))*(1/z);%%%��phidd(0)
phidd0=RESIDUE(H0);     %3.5714
[b1,b2]=meshgrid(-2:0.02:2,-1:0.01:1);

%%----�����ܺ�����Χ�߻��ֲ��ֳַ�����������������Χ�߻���----
rt1=0.5*(b1+sqrt(b1.^2+4*b2));%%rt1��rt2Ϊ�������rt3��rt4Ϊǰ����������
rt2=0.5*(b1-sqrt(b1.^2+4*b2));
% 
rt1=rt1.*(abs(rt1)<1)+100*(abs(rt1)>=1);
rt2=rt2.*(abs(rt2)<1)+100*(abs(rt2)>=1);

% rt3=1./rt1;rt4=1./rt2,������δ֪ϵͳ���ȶ�ϵͳ�����伫�㼴Ϊrt1��rt2��
%�ڵ�λԲ�ڣ��������ȶ�������rt3��rt4���ڵ�λԲ�⣬�������������㣻
res1=rt1./((1-rt1.*rt1).*(1-rt1.*rt2).*(rt1-rt2));%%����rt1��������
res2=rt2./((1-rt2.*rt1).*(1-rt2.*rt2).*(rt2-rt1));%%����rt2��������
%PART2---��ʽΪ-2*phidx(z)*H(z)/z
% H2=-2*z/((1-1.2*z+0.6*z^2).*(z^2-b1*z-b2));�����ĸ������ֱ�Ϊrt1,rt2,
%1+0.8165i,1-0.8156i,�����������ڵ�λԲ�⣬����
res3=-2*rt1./((1-1.2*rt1+0.6*rt1.^2).*(rt1-rt2));
res4=-2*rt2./((1-1.2*rt2+0.6*rt2.^2).*(rt2-rt1));
ksi=phidd0+res1+res2+res3+res4;     %%%%���ܱ��溯��
v=5*[0.1,0.3,0.5,0.7,0.9,0.99];%%ksiƫ�����ԭ����ΪPHIxx(z)=1��
figure(1);contour(b1,b2,ksi,v)
axis([-2 2 -1 1]);
hold on
plot(1.2,-0.6,'ro');%�������Ȩ
title('Ȩֵ��������','fontsize',15)


%% 2.GAѰ����ѽ�
popsize=50;     %��Ⱥ��ģ
pr=0.2;         %��Ⱥ��ѡ�����0.2  0.2
pc=0.6;         %��Ⱥ�Ľ������0.6  0.4
pm=0.1;         %��Ⱥ�ı������0.1  0.2
maxgen=100;     %��Ⱥ������������


%% 3.GAѰ�ţ��ɼ���LMS���ӣ�
pop0=rand(popsize,3);           %�����������(����)
fit=fitness(pop0);              %������Ⱥ����Ӧ��ֵ
wbest=zeros(maxgen,3);          %���Ա���ÿ�ε��������Ÿ��壨����Ȩ��
W=zeros(popsize,3,maxgen);
for k=1:maxgen
    pop1=select(pop0,fit,pr);       %pop0����ѡ���ƺ�õ�pop1
    pop2=cross(pop1,pc);            %pop1���������õ�pop2
    pop3=mutation(pop2,pm);         %pop2���������õ�pop3
    popnew=pop3;                    %��ѡ�񡢽��桢�����õ�������Ⱥpopnew
    
    %---------�Ŵ��㷨�м���LMS��SHARF����------------%
%     flag='SHARF';
%     ps=0.5;K=10;%���ӵ�ִ�и��ʼ���������
%     popnew=myoperator(popnew,flag,ps,K); 
    %--------------LMS����ִ�н���-END------------%
    
    pop0=popnew;
    fit=fitness(pop0); 
    %------------�ҵ�ÿ�ε���������Ȩ������----------%
    ind=find(fit==min(fit),1);
    wbest(k,:)=pop0(ind,:);   
    W(:,:,k)=[2*pop0(:,1),4*pop0(:,2)-2,2*pop0(:,3)-1];
end
toc

%------------�õ����һ�ε��������Ⱥ-------------------%
[~,index]=sort(fit);
pop=pop0(index,:);
a0c=pop(:,1);b1c=pop(:,2);b2c=pop(:,3);

a0=a0c*(2-0)+0;b1=b1c*(2+2)-2;b2=b2c*(1+1)-1;     %�����õ�a0  b1  b2����ʵֵ
w=[a0,b1,b2];
disp(w(1:5,:))           %��ʾ���һ�ε������ǰ5�����Ÿ���

%-----------------�õ�ÿ�ε���������Ȩ(�����)------------------%
wb(:,1)=wbest(:,1)*(2-0)+0;
wb(:,2)=wbest(:,2)*(2+2)-2;
wb(:,3)=wbest(:,3)*(1+1)-1;


%% ������ʾGAѰ�Ź���
% ---------����ÿ����Ⱥ�ı仯�ֲ�----------%
fig=figure(1);
h0=plot(wb(:,2),wb(:,3),'b*');
xlabel('b1','fontsize',15);ylabel('b2','fontsize',15);
% title('GA-LMSѰ��')
title('GAѰ��','fontsize',15)
legend('�ȸ���','���Ȩ','GA��Ⱥ')
% figure(1);k=10; title(['GA-LMSѰ��(k=' num2str(k) ')']);
% set(h0,'xdata',W(:,2,k));
%     set(h0,'ydata',W(:,3,k));

% for k=[1,10:10:maxgen]
%     set(h0,'xdata',W(:,2,k));
%     set(h0,'ydata',W(:,3,k));
%     title(['GA-LMSѰ��(k=' num2str(k) ')'],'fontsize',15);
%     pause(1);
% %     drawnow;
% end

aviobj = avifile('GA.avi','compression','None');
 for k=[(1:9),10:5:maxgen]
    set(h0,'xdata',W(:,2,k));
    set(h0,'ydata',W(:,3,k));
%     title(['GA-LMSѰ��(k=' num2str(k) ')']);
title(['GAѰ��(k=' num2str(k) ')'],'fontsize',15);
%     pause(1);
%     drawnow;
    
    F=getframe(fig);%��ȡ��ǰ����
    aviobj=addframe(aviobj,F);%����avi������
    im=frame2im(F);%ת��gifͼƬ��ֻ����256ɫ
    [I map]=rgb2ind(im,256);
    %=-д��GIF89a��ʽ�ļ�%----
    if k==1
        imwrite(I,map,'GA.gif','GIF', 'Loopcount',inf,'DelayTime',1);
    else
         imwrite(I,map,'GA.gif','GIF','WriteMode','append','DelayTime',1);
    end
 end
close(fig);
%�ر�avi����
aciobj=close(aviobj);



%% 4.ѧϰ���ߡ�Ȩֵ��������
ksi=myfunksi(wb);
figure(2)
plot(ksi)
xlabel('k','fontsize',15);ylabel('{\xi}','fontsize',15);title('ѧϰ����{\xi}-k','fontsize',15)
%---Ȩֵ��һ��������---%
figure(3)
subplot(311);plot(wb(:,1));xlabel('k','fontsize',15);ylabel('a0','fontsize',15);
title('Ȩֵ����ֵ','fontsize',15);hold on;plot(ones(1,maxgen),'r')
subplot(312);plot(wb(:,2));xlabel('k','fontsize',15);ylabel('b1','fontsize',15);
hold on;plot(1.2*ones(1,maxgen),'r')
subplot(313);plot(wb(:,3));xlabel('k','fontsize',15);ylabel('b2','fontsize',15);
hold on;plot(-0.6*ones(1,maxgen),'r')

% figure(3)
% k=1:maxgen;
% plot(k,wb(:,1),'r',k,1*ones(1,maxgen),'k',k,wb(:,2),'r',k,1.2*ones(1,maxgen),'k',...
%     k,wb(:,3),'r',k,-0.6*ones(1,maxgen),'k')
% set(gca,'ytick',[-0.6 0,1,1.2])
% xlabel('k','fontsize',15);ylabel('w','fontsize',15);title('Ȩֵ��������','fontsize',15)
% legend('Ȩ��������w(a0,b1,b2)','���Ȩֵwbest(a0,b1,b2)','location','best')



