

% �Ŵ��㷨ʾ������ %
% ����f(x)=10*sin(5x)+7*cos(4x) �ڷ�Χ x��[0,10]�ϵ����ֵ %
% �� x ��ֵ��һ��10λ�Ķ�ֵ��ʽ��ʾΪ��ֵ���⣬һ��10λ�Ķ�ֵ���ṩ�ķֱ�����ÿΪ (10-0)/(2^10-1)��0.01 
% �������� [0,10] ��ɢ��Ϊ��ֵ�� [0,1023], x=0+10*b/1023, ���� b �� [0,1023] �е�һ����ֵ���� %


%% ************�Ŵ��㷨������************ %%
%Ҫ�󾫶Ȳ�����0.01
clear
clc
popsize=20;                      %����Ⱥ���С,�൱��
chromlength=10;                  %�����ַ�������(���峤��)
pc=0.6;                          %���ý�����ʣ�ֻ���������С�ڸ�ֵʱ���Ż�������档
                                 %0.6�൱��60%
pm=0.001;                        %���ñ�����ʣ������0.001�൱��0.1%


pop=initpop(popsize,chromlength); %���������ʼȺ�壬pop�Ĵ�СΪ20x10��20��ʾ����Ⱥ��ĸ��������׹�ʽ�е�Qֵ��
                                  %��10��ʾ���������ö����ƶ�Ⱥ����б�ʾ������λ��
%% ************�����Ŵ��㷨�ĵ�����Ѱ������ֵ************ %%                                  
for i=1:20                               %20Ϊ�Ŵ����������Ŵ��㷨�ĵ���������Ϊ20�����൱�����׹�ʽ�е�nֵ
        [objvalue]=calobjvalue(pop);                  %����Ŀ�꺯��
        fitvalue=calfitvalue(objvalue);                 %����Ⱥ����ÿ���������Ӧ��

        [newpop]=selection(pop,fitvalue);                
%����
        [newpop1]=crossover(newpop,pc);              
%����
        [newpop2]=mutation(newpop1,pc);              
%����
        
        [objvalue]=calobjvalue(newpop2);                
%����Ŀ�꺯��
        fitvalue=calfitvalue(objvalue);                       
%����Ⱥ����ÿ���������Ӧ��
        
        [bestindividual,bestfit]=best(newpop2,fitvalue);     %���Ⱥ������Ӧֵ���ĸ��弰����Ӧֵ
        y(i)=bestfit;                                                              %���ص� y ������Ӧ��ֵ�����Ǻ���ֵ
        x(i)=decodechrom(bestindividual,1,chromlength)*10/1023;      %���Ա��������ʮ����
        pop=newpop2;
end
fplot('10*sin(5*x)+7*cos(4*x)',[0 10])
hold on
plot(x,y,'r*')                                          
hold on

[z index]=max(y);             %�������ֵ����λ��
x5=x(index)                     %�������ֵ��Ӧ��xֵ
ymax=z














