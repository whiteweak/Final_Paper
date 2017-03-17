

% 遗传算法示例程序 %
% 求函数f(x)=10*sin(5x)+7*cos(4x) 在范围 x∈[0,10]上的最大值 %
% 将 x 的值用一个10位的二值形式表示为二值问题，一个10位的二值数提供的分辨率是每为 (10-0)/(2^10-1)≈0.01 
% 将变量域 [0,10] 离散化为二值域 [0,1023], x=0+10*b/1023, 其中 b 是 [0,1023] 中的一个二值数。 %


%% ************遗传算法主程序************ %%
%要求精度不大于0.01
clear
clc
popsize=20;                      %设置群体大小,相当于
chromlength=10;                  %设置字符串长度(个体长度)
pc=0.6;                          %设置交叉概率，只有在随机数小于该值时，才会产生交叉。
                                 %0.6相当于60%
pm=0.001;                        %设置变异概率，这里的0.001相当于0.1%


pop=initpop(popsize,chromlength); %随机产生初始群体，pop的大小为20x10，20表示的是群体的个数即文献公式中的Q值，
                                  %而10表示的则是利用二进制对群体进行表示产生的位数
%% ************进行遗传算法的迭代，寻找最优值************ %%                                  
for i=1:20                               %20为遗传代数，即遗传算法的迭代次数设为20代，相当于文献公式中的n值
        [objvalue]=calobjvalue(pop);                  %计算目标函数
        fitvalue=calfitvalue(objvalue);                 %计算群体中每个个体的适应度

        [newpop]=selection(pop,fitvalue);                
%复制
        [newpop1]=crossover(newpop,pc);              
%交叉
        [newpop2]=mutation(newpop1,pc);              
%变异
        
        [objvalue]=calobjvalue(newpop2);                
%计算目标函数
        fitvalue=calfitvalue(objvalue);                       
%计算群体中每个个体的适应度
        
        [bestindividual,bestfit]=best(newpop2,fitvalue);     %求出群体中适应值最大的个体及其适应值
        y(i)=bestfit;                                                              %返回的 y 是自适应度值，而非函数值
        x(i)=decodechrom(bestindividual,1,chromlength)*10/1023;      %将自变量解码成十进制
        pop=newpop2;
end
fplot('10*sin(5*x)+7*cos(4*x)',[0 10])
hold on
plot(x,y,'r*')                                          
hold on

[z index]=max(y);             %计算最大值及其位置
x5=x(index)                     %计算最大值对应的x值
ymax=z














