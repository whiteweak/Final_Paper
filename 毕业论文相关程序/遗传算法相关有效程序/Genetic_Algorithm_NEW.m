

% 遗传算法示例程序 %
% 求函数f(x)=10*sin(5x)+7*cos(4x) 在范围 x∈[0,10]上的最大值 %
% 将 x 的值用一个10位的二值形式表示为二值问题，一个10位的二值数提供的分辨率是每为 (10-0)/(2^10-1)≈0.01 
% 将变量域 [0,10] 离散化为二值域 [0,1023], x=0+10*b/1023, 其中 b 是 [0,1023] 中的一个二值数。 %


%% ************遗传算法主程序************ %%
% 要求精度不大于0.01
clear
clc
popsize=20;                      % 设置群体大小,
chromlength=10;                  % 设置字符串长度(个体长度)
pc=0.6;                          % 设置交叉概率，只有在随机数小于该值时，才会产生交叉。
                                 % 0.6相当于60%
pm=0.01;                        % 设置变异概率，这里的0.001相当于0.1%


pop=rand(popsize,chromlength);   % rand函数随机产生大小介于{0,1}、行数为popsize、列数为chromlength的随机数矩阵，
                                 % 这时的矩阵元素的大小均为0-1之间；
pop=round(pop);                  % round函数则对矩阵的每个元素进行圆整即将其替换为与其最近的整数,进而产生初始
                                 % 种群，初始种群的值全为0或者1；
                                 
%% ************进行遗传算法的迭代，寻找最优值************ %%                                  
for j=1:20                               % 20为遗传代数，即遗传算法的迭代次数设为20代，相当于文献公式中的n值
        
       %% ************计算目标函数值 ************ %%
        % 这一段程序的目的是实现目标函数的计算，将二进制值域中的数转化为变量域的数（x），
        % 其公式采用本文示例仿真，可根据不同优化问题予以修改
        
        spoint=1;         % 参数spoint表示待解码的二进制串的起始位置
        length=10;        % 参数1ength表示所截取的二进制数据段的长度及需要进行二-十进制转换的数据段
                          % 由于要求精度不大于0.01，因此最小整数为1023，即需要10位二进制
        pop1=pop(:,spoint:spoint+length-1);    % pop1等价于pop，这里的目的是为了不改变原始pop的值
        [px,py]=size(pop);                     %求pop矩阵的行和列数
        for i=1:py        % 将二进制转换为10进制的循环                      
                pop1(:,i)=2.^(py-i).*pop(:,i);
        end
        pop2=sum(pop1,2);      % 求pop1的每行之和
        temp1=pop2;            % temp1表示的就是将二进制数转换成对应的十进制数的结果
        
        x=temp1*10/1023;       % x表示的则是经过转换和计算之后的20个位于0-10之间的可用于进行计算的数值
        
        objvalue=10*sin(5*x)+7*cos(4*x);   % 利用给出的函数表达式计算目标函数值
        
       %% ************ 计算初始群体中每个个体的适应度 ************ %% 
       
        [px,py]=size(objvalue);      % 目标值有正有负
        for i=1:px                   % 利用一种取值规则计算个体的适应度（这里的适应度是以函数值的大小来定义的）
                if objvalue(i)>0                    
                        temp=objvalue(i);          
                else
                        temp=0.0;
                end
                fitvalue(i)=temp;
        end        
        fitvalue=fitvalue';          % 计算群体中每个个体的适应度

       %% ************ 采用轮盘赌选择法进行选择复制操作 ************ %%
       
        % 轮盘赌选择法根据方程 pi=fi/∑fi=fi/fsum来进行，其具体选择步骤为：
        %（fit表示的是适应度值，pi表示的则是其概率，fsum表示的是所有适应度的总和）
        % 1） 在第 t 代，由（1）式计算 fsum 和 pi 
        % 2） 产生 {0,1} 的随机数 rand( .)，求 s=rand( .)*fsum
        % 3） 求 所有fi≥s 中最小的 k ，则第 k 个个体被选中
        % 4） 进行 N 次2）、3）操作，得到 N 个个体，成为第 t=t+1 代种群
        
        totalfit=sum(fitvalue);       % 求所有适应度值之和
        fitvalue=fitvalue/totalfit;   % 单个个体被选择的概率（此概率即为每个个体的适应度与适应度总和的比值）
        fitvalue=cumsum(fitvalue);    % 如cumsum函数的作用是用于计算一个数组的各行或各列的累加值，如：
                                      % fitvalue=[1 2 3 4]，则 cumsum(fitvalue)=[1 3 6 10] 
        [px,py]=size(pop);            % 计算pop矩阵的大小
        ms=sort(rand(px,1));          % 利用rand函数随机产生从小到大排列的随机数组
        fitin=1;                      % 设置循环初始变量
        newin=1;                      % 设置循环初始变量
        while newin<=px               % 循环选出20个新个体
                if(ms(newin))<fitvalue(fitin)
                        newpop(newin,:)=pop(fitin,:);
                        newin=newin+1;
                else
                        fitin=fitin+1;
                end
        end                
        
        %% ************ 利用交叉算法对群体的每个个体进行基因分裂与重组 ************ %%
        
        % 交叉(crossover)，群体中的每个个体之间都以一定的概率pc进行交叉，即两个个体从各自字符串的某一位置
        % （一般是随机确定）开始互相交换，这类似生物进化过程中的基因分裂与重组。例如，假设2个父代个体x1，x2为：
        % x1=0100110
        % x2=1010001
        % 从每个个体的第3位开始交叉，交又后得到2个新的子代个体y1，y2分别为：
        % y1＝0100001
        % y2＝1010110
        % 这样2个子代个体就分别具有了2个父代个体的某些特征。利用交又我们有可能由父代个体在子代组合成具有更高适合度的个体。
        % 交叉是遗传算法区别于其它传统优化方法的主要特点之一。
        
        [px,py]=size(newpop);
        newpop1=ones(size(newpop));
        for i=1:2:px-1     %步长为2，是将相邻的两个个体进行交叉                                               
                if(rand<pc)
                        cpoint=round(rand*py);
                        newpop1(i,:)=[newpop(i,1:cpoint),newpop(i+1,cpoint+1:py)];
                        newpop1(i+1,:)=[newpop(i+1,1:cpoint),newpop(i,cpoint+1:py)];
                else
                        newpop1(i,:)=newpop(i,:);
                        newpop1(i+1,:)=newpop(i+1,:);
                end
        end  
        
        %% ************ 利用突变算法对群体的每个个体进行基因突变 ************ %%        
        % 变异(mutation)，基因的突变普遍存在于生物的进化过程中。变异是指父代中的每个个体的每一位都以概率pm翻转，
        %即由“1”变为“0”，或由“0”变为“1”。
        %遗传算法的变异特性可以使求解过程随机地搜索到解可能存在的整个空间，因此可以在一定程度上求得全局最优解。
        
        [px,py]=size(newpop1);
        newpop2=ones(size(newpop1));
        for i=1:px
                if(rand<pm)
                        mpoint=round(rand*py);     % 产生的变异点在1-10之间 
                        if mpoint<=0
                                mpoint=1;          % 记录变异位置               
                        end
                        newpop2(i,:)=newpop1(i,:);
                        if any(newpop2(i,mpoint))==0
                                newpop2(i,mpoint)=1;
                        else
                                newpop2(i,mpoint)=0;
                        end
                else
                newpop2(i,:)=newpop1(i,:);
                end
        end
                
         %% ************计算目标函数值 ************ %%
         
        spoint=1;         % 参数spoint表示待解码的二进制串的起始位置
        length=10;        % 参数1ength表示所截取的二进制数据段的长度及需要进行二-十进制转换的数据段
                          % 由于要求精度不大于0.01，因此最小整数为1023，即需要10位二进制
        newpop3=newpop2(:,spoint:spoint+length-1);    % newpop3等价于newpop2，这里的目的是为了不改变原始newpop2的值
        [px,py]=size(newpop2);                        %求newpop2矩阵的行和列数
        for i=1:py        % 将二进制转换为10进制的循环                      
                newpop3(:,i)=2.^(py-i).*newpop2(:,i);
        end
        newpop4=sum(newpop3,2);      % 求newpop3的每行之和
        temp1=newpop4;               % temp1表示的就是将二进制数转换成对应的十进制数的结果
        
        x=temp1*10/1023;       % x表示的则是经过转换和计算之后的20个位于0-10之间的可用于进行计算的数值
        
        objvalue=10*sin(5*x)+7*cos(4*x);   % 利用给出的函数表达式计算目标函数值
        
       %% ************ 计算群体中每个个体的适应度 ************ %% 
       
        [px,py]=size(objvalue);      % 目标值有正有负
        for i=1:px                   % 利用一种取值规则计算个体的适应度（这里的适应度是以函数值的大小来定义的）
                if objvalue(i)>0                    
                        temp=objvalue(i);          
                else
                        temp=0.0;
                end
                fitvalue(i)=temp;
        end        
        fitvalue=fitvalue';          % 计算群体中每个个体的适应度
        
        %% ************ 找出群体中适应度值最大的个体及其适应度值 ************ %%
        
        [px,py]=size(newpop2);
        bestindividual=newpop2(1,:);
        bestfit=fitvalue(1);
        for i=2:px
                if fitvalue(i)>bestfit
                        bestindividual=newpop2(i,:);
                        bestfit=fitvalue(i);
                end
        end
        %% ************ 解码自变量x并更新父群体 ************ %%
        
        y(j)=bestfit;       %返回的 y 是自适应度值，而非函数值

        spoint=1;         % 参数spoint表示待解码的二进制串的起始位置
        length=10;        % 参数1ength表示所截取的二进制数据段的长度及需要进行二-十进制转换的数据段
                          % 由于要求精度不大于0.01，因此最小整数为1023，即需要10位二进制
        bestindividual_1=bestindividual(:,spoint:spoint+length-1);    % pop1等价于pop，这里的目的是为了不改变原始pop的值
        [px,py]=size(bestindividual);                     %求pop矩阵的行和列数
        for i=1:py        % 将二进制转换为10进制的循环                      
                bestindividual_1(:,i)=2.^(py-i).*bestindividual(:,i);
        end
        bestindividual_2=sum(bestindividual_1,2);      % 求pop1的每行之和
        temp1=bestindividual_2;            % temp1表示的就是将二进制数转换成对应的十进制数的结果
        x(j)=temp1*10/1023;
        pop=newpop2;
end
fplot('10*sin(5*x)+7*cos(4*x)',[0 10])
hold on
plot(x,y,'r*')                                          
hold on

[z index]=max(y);             %计算最大值及其位置
x5=x(index)                     %计算最大值对应的x值
ymax=z














