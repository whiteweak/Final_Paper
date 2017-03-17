function popnew=myoperator(pop3,flag,ps,K)
    %   pop3为变异后的种群
    %   flag为标志位，用以选择算子，可取LMS、SER、SHARF三种算子
    %   popnew为经算子寻优后的种群
    %   ps为算子执行的概率
    %   K为算子的迭代次数
    
    pnum=size(pop3,1);
    N=ps*pnum;
    N=round(N);          %需要执行算子的个体数目（取整数）
    
    %--------产生随机数确定执行算子的个体--------%
    r=rand(1,pnum);
    [~,index]=sort(r);
    popnew=pop3(index,:);
    
%     K=10;       %算子的迭代次数
%     flag='LMS';
    x=randn(1,200);
    M=diag([.05 .005 .0025]);
    for k=1:N
        %------------此处需解码参数------------%
        w0=[2*popnew(k,1),4*popnew(k,2)-2,2*popnew(k,3)-1];
%         w1=myLMS(flag,x,K,M,w0);
        w1=myfunIIR2(flag,x,K,M,w0);
        %----------将算子执行后得到的参数再次编码----------%
        popnew(k,:)=[w1(1)/2,(w1(2)+2)/4,(w1(3)+1)/2];
    end
   
end