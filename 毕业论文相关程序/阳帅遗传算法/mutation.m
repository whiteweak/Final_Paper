function pop3=mutation(pop2,pm)
    %变异函数，pop2为交叉后的种群，pm为变异概率
    %pop3为变异后的种群
    
%     pop3=zeros(popsize,3);
    pnum=size(pop2,1);  %得到种群大小
    N=pm*pnum*size(pop2,2);
%     N=pm*pnum;
    N=round(N);         %产生变异的个体数量（取整数）
    
    pop3=pop2;
    gene=pop3(:);       %得到基因序列
    r=rand(1,size(pop2,2)*pnum);    %通过随机数选择变异的基因
    [~,index]=sort(r);
    
    gene(index(1:N))=rand(1,N);     %得到变异后的基因
    pop3=reshape(gene,size(pop2));  %得到变异后的新种群
    
%     pop3=pop2(index,:);

%     for k=1:N
%         a=unidrnd(3);   %产生1-3之间的整数，选择需要变异的基因
%         pop3(k,a)=rand;
%         pop3(k,:)=rand(1,3);
%     end
    
end