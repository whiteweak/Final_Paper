function pop2=cross(pop1,pc)
    %交叉函数，pop1为经选择后的种群，pc为交叉概率
    %pop2为交叉后的新种群
    
    pop2=zeros(50,3);
    pnum=size(pop1,1);  %得到种群大小
    N=pc*pnum;          %N=40
    N=(mod(N,2)==0)*N+(mod(N,2)==1)*(N+1);  %保证交配个体数量为偶数
    
    r=rand(1,pnum);       %产生随机数选择交叉的染色体数目
    [~,index]=sort(r);
    pop1=pop1(index,:);
    
    for k=1:2:N
        a=rand;b=rand;c=rand;
        %-----染色体a0的交叉-------%
        pop2(k,1)=a*pop1(k,1)+(1-a)*pop1(k+1,1);
        pop2(k+1,1)=(1-a)*pop1(k,1)+a*pop1(k+1,1);
        %-----染色体b1的交叉-------%
        pop2(k,2)=b*pop1(k,2)+(1-b)*pop1(k+1,2);
        pop2(k+1,2)=(1-b)*pop1(k,2)+b*pop1(k+1,2);
        %-----染色体b2的交叉-------%
        pop2(k,3)=c*pop1(k,3)+(1-c)*pop1(k+1,3);
        pop2(k+1,3)=(1-c)*pop1(k,3)+c*pop1(k+1,3);
    end
    pop2(N+1:end,:)=pop1(N+1:end,:);
    
end