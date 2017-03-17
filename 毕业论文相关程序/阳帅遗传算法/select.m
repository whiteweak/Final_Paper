function pop1=select(pop0,fit,pr)
    %种群选择函数，pop0为选择前种群，pop1为选择后的种群%
    %fit为种群的适应度值，pr为种群的复制概率%
    [~,index]=sort(fit);
    pnum=size(pop0,1);
    N=round(pnum*pr);
    pop1(1:N,:)=pop0(index(1:N),:);%此处需保证N*pr为整数
    pop1(N+1:2*N,:)=pop0(index(1:N),:);
    pop1(2*N+1:pnum,:)=pop0(index(N+1:pnum-N),:);
end