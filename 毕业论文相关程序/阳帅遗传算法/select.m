function pop1=select(pop0,fit,pr)
    %��Ⱥѡ������pop0Ϊѡ��ǰ��Ⱥ��pop1Ϊѡ������Ⱥ%
    %fitΪ��Ⱥ����Ӧ��ֵ��prΪ��Ⱥ�ĸ��Ƹ���%
    [~,index]=sort(fit);
    pnum=size(pop0,1);
    N=round(pnum*pr);
    pop1(1:N,:)=pop0(index(1:N),:);%�˴��豣֤N*prΪ����
    pop1(N+1:2*N,:)=pop0(index(1:N),:);
    pop1(2*N+1:pnum,:)=pop0(index(N+1:pnum-N),:);
end