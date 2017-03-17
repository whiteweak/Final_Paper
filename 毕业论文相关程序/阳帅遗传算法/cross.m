function pop2=cross(pop1,pc)
    %���溯����pop1Ϊ��ѡ������Ⱥ��pcΪ�������
    %pop2Ϊ����������Ⱥ
    
    pop2=zeros(50,3);
    pnum=size(pop1,1);  %�õ���Ⱥ��С
    N=pc*pnum;          %N=40
    N=(mod(N,2)==0)*N+(mod(N,2)==1)*(N+1);  %��֤�����������Ϊż��
    
    r=rand(1,pnum);       %���������ѡ�񽻲��Ⱦɫ����Ŀ
    [~,index]=sort(r);
    pop1=pop1(index,:);
    
    for k=1:2:N
        a=rand;b=rand;c=rand;
        %-----Ⱦɫ��a0�Ľ���-------%
        pop2(k,1)=a*pop1(k,1)+(1-a)*pop1(k+1,1);
        pop2(k+1,1)=(1-a)*pop1(k,1)+a*pop1(k+1,1);
        %-----Ⱦɫ��b1�Ľ���-------%
        pop2(k,2)=b*pop1(k,2)+(1-b)*pop1(k+1,2);
        pop2(k+1,2)=(1-b)*pop1(k,2)+b*pop1(k+1,2);
        %-----Ⱦɫ��b2�Ľ���-------%
        pop2(k,3)=c*pop1(k,3)+(1-c)*pop1(k+1,3);
        pop2(k+1,3)=(1-c)*pop1(k,3)+c*pop1(k+1,3);
    end
    pop2(N+1:end,:)=pop1(N+1:end,:);
    
end