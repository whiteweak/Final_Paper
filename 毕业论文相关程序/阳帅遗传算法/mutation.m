function pop3=mutation(pop2,pm)
    %���캯����pop2Ϊ��������Ⱥ��pmΪ�������
    %pop3Ϊ��������Ⱥ
    
%     pop3=zeros(popsize,3);
    pnum=size(pop2,1);  %�õ���Ⱥ��С
    N=pm*pnum*size(pop2,2);
%     N=pm*pnum;
    N=round(N);         %��������ĸ���������ȡ������
    
    pop3=pop2;
    gene=pop3(:);       %�õ���������
    r=rand(1,size(pop2,2)*pnum);    %ͨ�������ѡ�����Ļ���
    [~,index]=sort(r);
    
    gene(index(1:N))=rand(1,N);     %�õ������Ļ���
    pop3=reshape(gene,size(pop2));  %�õ�����������Ⱥ
    
%     pop3=pop2(index,:);

%     for k=1:N
%         a=unidrnd(3);   %����1-3֮���������ѡ����Ҫ����Ļ���
%         pop3(k,a)=rand;
%         pop3(k,:)=rand(1,3);
%     end
    
end