function popnew=myoperator(pop3,flag,ps,K)
    %   pop3Ϊ��������Ⱥ
    %   flagΪ��־λ������ѡ�����ӣ���ȡLMS��SER��SHARF��������
    %   popnewΪ������Ѱ�ź����Ⱥ
    %   psΪ����ִ�еĸ���
    %   KΪ���ӵĵ�������
    
    pnum=size(pop3,1);
    N=ps*pnum;
    N=round(N);          %��Ҫִ�����ӵĸ�����Ŀ��ȡ������
    
    %--------���������ȷ��ִ�����ӵĸ���--------%
    r=rand(1,pnum);
    [~,index]=sort(r);
    popnew=pop3(index,:);
    
%     K=10;       %���ӵĵ�������
%     flag='LMS';
    x=randn(1,200);
    M=diag([.05 .005 .0025]);
    for k=1:N
        %------------�˴���������------------%
        w0=[2*popnew(k,1),4*popnew(k,2)-2,2*popnew(k,3)-1];
%         w1=myLMS(flag,x,K,M,w0);
        w1=myfunIIR2(flag,x,K,M,w0);
        %----------������ִ�к�õ��Ĳ����ٴα���----------%
        popnew(k,:)=[w1(1)/2,(w1(2)+2)/4,(w1(3)+1)/2];
    end
   
end