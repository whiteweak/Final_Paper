function W = myfunIIR2(flag,x,K,M,w0)
%UNTITLED5    基于IIR结构的LMS算法和SER算法   课本P113
%   flag=LMS,采用LMS算法，flag=SER,采用SER算法
%   x为输入序列，K为迭代次数，M为速率因子，W为迭代权
%   最后更新日期：2015-11-3
%   w0为初始权

W=zeros(K+1,3); 
W(3,:)=w0;
alpha0=zeros(1,K);beta1=zeros(1,K);beta2=zeros(1,K);
y=zeros(1,K);d=zeros(1,K);
switch flag
    case 'LMS'
        for k=3:K
            U=[x(k) y(k-1) y(k-2)];         %%%第k次迭代输入
            y(k)=W(k,:)*U';                 %%%第k次迭代输出
            d(k)=x(k)+1.2*y(k-1)-0.6*y(k-2);     %%k时刻未知系统的输出（模型的期望）
            alpha0(k)=x(k)+W(k,2)*alpha0(k-1)+W(k,3)*alpha0(k-2);
            beta1(k)=y(k-1)+W(k,2)*beta1(k-1)+W(k,3)*beta1(k-2);
            beta2(k)=y(k-2)+W(k,2)*beta2(k-1)+W(k,3)*beta2(k-2);
            delta=-2*(d(k)-y(k))*[alpha0(k),beta1(k),beta2(k)]';%%%第k次迭代梯度
            W(k+1,:)=W(k,:)-(M*delta)';    
        end
    case 'SER'
        %%---------求特征值均值--------%%
%         aveLambda=1/800*(x*x');
        dk(1)=x(1);dk(2)=x(2);
        for l=3:K   %迭代K次计算出期望值
             dk(l)=x(l)+1.2*dk(l-1)-0.6*dk(l-2);
        end
        Exx0=1/K*x(1:K)*x(1:K)';Exy1=1/K*x(1:K)*[0 dk(1:K-1)]';Exy2=1/K*x(1:K)*[0 0 dk(1:K-2)]';
        Eyy1=1/K*[0 dk(1:K-1)]*[0;0;dk(1:K-2)'];Eyy0=1/K*dk*dk';
        EYY1=1/K*[0 dk(1:K-1)]*[0 dk(1:K-1)]';%E[y(k-1)^2]
        EYY2=1/K*[0 0 dk(1:K-2)]*[0 0 dk(1:K-2)]';%E[y(k-2)^2]
        R=[Exx0 Exy1 Exy2;Exy1 EYY1 Eyy1;Exy2 Eyy1 EYY2];
        lambda=eig(R);
        aveLambda=mean(lambda);
%         disp(['特征值平均值：' num2str(aveLambda)])
        %%---------求特征值均值-end-------%%
        %%---------计算前2次迭代--------%%
        %aveLambda=1.38;
        alpha=0.93;q0=1;invQ(:,:,1)=diag([q0 q0 q0]);
        %----------k=1时刻的迭代----------%
        alpha0(1)=x(1);beta1(1)=0;beta2(1)=0;
        y(1)=0;d(1)=x(1);
        delta=-2*(d(1)-y(1))*[alpha0(1),beta1(1),beta2(1)]';
        S=invQ(:,:,1)*[x(2);y(1);0];                         %%k=2
        gamma=alpha+[x(2);y(1);0]'*S;
        invQ(:,:,2)=1/alpha*(invQ(:,:,1)-1/gamma*(S*S'));
        W(2,:)=W(1,:)-aveLambda*(1-alpha^2)/(1-alpha)*(M*invQ(:,:,1)*delta)';  
        %----------k=2时刻的迭代----------%
        alpha0(2)=x(2)+W(2,2)*x(1);beta1(2)=y(1);beta2(2)=0;
        y(2)=W(2,:)*[x(2),y(1),0]';
        d(2)=x(2)+1.2*y(1);
        delta=-2*(d(2)-y(2))*[alpha0(2),beta1(2),beta2(2)]';
        S=invQ(:,:,2)*[x(3);y(2);y(1)];
        gamma=alpha+[x(3);y(2);y(1)]'*S;
        invQ(:,:,3)=1/alpha*(invQ(:,:,2)-1/gamma*(S*S'));
        W(3,:)=W(2,:)-aveLambda*(1-alpha^3)/(1-alpha)*(M*invQ(:,:,2)*delta)';
        %----------k=3-K时刻的迭代----------%
        for k=3:K
            U=[x(k) y(k-1) y(k-2)];
            S=invQ(:,:,k-1)*U';
            gamma=alpha+U*S;
            invQ(:,:,k)=1/alpha*(invQ(:,:,k-1)-1/gamma*(S*S'));
            y(k)=W(k,:)*U';         %%%第k次迭代输出
            d(k)=x(k)+1.2*y(k-1)-0.6*y(k-2);     %%k时刻未知系统的输出（模型的期望）
            alpha0(k)=x(k)+W(k,2)*alpha0(k-1)+W(k,3)*alpha0(k-2);
            beta1(k)=y(k-1)+W(k,2)*beta1(k-1)+W(k,3)*beta1(k-2);
            beta2(k)=y(k-2)+W(k,2)*beta2(k-1)+W(k,3)*beta2(k-2);
            delta=-2*(d(k)-y(k))*[alpha0(k),beta1(k),beta2(k)]';
            W(k+1,:)=W(k,:)-aveLambda*(1-alpha^(k+1))/(1-alpha)*(M*invQ(:,:,k)*delta)';
        end     
    case 'SHARF'
        er=zeros(1,K);
        v=zeros(1,K);c=0.5;
        for k=3:K
            U=[x(k) y(k-1) y(k-2)];
            y(k)=W(k,:)*U'; 
            d(k)=x(k)+1.2*y(k-1)-0.6*y(k-2);     %%k时刻未知系统的输出（模型的期望）
            er(k)=d(k)-y(k);
            v(k)=er(k)+er(1:k)*(power(c,k:-1:1))';
            delta=-2*v(k)*[x(k) y(k-1) y(k-2)]';
            W(k+1,:)=W(k,:)-(M*delta)';
        end
    otherwise
        disp('请按照格式输入');
end
% W=W(1:K,:);
W=W(K,:);

end

