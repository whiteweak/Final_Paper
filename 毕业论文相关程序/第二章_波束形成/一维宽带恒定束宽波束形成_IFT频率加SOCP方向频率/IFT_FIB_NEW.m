        
    % ************************** 基于多维离散傅里叶变换的FIB **************************%  
 
 
        clc;
        clear all
%         close all

        %% ************************** 试验参数设置（与论文上参数设置一致）************************** %%
        array_number=21;           %设置阵元个数为21个
        fl=10000;                   %最低频率 
%         fl=6000;                   %最低频率 
        fh=20000;                  %最高频率
        fr=fl;                     %参考频率

        c=340;                     %声速
        d=c/fr/2;                  %阵元间距取最高频率对应波长的一半，避免出现重叠（满足采样定理）
        lamda=c/fr;                %参考频率对应的波长

        ang=-90:1:90;              %方位角范围为-90°- 90°
        
        

       %% **************************Dolph-Chebyshev加权设计期望波束响应 **************************%
        
        array_num=13;                          %设计期望波束响应采用的参考阵元数
        
        vec=-(array_num-1)/2:(array_num-1)/2;                        %期望波束响应时的导向矢量是以中间阵元为参考点来计算的，因为共有11个，故vec取为-5:5
        
       %% *************** 指定旁瓣级水平值或者主瓣宽度值 *************** %%
        
        % *************** 指定旁瓣级 *************** %
        
        sidelobe=-30;                         %旁瓣级水平设置为-30dB
        R=10.^(-sidelobe/10);            
        z=cosh(1/(array_num-1)*acosh(sqrt(R)));     %计算论文中式2-4中定义的z值

        % *************** 指定主瓣宽度 *************** %
%         mainbeam_wide=24;                     %主瓣宽度
%         sin_theta_NN=sin(mainbeam_wide/2*pi/180);       
%         z=cos(pi/(2*(array_num-1)))/cos(pi*d/lamda*sin_theta_NN);

        % *************** 计算切比雪夫加权系数 *************** %
        
        wh=zeros(array_num,1);              %初始化切比雪夫加权值矩阵为0矩阵

        wh(1)=(z^(array_num-1))/2;          %计算当m=1时的权值

        for m=2:floor(array_num/2+1)        %floor函数的作用是取该数组某一个元素的小于等于它的最大整数，这儿指的是2<=m<(M/2+1)的情况
                                            %temp1作为一个临时变量，用于存放计算出来的权值；
              temp=0;                       %temp为临时变量
              for i=1:(m-1)
                   temp=temp+0.5*(array_num-1)*factorial_0(m-2)*factorial_0(array_num-i-1)*(z^(array_num-2*m+1))*((z^2-1)^(m-i))/...
                      factorial_0(m-i)/factorial_0(i-1)/factorial_0(m-i-1)/factorial_0(array_num-m);  
              end                           %factorial_0的作用是求一个数的阶乘，这一句实现的是式2-4中2<=m<(M/2+1)的情况下权值的取值。
              
              wh(m)=temp;                   %给权值矩阵赋值
        end

        for m=ceil(array_num/2+1):array_num       
              wh(m)=wh(array_num+1-m);      %ceil函数取的是大于等于该数的最小整数，这一句实现的是m>=M/2+1情况下的权值的计算   
        end


        % *************** 计算切比雪夫加权矢量 *************** %
        
        beam_point=0;            %期望方位角设为0°

        v0=(exp(-1j*2*pi*fr*vec*d*sin(beam_point*pi/180)/c)).'/array_num;     %常规波束形成的加权矢量

        wd=wh.*v0;               %计算切比雪夫加权后的导向矢量

        % *************** 期望波束响应 *************** %
        
        hh=wd.';                              %计算加权导向矢量的非共轭转置，方便之后波束响应的计算

        beam_temp=zeros(1,length(ang));       %对加权后波束响应矩阵初始化

        for i=1:length(ang)                   %计算切比雪夫加权后的针对每一个方位角上的波束响应
                beam_temp(i)=exp(-1j*vec*pi*sin(ang(i)*pi/180))*(hh');
        end
        
        beam=(beam_temp)/max(abs(beam_temp));          %对波束响应进行归一化
        beam_temp=abs(beam);                           %此时的beam_temp为将波束响应进行归一化并取绝对值之后的波束响应矩阵

        % *************** 确认主瓣宽度 *************** %
        IndMin=find(diff(sign(diff(beam_temp)))>0)+1;  %diff函数的作用是对一个数组进行微分（差分）运算；
                                                       %sign函数的作用是判断一个值的正负情况：正数返回值为1，负数返回值为-1,0则返回0；
                                                       %find函数的作用是用于返回满足括号中条件的元素的所在位置(位置的判定：在矩阵中，第一列开始，
                                                       %自上而下，依次为1，2，3...,然后再从第二列，第三列依次往后数)
                                                       %此处整个语句的作用是：寻找波束响应的凹陷点即波束图中谷底的位置对应的角度下标；
                                                       
        temp=ang(IndMin);                              %这里的temp指代的是波束图凹陷点所对应的角度值；
        temp=temp-beam_point*ones(1,length(temp));     %temp与temp代表的都是波束凹陷点的角度值；   
        
        temp=sort(abs(temp));                          %对波束凹陷点角度值的绝对值进行从小到大的排序  
        
        mainbeam1=temp(1);                             
        mainbeam2=temp(2);                             %mainbeam1和mainbeam2代表的是波束主瓣宽度值（角度）的一半
        
        mainbeam_wide=mainbeam1+mainbeam2;             %mainbeam_wide表示的是波束主瓣宽度的值
        
        % *************** 确认旁瓣级 *************** %
        
        temp=find(ang<(beam_point-mainbeam_wide/2));   %这里的temp代表的是旁瓣区域的角度值所对应的波束响应中的下标值
        sidelobe_level=20*log10(max(beam_temp(temp))); %这里的旁瓣水平指的是旁瓣区域的最大值的dB值，用以确定旁瓣值是否为设定的-30dB
        
        
     

        %% *************** 画出用切比雪夫加权设计出的期望波束响应的波束图 *************** %%
        
        figure(1);hold on;box on;
        plot(ang,20*log10(beam_temp),'k-');
                
%         plot(beam_point*ones(1,81),-80:0,'k');
%         plot(ang,sidelobe_level*ones(size(ang)),'k--');

        axis([min(ang),max(ang),-80,0]);
        
        h = gca;
        set(h,'FontSize',10,'FontName','Times New Roman');
        set(h,'FontName','Times New Roman');
        
        xlabel('Azimuth Angle(Degrees)');
        ylabel('Beam Pattern Gain(dB)');
        title('Desired Beam Pattern');

%         legend('0度','20度','40度','60度');
%         legend('最低数字频率pi/2','最低数字频率pi/3');
%         legend('参考阵元数9','参考阵元数11');
%         legend('SR 方法','SOCP 方法','IFT 方法');
%         legend('SOCP 方法','IFT 方法');
%         legend('情形1','情形2');

        
        %% *************** 频率不变波束形成 *************** %%
        Omega_range=2*pi*linspace(fl,fh,25);                    %子带划分,得到w的范围为15个窄带子带
%         Omega_range=2*pi*linspace(fl,fh,35); 

%         Omega_range=2*pi*[linspace(fl,6000,4),linspace(6000,fh,12)];     %子带划分

        M=64;                                           %w1的离散点数，即为做IDFT的点数
        Omega1_range=-pi*ones(1,M)+2*pi/M*(0:(M-1));    %w1的范围为-pi到（M-1/M）*pi，以pi/M为间隔，共有64个值；

        % *************** 波束响应赋值 *************** %
        P=zeros(length(Omega1_range),length(Omega_range));% P存放的是线阵的波束响应，维度由w1和w的大小决定

        for i=1:length(Omega1_range)      %此循环的作用是通过变量替换求解出了每一个频率分量即每一个子带
                                          %所对应的波束响应，同时设置A（w1）为0以让形成的波束图更加光滑
                                          %对应于论文上公式3-29；
                for k=1:length(Omega_range)      %循环范围为1到w的长度
                        if abs(Omega1_range(i))<(Omega_range(k)*d/c)     
                                temp=(Omega1_range(i)*c/Omega_range(k)/d);  
                                P(i,k)=exp(-1j*vec*pi*temp)*(hh');
                                %当满足条件时，用(c*w1)/(w*d)代入波束响应方程即完成了变量替换,进而得到加权波束响应
                        else
                                 P(i,k)=0;     %其余情况的A（Omega1）设为0，以使得P尽量的平滑；
                        end
                end
        end


        % *************** 对波束响应做离散傅里叶逆变换得到加权系数 *************** %
        B=zeros(M,M);   %M为w1的离散点数，为做IDFT的点数，与公式中N~相对应

        for i=1:M
            B(i,:)=exp(-1j*(0:(M-1))*Omega1_range(i));    
        end

        % temp=abs(sum(diag(B)))/M;
        % B=B+temp/1000*eye(M,M);

        D=B\P;   %用B矩阵去左除P矩阵，得到的结果是inv（B）*P即P矩阵乘以B矩阵的逆

        w=[D(((M-(array_number-3)/2):M),:);D((1:((array_number-1)/2+1)),:)];

        %-----------------波束形成-----------------------------------%
        beam_FIB=zeros(length(Omega_range),length(ang));

        vec1=(-(array_number-1)/2):(array_number-1)/2;

        for i=1:length(Omega_range)
                for k=1:length(ang)                              
                        beam_FIB(i,k)=exp(-1j*vec1*Omega_range(i)*sin(ang(k)*pi/180)*d/c)*w(:,i);
                end
        end
        
        beam_FIB=beam_FIB/max(max(abs(beam_FIB)));
        beam_FIB_a=abs(beam_FIB);
         
        %---------------确定主瓣所在区域--------------------------%
        mainlobe=(beam_point-mainbeam1):1:(beam_point+mainbeam2); 
        ang_index=zeros(1,length(mainlobe));

        for i=1:length(mainlobe)
               ang_index(i)=find(ang==mainlobe(i));
        end
        %%-------------------画图---------------------------------------%%
        figure(2);grid on;
        [x,y]=meshgrid(ang,Omega_range);

        mesh(x,y/2/fh,20*log10(beam_FIB_a));

        axis([min(ang),max(ang),fl/fh*pi,pi,-80,0]);

        title('Frequency Invariant Beamforming Based On IFT');
        zlabel('Beam Pattern Gain(dB)');
        ylabel('Digital Frequency(pi)');    
        xlabel('Azimuth Angle(Degrees)');
        
        h = gca;
        set(h,'FontSize',10,'FontName','宋体');
        set(h,'FontName','Times New Roman');

        figure(3);hold on;box on;
        errors=zeros(1,length(Omega_range));  
        sl_FIB=-100;
        
        for i=1:length(Omega_range)
                %----------------确认旁瓣级----------------------------%
                temp=find(ang<(beam_point-mainbeam_wide/2-3));
                
                temp1=20*log10(max(beam_FIB_a(i,temp)));                  
                if temp1>sl_FIB
                      sl_FIB=temp1;  
                end
                
                %---------------求均方误差和---------------------------%
                errors(i)=sum(abs(beam_FIB(i,ang_index)-beam(1,ang_index)).^2);

                
        end
        
        %-----------------不同频率分量对应的主瓣均方根误差------------------------------%      
        plot(Omega_range/2/fh,10*log10(errors)-10*log10(length(mainlobe))*ones(size(errors)),'k--o');%+,*,。,x,s（方块),d（菱形）,p（五角形）
        
        h = gca;
        set(h,'FontSize',10,'FontName','宋体');
        set(h,'FontName','Times New Roman');
        
        title('Mean Square Errors of Mainlobe Associated With Different Frequencies');
        ylabel('Mean Square Errors of Mainlobe(dB)');
        xlabel('Digital Frequency(pi)');
        axis([1.4,3.2,-90,-60]);
%         axis([0.6,3.2,-90,-50]);
%         legend('SR 方法','SOCP 方法','IFT 方法');
%         legend('情形1','情形2');
%         legend('0度','20度','40度');
%         legend('最低数字频率pi/2','最低数字频率pi/3');
%         legend('参考阵元数9','参考阵元数11');
%         legend('SOCP 方法','IFT 方法');
     

