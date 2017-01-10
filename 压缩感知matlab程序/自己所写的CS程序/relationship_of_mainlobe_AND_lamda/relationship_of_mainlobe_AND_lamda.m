        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                 Dolph-Chebyshev加权方法                   %
        % copyright Liao Fengyi                                      %
        % last modified at 2015.09.22 （finished）                             %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        clc;
        clear all

        %% -----------------参数设置-----------------------------------%
        fh=16000;                       %最高频率
%         fl=4000;                        %最低频率
        f0=16000;                          %参考频率

        c=340;                          %声速
%         d=(c/fh)*2/3;                       %阵元间距为参考频率波长的一半
        lamda=c/fh;                     %参考频率对应的波长

        ang=-90:0.1:90;                 %方位角范围为-90°到90°


%         wide_c=[-20,-25,-30,-35,-40];         %旁瓣级取值
        wide_c=-30;        %旁瓣级取值

%         wide_c=[15,20,25,30,35,40,45,50,55];  %两零点间束宽取值   

%         wide_c=[0,5,10,15,20,25,30,35,40,45,50];        %波束指向取值

%         wide_c=[5,10,15,20,25];                       %阵元个数取值


        %% ----------------Dolph-Chebyshev加权方法设计窄带波束形成器-----------%%
        D=lamda/10:lamda/5:5*lamda/2;
%         D=8*lamda/5;
        ang_temp=zeros(1,length(D));
        sl_FIB=zeros(1,length(D));
        SL_wide=zeros(1,length(D));
        Ang_3dB_SOCP_A=zeros(1,length(D));
    for d_num=1:length(D)
        d=D(d_num);
        for coff_num=1:length(wide_c)   %系数个数设为coff_num
                array_num=11;           %阵元个数

                %% -------------dolph-Chebyshev加权--------------------------%%        
                %-------指定旁瓣级（指定旁瓣级和主瓣宽度只能二选一）-----------------------%
%                 sidelobe=-25;                              %旁瓣级，当其不为变量时取值为-25dB
                sidelobe=wide_c(coff_num);                %旁瓣级为变量
        
                R=10.^(-sidelobe/10); 
                z=cosh(1/(array_num-1)*acosh(sqrt(R)));  %论文中式2-4中z的定义
            
                %----------------------指定主瓣宽度--------------------------------------%
%                 mainlobe_wide=wide_c(coff_num);          %主瓣宽度            
%                 sin_theta_NN=sin(mainlobe_wide/2*pi/180);       
%         
%                 z=cos(pi/(2*(array_num-1)))/cos(pi*d/lamda*sin_theta_NN);

                %---------------------计算切比雪夫加权系数------------------------------------%
                wh=zeros(array_num,1);       %array_num代表的是阵元个数即式2-4中M的含义

                wh(1)=(z^(array_num-1))/2;

                for m=2:floor(array_num/2+1)  %floor函数的作用是取该数组某一个元素的小于等于它的最大整数
                        temp=0;               %temp作为一个临时变量，用于存放计算出来的权值；
                        for i=1:(m-1)
                                temp=temp+0.5*(array_num-1)*factorial_0(m-2)*factorial_0(array_num-i-1)*(z^(array_num-2*m+1))*((z^2-1)^(m-i))/...
                                        factorial_0(m-i)/factorial_0(i-1)/factorial_0(m-i-1)/factorial_0(array_num-m);  
                        end             %factorial_0的作用是求一个数的阶乘，这一句实现的是式2-4中2<=m<(M/2+1)的情况下权值的取值。

                        wh(m)=temp;
                end

                for m=ceil(array_num/2+1):array_num   %ceil函数取的是大于等于该数的最小整数
                        wh(m)=wh(array_num+1-m);      %wh（权值矩阵）是一个对称的矩阵，它的元素个数为10个，
                                                      %会有wh(1)=wh(10),wh(2)=wh(9),......,wh(5)=wh(6),是从中间对称的；
                end                                    

                %-----------------加波束指向的加权矢量---------------------------%
                beam_point=0;        %指定波束指向角度，这个角度也可以为变量；

                v0=(exp(-1j*2*pi*f0*(0:(array_num-1))*d*sin(beam_point*pi/180)/c)).'/array_num;  %.'的意思是非共轭转置；v0为导向矢量

                w=wh.*v0;                           %加波束指向的加权矢量

                %% ----------------波束形成-----------------------------%
                beam=zeros(1,length(ang));    %beam是波束响应矩阵

                for i=1:length(ang)
                        v=(exp(-1j*2*pi*f0*(0:(array_num-1))*d*sin(ang(i)*pi/180)/c)).';
                        beam(i)=w'*v;  %计算对应于每一个角度上的波束响应（阵增益）值。
                end

                beam_abs=abs(beam)/max(abs(beam));      %对波束响应取绝对值；
                                %% --------------画图-------------------------------%%
                
                figure(1);
%                 grid on;
                hold on;
                
                plot(ang,20*log10(beam_abs));    %画出波束图
                axis([-90,90,-80,5]);
%         
%                 plot(ang,-6*ones(size(ang)),'m--');              %画出旁瓣级水平线
% %         
%                 plot(45*ones(1,121),-100:20,'r');                     %画出波束指向角的位置
%         
%                 plot((beam_point-ang_temp(d_num)/2).*ones(1,81),-80:0,'r--');       
%          
%                 plot((beam_point+ang_temp(d_num)/2)*ones(1,81),-80:0,'r--');        %画出主瓣波束两零点的位置即标出波束主瓣的位置；
%         
%                 axis([min(ang),max(ang),-80,0]);
%                 
                xlabel('方位角/(Degrees)');
                ylabel('波束增益/(dB)');
%                 title('期望波束');

                %----------------y轴由……表示------------------------%
%                 coff_y(coff_num)=ang_temp(d_num);               %coff_y记录的是每一次不同的旁瓣级所对应的波束两零点间束宽水平
%                  
%                 %----------------x轴由……表示------------------------%
%                 coff_x(coff_num)=sidelobe_level;         % coff_x记录的是每一次的旁瓣级水平

%                 coff_x(coff_num)=array_num;
        
%                 coff_x(coff_num)=beam_point;

                %------------------确认波束3dB带宽-----------------%

                IndMin=find((abs(beam_abs-((0.5)*beam_abs(901))))<0.02);
                temp=ang(IndMin);
                Ang_3dB=min(abs(temp-beam_point))*2;
                Ang_3dB_SOCP_A(1,d_num)=Ang_3dB;
               

                
                
                %--------------确认两零点间束宽-------------------------%
%                 MainLobe_down=find(diff(sign(diff(beam_abs)))>0)+1;  %diff函数的作用是对一个数组进行微分（差分）运算；
%                                                               %sign函数的作用是判断一个值的正负情况：正数返回值为1，负数返回值为-1,0则返回0；
%                                                               %find函数的作用是用于返回满足括号中条件的元素的所在位置(位置的判定：在矩阵中，第一列开始，
%                                                               %自上而下，依次为1，2，3...,然后再从第二列，第三列依次往后数)
%                                                               %此处整个语句的作用是：寻找波束响应的凹陷点即波束图中谷底的位置对应的角度下标；
%                 temp=ang(MainLobe_down);                      %这里的temp指代的是波束图凹陷点所对应的角度值；
%                 temp=temp-beam_point*ones(1,length(temp));    
%                 ang_temp(d_num)=min(abs(temp))*2;                    %这里的ang_temp所对应的就是波束主瓣的宽度即波束两零点间束宽；
%                 
%                 %--------------确认旁瓣级-------------------------%
%                 temp=find(ang<(beam_point-ang_temp(d_num)/2));          %这里寻找的是角度小于-13.5°的角度的下标值
%                 sidelobe_level=20*log10(max(beam_abs(temp))) ;   %这里计算的是旁瓣水平即旁瓣级的最大值的分贝表示值；
%                 sl_FIB(d_num)=sidelobe_level;                           %sl_FIB指的是旁瓣水平的最大值
% 
%                 %--------------确认旁瓣级束宽-------------------------%
%                 temp=find((20*log10(beam_abs)-sidelobe_level*ones(size(beam_abs)))>0.001);%这里的0.001是人为设定的一个阈值，
%                                                                                           %语句的作用是用来寻找波束响应值比旁瓣级大的角度的下标；
%                 temp=ang(temp);
%                 SL_wide(d_num)=temp(end)-temp(1)     %SL_wide即为旁瓣级束宽，temp表示的是主瓣水平大于旁瓣级水平的角度的范围；
                


        end
    end
        figure(2);
        plot(D/lamda,Ang_3dB_SOCP_A,'r*:');hold on; 
        xlabel('阵元间距/{\lambda}');
        ylabel('波束主瓣3dB宽度/(Degrees)');
%         grid on;
%         plot(D,SL_wide,'bx-.');hold on; grid on;
    
%         figure(2);hold on;box on;
% 
%         plot(coff_x,coff_y,'k-o');
%         grid on;
%         h = gca;
%         set(h,'FontSize',10,'FontName','宋体');
%         set(h,'FontName','Times New Roman');
% 
%         xlabel('旁瓣水平(dB)');
% %         xlabel('波束指向角(度)');
% %         xlabel('阵元个数');
%         
%         ylabel('波束主瓣宽度(度)');
%         axis([-40.1,-20,26,44]);
%         title('主瓣宽度与旁瓣级的关系');
%         % title('主瓣宽度与阵元数的关系');

%         % legend('SOCP','Dolph-chebyshev');
%         % legend('间距为二分之一波长','间距为四分之一波长');
