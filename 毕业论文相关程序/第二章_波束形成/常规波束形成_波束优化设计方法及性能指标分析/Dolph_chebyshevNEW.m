        
        %************************************************************%
        %             Dolph-Chebyshev加权窄带波束形成方法             %
        %        copyright Liao Fengyi                               %
        %        last modified at 2015.11.24 （Done）                %
        %************************************************************%

        clc;
        clear all

                     %% ************************参数设置*************************** %%
%         fh=12000;                       %最高频率
%         fl=4000;                        %最低频率
        fr=4000;                        %参考频率

        c=340;                          %声速
        d=c/fr/4;                       %阵元间距为参考频率波长的一半
        lamda=c/fr;                     %参考频率对应的波长
        ang=-90:0.1:90;                 %方位角范围为-90°到90°

%         wide_c=lamda/8:lamda/8:lamda;
          wide_c=-10:-5:-50;                    %旁瓣级取值

%         wide_c=[15,20,25,30,35,40,45,50,55];        %两零点间束宽取值   

%         wide_c=0:5:50;    %波束指向取值

%         wide_c=5:5:30;                     %阵元个数取值
          
          
          N=10;
%           a=length(wide_c);
                     %% ************************Dolph-Chebyshev加权方法设计窄带波束形成器************************ %%
for  n=1:1:N                     
        for coefficient_number=1:length(wide_c)       %系数个数设为coefficient_number
%                 array_number=wide_c(coefficient_number);            %阵元个数
                
                array_number=15;                                    %阵元个数
                     %% ************************dolph-Chebyshev加权************************ %%        
                %**************指定旁瓣级（指定旁瓣级和主瓣宽度只能二选一）***************%
%                 sidelobe=-30;                                     %旁瓣级，当其不为变量时取值为-25dB
                sidelobe=wide_c(coefficient_number);                %旁瓣级取值
        
                R=10.^(-sidelobe/10); 
                z=cosh(1/(array_number-1)*acosh(sqrt(R)));          %论文中式2-4中z的定义
            
                %*************指定主瓣宽度************%
%                 mainlobe_wide=wide_c(coefficient_number);         %主瓣宽度取值            
%                 sin_theta_NN=sin(mainlobe_wide/2*pi/180);         
%         
%                 z=cos(pi/(2*(array_number-1)))/cos(pi*d/lamda*sin_theta_NN);

                    %% *************计算切比雪夫加权系数************* %%
                wh=zeros(array_number,1);             %array_number代表的是阵元个数即式2-4中M的含义，wh存放的是切比雪夫加权的系数，个数与阵元个数相同

                wh(1)=(z^(array_number-1))/2;         %当m=1时的加权值的计算

                for m=2:floor(array_number/2+1)       %floor函数的作用是取该数组某一个元素的小于等于它的最大整数，这儿指的是2<=m<(M/2+1)的情况
                        temp1=0;                      %temp1作为一个临时变量，用于存放计算出来的权值；
                        for i=1:(m-1)
                                temp1=temp1+0.5*(array_number-1)*factorial_0(m-2)*factorial_0(array_number-i-1)*(z^(array_number-2*m+1))*((z^2-1)^(m-i))/...
                                        factorial_0(m-i)/factorial_0(i-1)/factorial_0(m-i-1)/factorial_0(array_number-m);  
                        end                           %factorial_0的作用是求一个数的阶乘，这一句实现的是式2-4中2<=m<(M/2+1)的情况下权值的取值。

                        wh(m)=temp1;
                end

                for m=ceil(array_number/2+1):array_number       %ceil函数取的是大于等于该数的最小整数
                        wh(m)=wh(array_number+1-m);             %wh（权值矩阵）是一个对称的矩阵，它的元素个数为10个，会有wh(1)=wh(10),wh(2)=wh(9),......,wh(5)=wh(6),是从中间对称的；
                                                               
                end                                    

                   %% *************加波束指向的加权矢量**************** %%
%               beam_point=wide_c(coefficient_number);
                beam_point=0;               %指定波束指向角度，这个角度也可以为变量；
%                 d=wide_c(coefficient_number);

                v0=(exp(-1j*2*pi*fr*(0:(array_number-1))*d*sin(beam_point*pi/180)/c)).'/array_number;  %.'的意思是非共轭转置；v0为阵列的导向矢量

                w=wh.*v0;                   %加波束指向的加权矢量

                  %% **************波束形成*************** %%
                beam=zeros(1,length(ang));             %beam是波束响应矩阵

                for i=1:length(ang)                    %计算对应于每一个角度上的波束响应（阵增益）值。
                        v=(exp(-1j*2*pi*fr*(0:(array_number-1))*d*sin(ang(i)*pi/180)/c)).';
                        beam(i)=w'*v;             
                end

                beam_abs=abs(beam)/max(abs(beam));     %对波束响应取绝对值并归一化

                  %% **************确认两零点间束宽************** %%
                MainLobe_down=find(diff(sign(diff(beam_abs)))>0)+1;  %diff函数的作用是对一个数组进行微分（差分）运算；
                                                                     %sign函数的作用是判断一个值的正负情况：正数返回值为1，负数返回值为-1,0则返回0；
                                                                     %find函数的作用是用于返回满足括号中条件的元素的所在位置(位置的判定：在矩阵中，第一列开始，
                                                                     %自上而下，依次为1，2，3...,然后再从第二列，第三列依次往后数)
                                                                     %此处整个语句的作用是：寻找波束响应的凹陷点即波束图中谷底的位置对应的角度下标；
                temp2=ang(MainLobe_down);                            %这里的temp2指代的是波束图凹陷点所对应的角度值；
                temp3=temp2-beam_point*ones(1,length(temp2));        %temp3与temp2代表的都是波束凹陷点的角度值；
                ang_temp=min(abs(temp3))*2;                          %这里的ang_temp所对应的就是波束主瓣的宽度即波束两零点间束宽；
                
                 %% **************确认旁瓣级************** %%
                temp4=find(ang<(beam_point-ang_temp/2));             %这里寻找的是角度小于波束两零点间束宽一半宽度的角度的下标值
                sidelobe_level=20*log10(max(beam_abs(temp4))) ;      %这里计算的是旁瓣水平即旁瓣级的最大值的分贝表示值；
                sidelobe_max=sidelobe_level;                         %sidelobe_max指的是旁瓣水平的最大值

                 %% **************确认旁瓣级束宽************** %%
                temp5=find((20*log10(beam_abs)-sidelobe_level*ones(size(beam_abs)))>0.0001);%这里的0.001是人为设定的一个阈值，语句的作用是用来寻找波束响应值比旁瓣级大的角度的下标；
                                                                                          
                temp6=ang(temp5);
                SL_wide=temp6(end)-temp6(1);                          %SL_wide即为旁瓣级束宽，temp6表示的是主瓣水平大于旁瓣级水平的角度的范围；
                
                 %% **************画图************** %%
                
%                 figure(1);grid on;
%                 hold on;
%                 
%                 plot(ang,20*log10(beam_abs));    %画出波束图
%         
%                 plot(ang,sidelobe_level*ones(size(ang)),'m--');              %画出旁瓣级水平线
%         
%                 plot(beam_point*ones(1,81),-80:0,'r');                     %画出波束指向角的位置
%         
%                 plot((beam_point-ang_temp/2).*ones(1,81),-80:0,'r--');       
%          
%                 plot((beam_point+ang_temp/2)*ones(1,81),-80:0,'r--');        %画出主瓣波束两零点的位置即标出波束主瓣的位置；
%         
%                 axis([min(ang),max(ang),-80,0]);
%                 
%                 xlabel('方位/度');
%                 ylabel('波束/dB');
%                 title('期望波束');

                %**************y轴由……表示**************%
                coff_y(coefficient_number,n)=ang_temp;               %coff_y记录的是每一次不同的旁瓣级所对应的波束两零点间束宽水平
                 
                %**************x轴由........表示**************%
                                           % coff_x记录的是每一次的旁瓣级水平/阵元个数/波束指向角
                coff_x(coefficient_number,n)=sidelobe_level;         

%                 coff_x(coefficient_number,n)=array_number;
        
%                 coff_x(coefficient_number,n)=beam_point;
                
%                 coff_x(coefficient_number,n)=d/lamda;

        end
        figure(1);grid on;
                hold on;
                
                plot(ang,20*log10(beam_abs));    %画出波束图
        
%                 plot(ang,sidelobe_level*ones(size(ang)),'m--');              %画出旁瓣级水平线
        
                plot(beam_point*ones(1,81),-80:0,'r');                     %画出波束指向角的位置
        
                plot((beam_point-ang_temp/2).*ones(1,81),-80:0,'r--');       
         
                plot((beam_point+ang_temp/2)*ones(1,81),-80:0,'r--');        %画出主瓣波束两零点的位置即标出波束主瓣的位置；
        
                axis([min(ang),max(ang),-80,0]);
                
                xlabel('Angle/(Degree)');
                ylabel('Beam Gain/(dB)');
                title('Desired Beam Pattern');
                h = gca;
                set(h,'FontSize',10,'FontName','Times New Roman');
                set(h,'FontName','Times New Roman');
     n=n+1;
end        
        coff_Y=sum(coff_y,2);
        coff_Y=coff_Y./N
        coff_X=sum(coff_x,2);
        coff_X=coff_X./N
        figure(2);hold on;box on;

        plot(coff_X,coff_Y,'k-*');
%         plot(coff_X,coff_Y,'r-*');
%         grid on;
        h = gca;
        set(h,'FontSize',10,'FontName','Times New Roman');
        set(h,'FontName','Times New Roman');

        xlabel('Sidelobe Level(dB)');
%         xlabel('Sidelobe Level(Real Value)');
%         xlabel('Beam Point(Degree)');
%         xlabel('Array Number');
%           xlabel('Array Element Space');
        
        ylabel('Mainlobe Width(Degree)');
%         axis([-40.1,-20,26,44]);
%         axis([-51,-9,9,26]);
        title('Relationship Between Mainlobe Width And Sidelobe Level');
%         title('Relationship Between Mainlobe Width And Beam Point');
%         title('Relationship Between Mainlobe Width And Array_number');
%         axis([7,16,18,42]);
        % legend('SOCP','Dolph-chebyshev');
        % legend('间距为二分之一波长','间距为四分之一波长');
