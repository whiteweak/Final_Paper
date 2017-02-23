
       %*************************************************************%
       %          基于二阶锥规划（SOPC）的波束优化设计               %
       %          波束的设计原理是稳健旁瓣控制高增益波束设计         %
       %   copyright Liao Fengyi                                     %
       %   last modified at 2015.11.24 (done)                        %

        clc;
        clear;

                     %% ****************参数设置**************** %%
%         fl=4000;                                %最低频率
%         fh=12000;                               %最高频率
        fr=20000;                                %参考频率
        C=340;                                  %声速

%         wide_c=0:5:20;                                          %波束指向角度取值
%         wide_c=8:1:15;                                          %阵元个数取值
%         wide_c=-14:-1:-32;               %旁瓣水平取值
        wide_c=-40;               %旁瓣水平取值

%         mainlobe_c=20.2:0.8:34.6;                      %半波长时对应的旁瓣级束宽
        mainlobe_c=25;                      %半波长时对应的旁瓣级束宽
        
%         mainlobe_c=37:2.8:65;                     %四分之一波长时对应的旁瓣级束宽
                     %% ****************波束优化设计**************** %% 
        for coefficient_number=1:length(wide_c)               %coefficient_number的数值为1到5，为旁瓣水平的取值数
            sidelobe_level=wide_c(coefficient_number);        %旁瓣级水平取值为变量
%             sidelobe_level=-25;
            
            array_number=15;                                  %阵元数 
%             array_space=C/4/fr; 
            array_space=C/2/fr;                               %阵元间距为四分之一波长，波长为C/fr

            beam_point=0;
%             beam_point=wide_c(coefficient_number);                                    %波束指向角度设置为0°

%             mainbeam_wide=29;
            mainbeam_wide=mainlobe_c(coefficient_number);     %旁瓣级束宽取值即主瓣宽度

            silelobe=[-90:1:(beam_point-mainbeam_wide/2),(beam_point+mainbeam_wide/2):1:90];  %以1度的间距对旁瓣区域进行角度上的离散化(针对方位角)
 
            silelobe_number=length(silelobe);        %旁瓣离散化后的数目，这里的sidelobe_num代表的是SOCP中的Nsl

            ang=-90:0.1:90;                          %设定方位角的变化（取值）范围

            beam=zeros(1,length(ang));               %建立波束矩阵，矩阵维数为1x361，这是针对于每一个离散角度（0.5度）来说的

            matrix_zero=zeros(array_number*2,1);     %建立一个（2倍阵元个数）x1的0矩阵，对应的是公式中的零矩阵
            
            
            
                     %% ***************计算用SeDuMi工具箱来计算SOCP加权值所需的相关参数，为整个算法的核心部分***************** %%
                     
            %*************矩阵b*************%
            b=([-1,matrix_zero.']).';                %建立一个（2倍阵元个数+1）x1的矩阵（列向量），作为公式中的b矩阵

            %*************二阶锥约束的维数*************%
            K.f=2;                                   % K.f指的是sedumi工具中，K结构中描述的对应于“无约束的”的维数个数，
                                                     % 其值表示的是整个约束条件中有几个等式约束

            K.q=[3*ones(1,silelobe_number),(2*array_number+1)];      % K.q确定的是sedumi工具中，洛伦兹约束的维度，这个值是与公式中

            %*************矩阵C的转置*************%
            c1=([1;0]).';                            %“.'”的作用是求解该矩阵的非共轭转置,c1的维度为1x2；

            c2=zeros(1,3*silelobe_number);           %定义一个c2矩阵，出初始元素为0
            for i=1:silelobe_number                  %此循环的作用是对c2进行相应的赋值操作
                    c2(1,1+3*(i-1):3*i)=([10^(sidelobe_level/20);0;0]).';   %给c2矩阵赋值,c2矩阵指代的是公式中C2+i（i=0,1,2,...Nsl）矩阵      
            end
                
          
            c3=([0;matrix_zero]).';                 %给c3矩阵赋值

            %*************矩阵A*************%
            a3=([-1,matrix_zero.';matrix_zero,-eye(array_number*2,array_number*2)]).';       % eye产生一个22x22的单位阵，本语句执行完后会产生一个标准的负的单位阵即a3；

            v_0=exp(-1j*2*pi*fr*(0:(array_number-1))*array_space*sin(beam_point*pi/180)/C);  %array_space为阵元间距，相当于d；v_0为波束指向角为0度时的导向矢量矩阵，
                                                                                             %此处为一1x11的全1矩阵；
            va_0=[real(v_0),imag(v_0)].';      
            vb_0=[-imag(v_0),real(v_0)].';                      %分别构造由导向矢量的实部和虚部构成的矩阵va_0和vb_0

            a1=([0,va_0.';0,-vb_0.']).';                        %定义a1矩阵

            a2=zeros((array_number*2+1),3*silelobe_number);     %定义a2矩阵

            for i=1:silelobe_number      % i从1到116
                    v_i=exp(-1j*2*pi*fr*(0:(array_number-1))*array_space*sin(silelobe(i)*pi/180)/C);          %计算经过离散化之后的旁瓣区域
                                                                                                              %各角度上的导向矢量v_i；
                    va_i=[real(v_i),imag(v_i)].';
                    vb_i=[-imag(v_i),real(v_i)].';              %分别构造由导向矢量的实部和虚部构成的矩阵va_i和vb_i；

                    a2(1:(array_number*2+1),1+3*(i-1):3*i)=([0,matrix_zero.';0,-va_i.';0,-vb_i.']).';         %构造a2矩阵；
            end

            %*************二阶锥规划问题的求解（利用matlab的工具箱）*************%
            C_col=([c1,c2,c3]).';

            A_col=([a1,a2,a3]).';

            [x,y,info]=sedumi(A_col,b,C_col,K);    %利用sedumi工具求解二阶锥规划问题的相关参数

            w=y(2:array_number+1)+1j*y((array_number+2):(array_number*2+1));   %w是计算出来的权值矩阵，它的取值是y矩阵的2—Nsl+1个元素
            
            
            
                     %% ****************计算利用设计出来的权值加权的波束形成***************       
            %*************波束形成*************%
            for i=1:length(ang)
                    beam(i)=conj(exp(-1j*2*pi*fr*(0:(array_number-1))*array_space*sin(ang(i)*pi/180)/C))*w;           %beam为加权波束响应矩阵
            end
%             beam1=max(abs(beam));
            beam_abs=abs(beam)/max(abs(beam));              %对波束响应进行归一化
            beam_abs_1=beam_abs;
            %*************确认两零点间束宽*************%
            MainLobe_down=find(diff(sign(diff(beam_abs)))>0)+1;    %diff函数的作用是对一个数组进行微分（差分）运算；
                                                            %sign函数的作用是判断一个值的正负情况：正数返回值为1，负数返回值为-1,0则返回0；
                                                            %find函数的作用是用于返回满足括号中条件的元素的所在位置(位置的判定：在矩阵中，第一列开始，
                                                            %自上而下，依次为1，2，3...,然后再从第二列，第三列依次往后数)
                                                            %此处整个语句的作用是：寻找波束响应的凹陷点即波束图中谷底的位置对应的角度；
            temp1=ang(MainLobe_down);
            temp2=temp1-beam_point*ones(1,length(temp1));    
            ang_temp=min(abs(temp2))*2;                      %这里的ang_temp所对应的就是波束主瓣的宽度即波束两零点间束宽；
            
            %*************旁瓣级与设计指标的逼近程度*************%
            temp3=find(ang<(beam_point-ang_temp/2));
            sl_FIB=20*log10(max(beam_abs(temp3)));  
%             sl_FIB-sidelobe_level
                   
            %% ****************画图**************** %%
            
            figure(1);hold on;grid on;
            
%             plot((beam_point-ang_temp/2).*ones(1,81),-80:0,'r--');
%             plot((beam_point+ang_temp/2)*ones(1,81),-80:0,'r--');

            plot(ang,20*log10(beam_abs));

%             plot(beam_point*ones(1,81),-80:0,'m--');
%             plot(ang,sl_FIB*ones(size(ang)),'m--');
            
            axis([min(ang),max(ang),-100,0]);

            xlabel('方位/度');
            ylabel('波束/dB');
            title('期望波束');
         
            beam_abs=beam_abs*max(abs(beam)); 
            coff_y(coefficient_number)=20*log10(max(beam_abs));
%             coff_y(coefficient_number)=max(beam);
            sidelobe_level=10^(sidelobe_level/20);
            coff_x(coefficient_number)=20*log10(sidelobe_level*max(abs(beam)));
%             coff_y(coefficient_number)=ang_temp;
%             coff_x(coefficient_number)=wide_c(coefficient_number);

        end
        dlmwrite('coff_w_socp.txt',w,';'); 
        beam_abs_1_max=max(beam_abs_1);
        indmin=find(0.5<(beam_abs_1_max-beam_abs_1)<0.6);
        indmin=max(indmin);
        beam_abs_3dbwidth=(ang(indmin))*2;
%         figure(2);hold on;box on;
% 
%         plot(coff_x,coff_y,'r-^');
%         grid on;
%         xlabel('旁瓣级(dB)');
% %         ylabel('主瓣宽度（度）');
% %         title('主瓣宽度与旁瓣级的关系');
%         ylabel('主瓣增益(dB)');
%         title('主瓣最大增益与旁瓣级的关系');
% %         axis([min(wide_c),max(wide_c),-10,10]);