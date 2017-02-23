       

        %% ****************************************************** %%
         %                 方向不变恒定束宽窄带一维波束形成  
         %  此程序是基于二阶锥规划（SOCP）的方向不变恒定束宽窄带一维波束形成，所利用的方
         %  法是文献26中提出的，它的基本思想是：
         %    在保证各个指向上所设计波束的主瓣与相应的参考波束的主瓣之间的均方误差最小的
         %    情况下，使得所设计的波束的旁瓣水平最低，从而实现在各个方向上具有恒定束宽。
         %    这个方法主要针对窄带波束形成器的波束主瓣宽度会随波束指向角的增大而逐渐变宽
         %    的问题，能够使在一定波束指向角范围内的波束主瓣宽度保持恒定，聚焦于一维波束
         %    形成即指向角为方位角。
         %      
        %% ****************************************************** %%         

        clc;
        clear;
%         close all

        %% ******************* 参数设置 ******************* %%
        array_num=13;                   %设置阵元数

        f0=20000;                       %设置最高频率即参考频率，由于这里只针对方向，故不考虑频率不变性
        C=340;                          %声速为340m/s
        array_space=C/2/f0;             %阵元间距取为参考频率对应波长的一半，保证满足采样定理
        
        sidelobe_level=-35;             %参考波束旁瓣级

        ang_reference=0;               %设置参考波束指向角
        mainbeam_width=25;              %设置参考波束主瓣宽度，这里用的是旁瓣级束宽

        sidelobe=[-90:1:(ang_reference-mainbeam_width/2),(ang_reference+mainbeam_width/2):1:90];  
                                        %确定旁瓣区域取值范围并将其离散化

        silelobe_num=length(sidelobe);  %旁瓣离散化之后的个数

        ang=-90:0.1:90;                   %设置方位角的范围

        beam=zeros(1,length(ang));      %初始化波束响应矩阵

        %% ******************* 基于SOCP设计参考波束响应 ******************* %%
        
        matrix_zero=zeros(array_num*2,1);           %设置一个初始26*1维零矩阵

        % ******************* 设置矩阵b ******************* %
        
        b=([-1,matrix_zero.']).';        %设置矩阵b，为一个27*1维的矩阵

        % ******************* 设置二阶锥约束的维数 ******************* %
        K.f=2;                          
        K.q=[3*ones(1,silelobe_num),(2*array_num+1)];

        % ******************* 矩阵C的转置 ******************* %
        
        c1=([1;0]).';                  %c1为1*2维 
        c2=zeros(1,3*silelobe_num);    %c2为1*456维，初始化为0矩阵 
        for i=1:silelobe_num
                c2(1,1+3*(i-1):3*i)=([10^(sidelobe_level/20);0;0]).';
        end                            %这一句对c2矩阵进行赋值，它的每三个值会重复

        c3=([0;matrix_zero]).';        %c3为1*27的零矩阵

        % ******************* 计算矩阵A ******************* %
        
        a3=([-1,matrix_zero.';matrix_zero,-eye(array_num*2,array_num*2)]).'; %a3为27*27的负单位矩阵
                                            
        v_0=exp(-1j*2*pi*f0*(0:(array_num-1))*array_space*sin(ang_reference*pi/180)/C);
                                             %V_0表示在参考频率下参考方向上的阵列的导向矢量，为1*13维
        va_0=[real(v_0),imag(v_0)].';        %va_0为V_0的实部与虚部所构成的矩阵，为26*1维
        
        vb_0=[-imag(v_0),real(v_0)].';       %vb_0为V_0的实部与虚部所构成的矩阵，为26*1维，与va_0有些微区别 

        a1=([0,va_0.';0,-vb_0.']).';         %a1的维度为27*2，a1对应于参考方向上的波束的形成 

        a2=zeros((array_num*2+1),3*silelobe_num);     %a2为27*456维的0矩阵 

        for i=1:silelobe_num
                v_i=exp(-1j*2*pi*f0*(0:(array_num-1))*array_space*sin(sidelobe(i)*pi/180)/C);
                va_i=[real(v_i),imag(v_i)].';
                vb_i=[-imag(v_i),real(v_i)].';

                a2(1:(array_num*2+1),1+3*(i-1):3*i)=([0,matrix_zero.';0,-va_i.';0,-vb_i.']).';
        end                                         %a2为对应于旁瓣离散化区域每一个方位角度的波束的形成 

        % ******************* 问题的求解 ******************* %
                      % 利用前面赋值好的相关的a,b和c矩阵以及sedumi参数，调用sedumi工具箱（函数）进行权值的求解 
        
        C_col=([c1,c2,c3]).';       %定义C_col矩阵

        A_col=([a1,a2,a3]).';       %定义A_col矩阵
        
        [x,y,info]=sedumi(A_col,b,C_col,K);     %调用sedumi工具箱（函数）求解上述的二阶锥规划（凸优化）问题，
                                                %参考波束的加权值存于y中，输出的y为27*1维

        wr=y(2:array_num+1)+1j*y((array_num+2):(array_num*2+1));    %得到参考方向上的波束加权值，是从y中取出的
        
        % ******************* 波束形成 ******************* %
        
        for i=1:length(ang)
                beam(i)=conj(exp(-1j*2*pi*f0*(0:(array_num-1))*array_space*sin(ang(i)*pi/180)/C))*wr;
        end                                     %计算在每一个离散角度上的波束响应值

        beam_abs=abs(beam)/max(abs(beam));      %对求得的波束响应矩阵求绝对值并且利用其最大值进行归一化，以方便画图

        % ******************* 画图 ******************* %   
        
        figure(1);hold on;box on;
        
        plot(ang,20*log10(beam_abs),'r:');

        plot(ang_reference*ones(1,81),-80:0,'g');
        plot(ang,sidelobe_level*ones(size(ang)),'k--');
        
        axis([min(ang),max(ang),-80,0]);

        xlabel('方位角(度)');
        ylabel('波束(dB)');
 %         title('期望波束');   
 
        h = gca;
        set(h,'FontSize',10,'FontName','宋体');
        set(h,'FontName','Times New Roman');
        
%         legend('40度','60度')

       %% ******************* 方向不变波束形成设计 ******************* %%
        
        % ******************* 初始化对应于25个波束指向角的加权矢量和波束响应矩阵
        ang_range=-90:5:90;                                %给定波束指向角的范围，每5°取一个值
        multi_w=zeros(array_num,length(ang_range));        %不同波束指向角对应的加权矢量矩阵，为13*25维
        beam_multi=zeros(length(ang_range),length(ang));   %不同波束指向角对应的波束响应矩阵，为25*181维
        MainWidth=zeros(1,length(ang_range));
        % ******************* 角度信息 ******************* %
        
        for ang_num=1:length(ang_range)           %循环从1—25,本循环实现的是对所选取的-60到60度范围内
                                                  %每一个波束指向角上的波束形成，方法仍然是SOCP
                beam_point=ang_range(ang_num);    %波束指向角分别取-60：5：60°，共25个波束指向角
                
                mainlobe=(beam_point-mainbeam_width/2):1:(beam_point+mainbeam_width/2);  
                                                  %将主瓣以1度进行离散化，以便于后面将平移后的波束主瓣进行一一对比                
                sidelobe=[-90:1:(beam_point-mainbeam_width/2),(beam_point+mainbeam_width/2):1:90];  
                                                  %对每一个主瓣所对应的旁瓣以1度进行离散化

                mainlobe_num=length(mainlobe(1,:));   %确定主瓣离散数
                silelobe_num=length(sidelobe);        %确定旁瓣离散数

                % ******************* 主瓣期望响应 ******************* %
                beam_d=zeros(1,mainlobe_num);               %主瓣期望响应
               
                for i=1:length(mainlobe)   %这个循环的作用是求解出在该波束指向角每一个主瓣离散角度上的波束响应值
                                           %用来作为之后约束条件的书写条件
                       beam_d(i)=conj(exp(-1j*2*pi*f0*(0:(array_num-1))*array_space*sin((mainlobe(i)-(beam_point-ang_reference))*pi/180)/C))*wr;
                end
                beam_d=beam_d/max(abs(beam_d));    %对该波束指向角上的波束响应进行归一化

        %% ******************* 二阶锥规划计算针对该波束指向角的加权值 ******************* %%
                % ******************* 控制波束指向角上的波束响应为1 ******************* %
                cc1=([1;0]).';                      %确定cc1
                 
                v_0=exp(-1j*2*pi*f0*(0:(array_num-1))*array_space*sin(beam_point*pi/180)/C);
                                                    %v_0表示该波束指向角度上的导向矢量  
                va_0=[real(v_0),imag(v_0)].';
                vb_0=[-imag(v_0),real(v_0)].';      %va_0和vb_0都由v_0的实部和虚部的组合构成
                
                aa1=([0,zeros(1,mainlobe_num),va_0.';0,zeros(1,mainlobe_num),-vb_0.']).';   %构造aa1矩阵
                
                % ******************* 控制主瓣的误差和小于某个数 ******************* %
                cc2=10;       
                aa2=([0,ones(1,mainlobe_num),zeros(1,2*array_num)]).';       %构造aa2矩阵
               
                % ******************* 控制主瓣 ******************* %
                cc3=zeros(1,4*mainlobe_num);         
                aa3=zeros((array_num*2+mainlobe_num+1),4*mainlobe_num);

                % ******************* 控制旁瓣 ******************* %
                cc4=zeros(1,3*silelobe_num);
                aa4=zeros((array_num*2+mainlobe_num+1),3* silelobe_num);
                
                % ******************* 控制稳健性 ******************* %
                cc5=([0.4;matrix_zero]).';
                aa5=([zeros(1,mainlobe_num+1),matrix_zero.';zeros(2*array_num,mainlobe_num+1),-eye(array_num*2,array_num*2)]).';   
                        
                % ******************* 矩阵B ******************* %
                bb=([-1,zeros(1,mainlobe_num),matrix_zero.']).';

                % ******************* 二阶锥约束的维数 ******************* %
                Q.f=2;
                Q.l=1;
                Q.q=[4*ones(1,mainlobe_num),3*ones(1,silelobe_num)];

                           
                % ******************* 矩阵C的转置和矩阵A ******************* %

                 for i=1:mainlobe_num
                      cc3(1,1+4*(i-1):4*i)=([1;2*real(beam_d(i));2*imag(beam_d(i));-1]).';

                      v_i=exp(-1j*2*pi*f0*(0:(array_num-1))*array_space*sin(mainlobe(i)*pi/180)/C);
                      va_i=[real(v_i),imag(v_i)].';
                      vb_i=[-imag(v_i),real(v_i)].';

                      q_one=zeros(mainlobe_num,1);
                      q_one(i)=1;

                      aa3(1:(array_num*2+mainlobe_num+1),1+4*(i-1):4*i)=([0,-q_one.',matrix_zero.';0,zeros(1,mainlobe_num),2*va_i.';0,zeros(1,mainlobe_num),2*vb_i.';0,-q_one.',zeros(1,2*array_num)]).';
                 end


                 for i=1:silelobe_num
%                        cc4(1,1+3*(i-1):3*i)=([10^(-5/20);0;0]).';
                       cc4(1,1+3*(i-1):3*i)=([0;0;0]).';
                       v_i=exp(-1j*2*pi*f0*(0:(array_num-1))*array_space*sin(sidelobe(i)*pi/180)/C);
                       va_i=[real(v_i),imag(v_i)].';
                       vb_i=[-imag(v_i),real(v_i)].';

                       aa4(1:(array_num*2+mainlobe_num+1),1+3*(i-1):3*i)=([-1,zeros(1,mainlobe_num),matrix_zero.';0,zeros(1,mainlobe_num),-va_i.';0,zeros(1,mainlobe_num),-vb_i.']).';
                 end

                % ******************* 问题的求解 ******************* %
                CC_col=([cc1,cc2,cc3,cc4]).';

                AA_col=([aa1,aa2,aa3,aa4]).';

                [x,y,info]=sedumi(AA_col,bb,CC_col,Q);

                ww=y(mainlobe_num+2:mainlobe_num+array_num+1)+1j*y((mainlobe_num+array_num+2):(mainlobe_num+array_num*2+1));
                      
                multi_w(:,ang_num)=ww;
                
                % ******************* 波束形成 ******************* %     
                for i=1:length(ang)
                      beam(i)=conj(exp(-1j*2*pi*f0*(0:(array_num-1))*array_space*sin(ang(i)*pi/180)/C))*ww;
                end
                
                beam_abs=abs(beam)/max(abs(beam));  %归一化
        
                beam_multi(ang_num,:)=beam_abs;
               
                %-----------------确定两零点束宽----------------------------%
                IndMin=find(diff(sign(diff(beam_abs)))>0)+1;
                temp=ang(IndMin);
                temp=temp-beam_point*ones(1,length(temp));    
                mainlobe_wide=min(abs(temp))*2;
                MainWidth(1,ang_num)=mainlobe_wide;
                
                %----------------------确认旁瓣级----------------------------%
%                 temp=find(ang>(beam_point+mainlobe_wide/2));
                temp=find(ang<(beam_point-mainlobe_wide/2));
                sidelobe=20*log10(max(beam_abs(temp)));   
     
               
        end
        %% -------------------画图--------------------------------%   
        figure(2);hold on;box on;

        plot(ang,20*log10(beam_abs),'k:');

        plot(beam_point*ones(1,81),-80:0,'k--');
        plot(ang,sidelobe*ones(size(ang)),'k--');

        axis([min(ang),max(ang),-80,0]);
        
        xlabel('方位角(度)');
        ylabel('波束(dB)');
        title('恒定束宽波束');
        
        h = gca;
        set(h,'FontSize',10,'FontName','宋体');
        set(h,'FontName','Times New Roman');
        
%         legend('40度','60度')
        

        figure(3);
        [x,y]=meshgrid(ang,ang_range);
        mesh(y,x,20*log10(beam_multi));
        
        xlabel('波束指向角(度)');
        ylabel('方位角(度)');
        zlabel('波束(dB)');
        
        h = gca;
        set(h,'FontSize',10,'FontName','宋体');
        set(h,'FontName','Times New Roman');   
        
        figure(4);
        plot(ang_range,MainWidth,'r*-');
        axis([-65,65,20,40]);
        
        dlmwrite('multicoff.txt',multi_w,';');   %需保存加权矢量时才采用
         
