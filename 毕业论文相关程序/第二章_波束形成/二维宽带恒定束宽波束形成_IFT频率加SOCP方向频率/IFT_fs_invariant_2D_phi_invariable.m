        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %    基于多维傅里叶逆变换的宽带频率-方向不变恒定束宽二维波束形成             %
        %    本程序以基于SOCP的二维方向不变恒定束宽波束形成所设计出的波束            %
        %    作为后续基于离散傅里叶逆变换（IFT）的恒定束宽频率不变波束形             %
        %    成的参考波束，进而联合实现宽带方向—频率不变波束形成。
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        clc;
        clear all
        % close all

        %% ************************** 参数设置 ************************** %%
        
        array_numx=15;                  %平面阵X轴上的阵元个数
        array_numy=15;                  %平面阵Y轴上的阵元个数

        fh=20000;                       %参考频率
        fl=fh/2;                      %最低频率
        f0=fh;                          %最高频率

        c=340;                          %传播速度
        d=c/f0/2;                       %阵元间距
        lamda=c/f0;                     %参考频率对应的波长

        theta=pi/180*(0:0.2:90);          %俯仰角扫描范围
        phi=pi/180*(-180:0.2:180);        %方位角扫描范围

        num=15;                         %期望波束采用的阵元数为num*num
        vec=-(num-1)/2:(num-1)/2;


               %% ************************** 频率不变波束形成 ************************** %%
        theta_point=pi/180*60;                       %俯仰角波束指向60度
        
        phi_point=pi/180*0;              %方位角波束指向范围

        multicoff=dlmread('coeffcient_phi0_theta60_all.txt');             %coff_w.txt文件中存放的是在steer_invariant_constent_width_2D_BF程序中
                                                                   %已经计算好的期望波束响应加权值，这里是为了简便运算
                                                                   
               %% ************************** 初始化参考方向不同频率波束响应矩阵以及主瓣宽度矩阵 ************************** %%                                             
        
        OMEGA_range=2*pi*(linspace(fl,fh,20));    %划分子带范围

        theta_beam_plus60=zeros(length(OMEGA_range),length(theta));   %俯仰角方向上的波束响应值
        phi_beam_plus60=zeros(length(OMEGA_range),length(phi));       %方位角方向上的波束响应值
        
        beam_theta_frequency_plus60=zeros(length(theta_point),length(OMEGA_range));
                           %初始化当方向为参考方向时不同频率对应的波束响应矩阵，为31x20维，记录的是频率变化时得到的波束图中每一个子频带下位于参考
                           %波束方向（方位角0°，俯仰角46°）位置的波束响应值
%         beam_theta_frequency_plus40_mainlobewidth=zeros(length(theta_point),length(OMEGA_range));
                           %记录beam_phi_frequency中每个波束的主瓣宽度值，方便后续画图
        

        for angle_num=1:length(theta_point)        
                h1=multicoff(:,angle_num);           %取该矩阵的每一列
                h1=h1';                              %为*15x15）*1的列向量

               % ************************** 计算期望波束响应 ************************** %
                beam_temp=zeros(length(theta),length(phi));   
                for k=1:length(theta)
                        for i=1:length(phi)
                                v_x=(exp(-1j*(vec)*pi*sin(theta(k))*cos(phi(i)))).';
                                v_y=(exp(-1j*(vec)*pi*sin(theta(k))*sin(phi(i)))).';
                                v=kron(v_y,v_x);
                                beam_temp(k,i)=(h1*v);

                        end
                end
                beam_temp_original_plus60=abs(beam_temp);
                beam_temp=abs(beam_temp)/max(max(abs(beam_temp)));  
                      %这里形成的波束已经实现了基于切比雪夫加权的窄带二维参考波束设计以及实现了该波束基于SOCP方法的方向不变性，
                      %用来作为接下来的基于傅里叶逆变换的频率不变波束形成的参考波束，以实现针对每一个方向上的频率不变波束形成
                      %的加权值的设计
                                                                          

                %% ************************** 直角坐标系的波束图 ************************** %
                
                [x,y]=meshgrid(phi*180/pi,theta*180/pi);
                figure(1);           
                mesh(x,y,20*log10(beam_temp));hold on;grid on;

                xlabel('方位角/度');
                ylabel('俯仰角/度');
                zlabel('波束/dB');
                title('期望波束图');

                %% ************************** 球坐标系的波束图 ************************** %
                figure(2);
                Patt3d(beam_temp.',3);hold on;grid on;        

%                 %% ************************** 方位角方向上的波束图，俯仰角为波束指向角 ************************** %
%                 figure(3);       
%                 plot(phi*180/pi,20*log10(beam_temp(find(theta==theta_point),:)),'k:');hold on;box on;
%               
%                 axis([-180,180,-80,0]);
%         
%                 h = gca;
%                 set(h,'FontSize',10,'FontName','宋体');
%                 set(h,'FontName','Times New Roman');
%                 
%                 xlabel('方位角(度)');
%                 ylabel('波束(dB)');
%         %         title('方位角方向上的期望波束图');
%         
%         %         legend('(30,0)','(30,120)');

                %% ************************** 俯仰角方向上的波束图，方位角为波束指向角 ************************** %
                figure(4);
                plot(theta*180/pi,20*log10(beam_temp(:,find(phi==phi_point))),'k:');hold on;box on;
                
                axis([0,90,-80,0]);
                
                h = gca;
                set(h,'FontSize',10,'FontName','宋体');
                set(h,'FontName','Times New Roman');
                
                xlabel('俯仰角(度)');
                ylabel('波束(dB)');
        %         title('俯仰角方向上的期望波束图');
        
        %         legend('(30,0)','(30,120)');


               %% ************************** 基于多维傅里叶逆变换的FIB ************************** %%
               
               
                for OMEGA_num=1:length(OMEGA_range) 

                        Omega=OMEGA_range(OMEGA_num);     

                        M=64;                                         %频率离散点数
                        Omega1=-pi*ones(1,M)+2*pi/M*(0:(M-1));        %频率w1的范围
                        Omega2=-pi*ones(1,M)+2*pi/M*(0:(M-1));        %频率w2的范围

                        % ************************** 波束响应赋值 ************************** %
                        P=zeros(length(Omega1),length(Omega2));    %初始化波束形成矩阵

                        for i=1:length(Omega1)         %波束形成循环
                                for n=1:length(Omega2)
                                                if ((c*Omega1(i)/d/Omega)^2+(c*Omega2(n)/d/Omega)^2)<=1
                                                        v_x=(exp(-1j*(vec)*pi*(c*Omega1(i)/d/Omega))).';
                                                        v_y=(exp(-1j*(vec)*pi*(c*Omega2(n)/d/Omega))).';
                                                        v=kron(v_y,v_x);

                                                        P(i,n)=h1*v;

                                                else
                                                        P(i,n)=0;
                                                end                
                                end
                        end

                        % ************************** 对波束响应做逆离散傅里叶变换得到加权系数 ************************** %
                        B1=zeros(M,M);

                        for i=1:M
                            B1(i,:)=exp(-1j*(0:(M-1))*Omega1(i));    
                        end

                        B2=zeros(M,M);

                        for i=1:M
                            B2(:,i)=(exp(-1j*(0:(M-1))*Omega2(i))).';    
                        end

                        D1=B1\P/B2;

                        temp=[D1(((M-(array_numx-3)/2):M),:);D1((1:((array_numx-1)/2+1)),:)];
                        D=[temp(:,(M-(array_numy-3)/2):M),temp(:,(1:((array_numy-1)/2+1)))];

                        w=reshape(D,array_numx*array_numy,1);  %得到加权系数矩阵

                        % ************************** 波束形成 ************************** %
                        
                        vec_beam=((1-array_numx)/2):((array_numx-1)/2);

                        beam_FIB=zeros(length(theta),length(phi));
                        for k=1:length(theta)
                               for n=1:length(phi)
                                       v_x=(exp(-1j*vec_beam*Omega*d/c*sin(theta(k))*cos(phi(n))));
                                       v_y=(exp(-1j*vec_beam*Omega*d/c*sin(theta(k))*sin(phi(n))));
                                       v=kron(v_y,v_x);

                                       beam_FIB(k,n)=v*w;
                               end
                        end
                        beam_FIB_original=abs(beam_FIB);
                        beam_FIB=abs(beam_FIB)/max(max(abs(beam_FIB)));

                        theta_beam_plus60(OMEGA_num,:)=beam_FIB(:,find(phi==phi_point));  %俯仰角方向上的波束响应值
                        phi_beam_plus60(OMEGA_num,:)=beam_FIB(find(theta==theta_point(angle_num)),:);   %方位角方向上的波束响应值
                        
                        beam_theta_frequency_plus60(angle_num,OMEGA_num)=beam_FIB_original(find(theta==theta_point(angle_num)),find(phi==phi_point));
                                                %计算当俯仰角固定时不同方位角下各频率上的波束响应值

                        % ************************** 画图 ************************** %                 
                        figure(5); 
                        [x,y]=meshgrid(phi,theta); 
                        mesh(x*180/pi,y*180/pi,20*log10(beam_FIB));
                        hold on;
                        grid on;

                        ylabel('方位角/度');
                        xlabel('俯仰角/度');
                        zlabel('波束/dB');

                        title('基于离散傅里叶变换的FIB');

                        figure(6);
                        Patt3d(beam_FIB.',3);
                        hold on;
                        grid on;            

                end
        end

        figure(7);
        hold on; 
        box on;
        [x,y]=meshgrid(theta*180/pi,OMEGA_range/2/fh);  
        
        mesh(x,y,20*log10(theta_beam_plus60));
        
        xlabel('俯仰角/度');
        ylabel('数字频率');
        zlabel('波束/dB');
        
        title('不同频率对应的俯仰角上的波束响应');
        
%         figure(8); 
%         hold on; 
%         box on;
%         [x,y]=meshgrid(phi*180/pi,OMEGA_range/2/fh);
%         
%         mesh(x,y,20*log10(phi_beam_plus40));
%         
%         xlabel('方位角/度');
%         ylabel('数字频率');
%         zlabel('波束/dB');
%         
%         title('不同频率对应的方位角上的波束响应');