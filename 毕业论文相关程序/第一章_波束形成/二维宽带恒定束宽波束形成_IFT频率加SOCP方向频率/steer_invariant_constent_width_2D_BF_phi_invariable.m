        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                 方向不变恒定束宽窄带二维波束形成             %
        %        用以实现基于矩形平面阵的二维方向不变恒定束宽波束形成   %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        clc;
        clear all;
        % close all
        %% ****************** 参数设置 ****************** %%
        array_numx=15;  %平面阵X轴上的阵元个数
        array_numy=15;  %平面阵Y轴上的阵元个数

        f0=20000;       %参考频率

        c=340;          %传播速度
        d=c/f0/2;       %阵元间距
        lamda=c/f0;     %最小波长

        theta=(0:0.2:90);          %俯仰角范围
        phi=(-180:0.2:180);        %方位角范围
        %% ****************** 基于SOCP优化加权法设计参考波束响应 ****************** %%
        
        num=15;                  %设计参考波束响应对应的阵元数
        vec=-(num-1)/2:(num-1)/2;
        
        wh=dlmread('coff_w_socp.txt');   %读入利用SOCP方法设计出来的权值
        wh_x=wh.';                       %x轴向上的加权值；
        wh_y=wh.';                       %y轴向上的加权值； 
        
        % ****************** 参考波束响应 ******************%
        theta_p=46;  %俯仰角参考波束指向
        phi_p=0;     %方位角参考波束指向

        v_x=conj(wh_x.*exp(-1j*(vec)*pi*sin(theta_p*pi/180)*cos(phi_p*pi/180)));  %参考俯仰角上x轴向的导向矢量；
        v_y=conj(wh_y.*exp(-1j*(vec)*pi*sin(theta_p*pi/180)*sin(phi_p*pi/180)));  %参考俯仰角上y轴向的导向矢量；

        h1=kron(v_y,v_x);     %含波束指向的加权矢量，为两个方向上导向矢量的克罗内克积  
     
        beam_temp=zeros(length(theta),length(phi));   %初始化参考波束响应矩阵

        for k=1:length(theta)         %此循环的作用是对俯仰角和方位角观察角度范围内的所有角度进行波束形成，即进行扫描
                for i=1:length(phi)
                        v_x=(exp(-1j*(vec)*pi*sin(theta(k)*pi/180)*cos(phi(i)*pi/180))).';
                        v_y=(exp(-1j*(vec)*pi*sin(theta(k)*pi/180)*sin(phi(i)*pi/180))).';
                        v=kron(v_y,v_x);

                        beam_temp(k,i)=(h1*v);     %计算波束响应
                end
        end
        beam_temp=abs(beam_temp)/max(max(abs(beam_temp)));      %对波束响应进行取绝对值和归一化处理

        %% ****************** 画球坐标系下的波束图 ****************** %%
        
         figure(1);
         Patt3d(beam_temp.',3);     %画出球坐标下的波束图
         hold on;
         grid on;   
        
         
        %% ****************** 画直角坐标系下的波束图 ****************** %%
        
        figure(2);    
        [x,y]=meshgrid(phi,theta);       %对参数进行分割形成网格，便于画图         
        mesh(x,y,20*log10(beam_temp));   %画出直角坐标系下的波束图
        
        xlabel('方位角/度');
        ylabel('俯仰角/度');
        zlabel('波束/dB');
        h = gca;
        set(h,'FontSize',10.5,'FontName','宋体');
        set(h,'FontName','Times New Roman');  
        title('参考波束图');
         
        %% ****************** 方位角方向上的波束图，俯仰角为波束指向角 ****************** %%
        
        figure(3);
        subplot(2,1,1);
        plot(phi,20*log10(beam_temp(find(theta==theta_p),:)),'r');   %画出方位角方向上，参考俯仰角方向的波束图
        hold on;
        grid on;
        
        beam_abs=beam_temp(find(theta==theta_p),:);    %寻找当俯仰角确定时的各个方位角上的波束响应值
        
        % ****************** 确认两零点间束宽 ****************** %
        
        IndMin=find(diff(sign(diff(beam_abs)))>0)+1;
        temp=phi(IndMin);
        temp=temp-phi_p*ones(1,length(temp));    
        phi_wide1=min(abs(temp))*2;                  %确定方位角方向上两零点间束宽

        % ****************** 确认旁瓣级 ****************** %
        
        temp=find(phi<(phi_p-phi_wide1/2));
        sidelobe1=20*log10(max(beam_abs(temp)));     %确定方位角方向上最高旁瓣水平即旁瓣级
        % sidelobe1=sidelobe;
        % sl_phi=-25;

        % ****************** 确认旁瓣级束宽 ****************** %
        
        temp=find((20*log10(beam_abs)-sidelobe1*ones(size(beam_abs)))>0.1); 
        SL_phiwide1=phi(temp(end));
        SL_phiwide2=phi(temp(1));
        SL_phiwide3=SL_phiwide1-SL_phiwide2;   %计算方位角方向上旁瓣级束宽      

        xlabel('方位角(度)');
        ylabel('波束(dB)');
        % title('参考波束图（方位角）');
        % legend('(30,0)','(50,0)')

        %% ****************** 俯仰角方向上的波束图，方位角为波束指向角 ****************** %
        
        subplot(2,1,2);
        plot(theta,20*log10(beam_temp(:,find(phi==phi_p))),'b');  %画出俯仰角方向上，参考方位角方向的波束图
        hold on;
        grid on;

        beam_abs=beam_temp(:,find(phi==phi_p));   %寻找当方位角确定时的各个俯仰角上的波束响应值
        
        % ****************** 确认两零点间束宽 ****************** %
        
        IndMin=find(diff(sign(diff(beam_abs)))>0)+1;
        temp=theta(IndMin);
        temp=temp-theta_p*ones(1,length(temp));    
        theta_wide1=min(abs(temp))*2;                  %确定俯仰角方向上两零点间束宽

         % ****************** 确认旁瓣级 ****************** %
        temp=find(theta>(theta_p+theta_wide1/2));
        sidelobe2=20*log10(max(beam_abs(temp)));       %确定俯仰角方向上最高旁瓣水平即旁瓣级
        % sidelobe2=sidelobe;
        % sl_theta=-25;

         %------------确认旁瓣级束宽----------------------------%
        temp=find((20*log10(beam_abs)-sidelobe2*ones(size(beam_abs)))>0.2); 
        SL_thetawide1=theta(temp(end));
        SL_thetawide2=theta(temp(1));
        SL_thetawide3=SL_thetawide1-SL_thetawide2      %计算俯仰角方向上旁瓣级束宽

        xlabel('俯仰角(度)');
        ylabel('波束(dB)');
        % title('参考波束图（俯仰角）');

        sidelobe=max(sidelobe1,sidelobe2);

        %% ****************** 基于SOCP设计不同波束指向下的加权矢量 ****************** %%
        
        phi_range=0;                        %方位角扫描范围，这里设为固定值45°，减少计算量便于画图
        theta_range=20;                 %俯仰角扫描范围
        w_coff=zeros(num*num,length(phi_range)*length(theta_range));      %初始化权值矩阵
        
        theta_beam_mainlobewidth=zeros(1,length(theta_range));
        beam_theta_plus60=zeros(length(theta_range),length(theta));
                 %此程序中由于俯仰角固定，故存储的是俯仰角固定为46度时，方位角方向上不同角度的波束响应值用于画图和计算主瓣宽度
                 
        for beamnum_phi=1:length(phi_range)
                for beamnum_theta=1:length(theta_range)
                        phi_point=phi_range(beamnum_phi);     %方位角波束指向
                        theta_point=theta_range(beamnum_theta);                   %俯仰角波束指向

                        temp1=SL_phiwide3;        %方位角方向上的旁瓣级束宽
                        temp2=SL_thetawide3;      %俯仰角方向上的旁瓣级束宽

                        mainlobe_theta=(theta_point-temp2/2):2:(theta_point+temp2/2);  %主瓣离散化
                        mainlobe_phi=(phi_point-temp1/2):5:(phi_point+temp1/2);  %主瓣离散化

                        silelobe_theta=[0:2:(theta_point-temp2/2),(theta_point+temp2/2):2:90];  %旁瓣离散化
                        silelobe_phi=[-180:10:(phi_point-temp1/2),(phi_point+temp1/2):10:180];  %旁瓣离散化

                        slphi_num=length(silelobe_phi);  
                        sltheta_num=length(silelobe_theta);
                        silelobe_num=slphi_num*sltheta_num;   %旁瓣离散数

                        mlphi_num=length(mainlobe_phi);  
                        mltheta_num=length(mainlobe_theta);
                        mainlobe_num=mlphi_num*mltheta_num;  %旁瓣离散数

                % ****************** 初始化 ****************** %
                        matrix_zero=zeros(num*num*2,1);
                 
                        beam_d=zeros(mltheta_num,mlphi_num);               %参考波束响应
                
                        phi_temp=mainlobe_phi-(phi_point-phi_p)*ones(size(mainlobe_phi));         %参考波束对应的方位角范围
                        theta_temp=mainlobe_theta-(theta_point-theta_p)*ones(size(mainlobe_theta));   %参考波束对应的俯仰角范围

                % ****************** 主瓣期望响应 ****************** %
                        for k=1:mltheta_num
                                for i=1:mlphi_num
                                        v_x=(exp(-1j*(vec)*pi*sin(theta_temp(k)*pi/180)*cos(phi_temp(i)*pi/180))).';
                                        v_y=(exp(-1j*(vec)*pi*sin(theta_temp(k)*pi/180)*sin(phi_temp(i)*pi/180))).';
                                        v=kron(v_y,v_x);  

                                        beam_d(k,i)=(h1*v);                          
                                end
                        end
                
                        beam_d=beam_d/max(max(abs(beam_d)));  %归一化


              %% ****************** 二阶锥规划设计方向不变恒定束宽窄带二维波束响应 ****************** %%
                
                % ****************** 控制波束指向的响应为1 ****************** %
                        cc1=([1;0]).';

                        v_x=(exp(-1j*(vec)*pi*sin(theta_point*pi/180)*cos(phi_point*pi/180))).';
                        v_y=(exp(-1j*(vec)*pi*sin(theta_point*pi/180)*sin(phi_point*pi/180))).';
                        v_0=kron(v_y,v_x);  

                        va_0=[real(v_0);imag(v_0)];
                        vb_0=[-imag(v_0);real(v_0)];
                        aa1=([0,zeros(1,mainlobe_num),va_0.';0,zeros(1,mainlobe_num),-vb_0.']).';

                % ****************** 控制主瓣的误差和小于某个数 ****************** %
                        cc2=1;  %0.4
                        aa2=([0,ones(1,mainlobe_num),zeros(1,2*num*num)]).';

                % ****************** 控制主瓣宽度 ****************** %
                        cc3=zeros(1,4*mainlobe_num);
                        aa3=zeros((num*num*2+mainlobe_num+1),4*mainlobe_num);

                % ****************** 控制旁瓣水平 ****************** %
                        cc4=zeros(1,3*silelobe_num);
                        aa4=zeros((num*num*2+mainlobe_num+1),3* silelobe_num);

                % ****************** 控制波束稳健性 ****************** %
                        cc5=([0.1;matrix_zero]).';
                        aa5=([zeros(1,mainlobe_num+1),matrix_zero.';zeros(2*num*num,mainlobe_num+1),-eye(num*num*2,num*num*2)]).';  
                
                % ****************** 构建矩阵B ****************** %
                        bb=([-1,zeros(1,mainlobe_num),matrix_zero.']).';

                % ****************** 二阶锥约束的维数 ****************** %
                        Q.f=2;
                        Q.l=1;
                        Q.q=[4*ones(1,mainlobe_num),3*ones(1,silelobe_num),num*num*2+1];

                % ****************** 矩阵C的转置和矩阵A ****************** %
                        for k=1:mltheta_num
                                for i=1:mlphi_num
                                        temp=(((k-1)*mlphi_num+(i-1))*4+1):(((k-1)*mlphi_num+i)*4);
                                        cc3(1,temp)=([1;2*real(beam_d(k,i));2*imag(beam_d(k,i));-1]).';   

                                        v_x=(exp(-1j*(vec)*pi*sin(mainlobe_theta(k)*pi/180)*cos(mainlobe_phi(i)*pi/180))).';
                                        v_y=(exp(-1j*(vec)*pi*sin(mainlobe_theta(k)*pi/180)*sin(mainlobe_phi(i)*pi/180))).';
                                        v_i=kron(v_y,v_x);  

                                        va_i=[real(v_i);imag(v_i)];
                                        vb_i=[-imag(v_i);real(v_i)];

                                        q_one=zeros(mainlobe_num,1);
                                        temp1=(k-1)*mlphi_num+(i-1)+1;
                                        q_one(temp1)=1;   

                                        aa3(:,temp)=([0,-q_one.',matrix_zero.';0,zeros(1,mainlobe_num),2*va_i.';0,zeros(1,mainlobe_num),2*vb_i.';0,-q_one.',zeros(1,2*num*num)]).';
                                end
                        end


                        for k=1:sltheta_num
                                for i=1:slphi_num

                                        temp=(((k-1)*slphi_num+(i-1))*3+1):(((k-1)*slphi_num+i)*3);
                                        cc4(1,temp)=([0;0;0]).';
     
                                        v_x=(exp(-1j*(vec)*pi*sin(silelobe_theta(k)*pi/180)*cos(silelobe_phi(i)*pi/180))).';
                                        v_y=(exp(-1j*(vec)*pi*sin(silelobe_theta(k)*pi/180)*sin(silelobe_phi(i)*pi/180))).';
                                        v_i=kron(v_y,v_x);  

                                        va_i=[real(v_i);imag(v_i)];
                                        vb_i=[-imag(v_i);real(v_i)];

                                        aa4(:,temp)=([-1,zeros(1,mainlobe_num),matrix_zero.';0,zeros(1,mainlobe_num),-va_i.';0,zeros(1,mainlobe_num),-vb_i.']).';
                                end
                        end

                % ****************** 问题的求解 ****************** %
                        CC_col=([cc1,cc2,cc3,cc4,cc5]).';

                        AA_col=([aa1,aa2,aa3,aa4,aa5]).';

                        [x,y,info]=sedumi(AA_col,bb,CC_col,Q);
                        y(1)

                        w=y(mainlobe_num+2:mainlobe_num+num*num+1)+1j*y((mainlobe_num+num*num+2):(mainlobe_num+num*num*2+1));
                        w=w/max(abs(w));

                % ****************** 方向不变恒定束宽波束形成 ****************** %               
                        for k=1:length(theta)
                                for i=1:length(phi)
                                        v_x=(exp(-1j*(vec)*pi*sin(theta(k)*pi/180)*cos(phi(i)*pi/180))).';
                                        v_y=(exp(-1j*(vec)*pi*sin(theta(k)*pi/180)*sin(phi(i)*pi/180))).';
                                        v=kron(v_y,v_x);


                                        beam_temp(k,i)=v'*w;
                                end
                        end

                        beam_temp=abs(beam_temp)/max(max(abs(beam_temp)));      %归一化

                        w_coff(:,(beamnum_phi-1)*length(theta_range)+beamnum_theta)=w;
                        
                        beam_theta_plus60(beamnum_theta,:)=beam_temp(:,find(phi==phi_point));
                        
                        %% ****************** 画出方位角固定时波束主瓣宽度随俯仰角变化曲线 ****************** %%
%                           
%                 % ****************** 确认两零点间束宽 ****************** %
%                 
%                         IndMin=find(diff(sign(diff(beam_theta_plus50(beamnum_theta,:))))>0)+1;
%                         temp=theta(IndMin);
%                         temp=temp-theta_range(beamnum_theta)*ones(1,length(temp));    
%                         theta_wide1=min(abs(temp))*2;                  %确定方位角方向上两零点间束宽
% 
%                 % ****************** 确认旁瓣级 ****************** %
%         
%                         temp=find(theta<(theta_range(beamnum_theta)-theta_wide1/2));
%                         sidelobe1=20*log10(max(beam_theta_plus50(temp)));     %确定方位角方向上最高旁瓣水平即旁瓣级
%                         % sidelobe1=sidelobe;
%                         % sl_phi=-25;
% 
%                 % ****************** 确认旁瓣级束宽 ****************** %
%         
%                         temp=find((20*log10(beam_theta_plus50(beamnum_theta,:))-sidelobe1*ones(1,length(theta)))>0.1); 
%                         SL_thetawide1=theta(temp(end));
%                         SL_thetawide2=theta(temp(1));
%                         SL_thetawide3=SL_thetawide1-SL_thetawide2;   
%                                     %计算方位角方向上旁瓣级束宽,SL_phiwide3即为在固定俯仰角时，各个方位角上的波束主瓣旁瓣级束宽值  
%                                     
%                           %% ****************** 存储计算出来的波束主瓣宽度值（旁瓣级束宽或者波束两零点间束宽） ****************** %%            
% %                         phi_beam_mainlobewidth(1,beamnum_phi)=SL_phiwide3;  
%                                     %此时phi_beam_sidelobewidth矩阵存储的是方位角方向上主瓣旁瓣级束宽值，用于后续画图
%                         theta_beam_mainlobewidth(1,beamnum_theta)=theta_wide1;  
%                                     %此时phi_beam_sidelobewidth矩阵存储的是方位角方向上主瓣波束两零点间束宽值，用于后续画图


                %% ****************** 画图 ****************** %
                
%                         figure(4);   %画出球坐标系下的波束图
%                         Patt3d(beam_temp.',3);

                        figure(5);    %画出直角坐标系下的波束图
                        [x,y]=meshgrid(phi,theta);
                        mesh(x,y,20*log10(beam_temp));

                        xlabel('方位角(度)');
                        ylabel('俯仰角(度)');
                        zlabel('波束(dB）');
                
                        h = gca;
                        set(h,'FontSize',10,'FontName','宋体');
                        set(h,'FontName','Times New Roman'); 
                
                        title('基于SOCP的方向不变波束形成');
              
                        figure(6);
                % ****************** 方位角方向上的波束图，俯仰角为参考波束指向俯仰角 ****************** %
%                         subplot(2,1,1);
                        plot(phi,20*log10(beam_temp(find(theta==theta_point),:)),'k');hold on;grid on;

                        h = gca;
                        set(h,'FontSize',10,'FontName','宋体');
                        set(h,'FontName','Times New Roman');  

                        xlabel('方位角(度)');
                        ylabel('波束(dB）');
                        axis([-180,180,-80,0]);
        
                % ****************** 俯仰角方向上的波束图，方位角为参考波束指向方位角 ****************** %
                        figure(7);
                        plot(theta,20*log10(beam_temp(:,find(phi==phi_point))),'k');hold on;grid on;

                        h = gca;
                        set(h,'FontSize',10,'FontName','宋体');
                        set(h,'FontName','Times New Roman'); 
                
                        xlabel('俯仰角(度)');
                        ylabel('波束(dB）');
%                         axis([0,90,-50,0]);
               

                end
        end

%                 dlmwrite('coeffcient_phi0_theta60_all.txt',w_coff,';');       %将权值输出为.txt文件，便与后续的调用与修改

%                 figure(8);
%                 theta=(30:2:66);
%                 plot(theta,theta_beam_mainlobewidth(1,:),'r*');grid on;
