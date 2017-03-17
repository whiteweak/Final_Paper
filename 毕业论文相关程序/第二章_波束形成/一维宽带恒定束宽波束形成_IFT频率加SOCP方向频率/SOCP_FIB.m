%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             基于二阶锥规划（SOPC）的FIB                              %
% 期待波束的设计原理是稳健旁瓣控制高增益波束设计，参考频率为最低频率。    %
% FIB设计原理是以最小均方误差逼近于期望波束响应的波束。                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        clc;
        clear;
%         close all

        %% ----------------参数设置------------------------%%
        array_num=21;           %阵元数

        fl=10000;               %最低频率
%         fl=6000;               %最低频率
        fh=20000;               %最高频率
        fr=fl;                  %参考频率
        
        C=340;                  %声速
        array_space=C/2/fh;     %阵元间距
        
        sidelobe_level=-30;     %旁瓣级

        beam_point=0;           %波束指向
        mainbeam_wide=28;       %旁瓣级束宽

        silelobe=[-90:1:(beam_point-mainbeam_wide/2),(beam_point+mainbeam_wide/2):1:90];  %旁瓣离散化
        silelobe_num=length(silelobe);  %旁瓣离散数

       %% ------------期望波束响应设计-------------------%%
        matrix_zero=zeros(array_num*2,1);

        %--------------矩阵b------------------------%
        b=([-1,matrix_zero.']).';

        %---------------二阶锥约束的维数------------------%
        K.f=2;

        K.q=[3*ones(1,silelobe_num),(2*array_num+1)];

        %---------------矩阵C的转置--------------------------%
        c1=([1;0]).';

        c2=zeros(1,3*silelobe_num); 
        for i=1:silelobe_num
                c2(1,1+3*(i-1):3*i)=([10^(sidelobe_level/20);0;0]).';
        end

        c3=([0;matrix_zero]).';

        %---------------矩阵A--------------------------%
        v_0=exp(-1j*2*pi*fr*(0:(array_num-1))*array_space*sin(beam_point*pi/180)/C);
        va_0=[real(v_0),imag(v_0)].';
        vb_0=[-imag(v_0),real(v_0)].';

        a1=([0,va_0.';0,-vb_0.']).';

        a2=zeros((array_num*2+1),3*silelobe_num);

        for i=1:silelobe_num
                v_i=exp(-1j*2*pi*fr*(0:(array_num-1))*array_space*sin(silelobe(i)*pi/180)/C);
                va_i=[real(v_i),imag(v_i)].';
                vb_i=[-imag(v_i),real(v_i)].';

                a2(1:(array_num*2+1),1+3*(i-1):3*i)=([0,matrix_zero.';0,-va_i.';0,-vb_i.']).';
        end
        
        a3=([-1,matrix_zero.';matrix_zero,-eye(array_num*2,array_num*2)]).';
        % --------------问题的求解---------------------%
        C_col=([c1,c2,c3]).';

        A_col=([a1,a2,a3]).';

        [x,y,info]=sedumi(A_col,b,C_col,K);

        w=y(2:array_num+1)+1j*y((array_num+2):(array_num*2+1));
        %------------------波束形成---------------------%     
        ang=-90:1:90;   %方位角范围
        beam=zeros(1,length(ang));
        
        for i=1:length(ang)
                beam(i)=conj(exp(-1j*2*pi*fr*(0:(array_num-1))*array_space*sin(ang(i)*pi/180)/C))*w;
        end

        beam=(beam)/max(abs(beam));
        beam_abs=abs(beam);
        
        %--------------确认主瓣宽度---------------------------------------%
        IndMin=find(diff(sign(diff(beam_abs)))>0)+1;
        temp=ang(IndMin);
        temp=temp-beam_point*ones(1,length(temp));    
        ang_wide=min(abs(temp))*2;
        
        %-------------------确认旁瓣级-----------------------------------%
        temp=find(ang<(beam_point-ang_wide/2));
        sidelobe_level1=20*log10(max(beam_abs(temp)));  
        sidelobe_level1-sidelobe_level;
        
        %-------------------画图--------------------------------%       
        figure(1);hold on;box on;
        plot(ang,20*log10(beam_abs),'k-');

%         plot(beam_point*ones(1,81),-80:0,'k--');
%         plot(ang,sidelobe_level*ones(size(ang)),'k--');
                     
        h = gca;
        set(h,'FontSize',10,'FontName','宋体');
        set(h,'FontName','Times New Roman');
        
        axis([min(ang),max(ang),-80,0]);
        xlabel('Azimuth Angle(Degrees)');
        ylabel('Beam Pattern Gain(dB)');
        title('Desired Beam Pattern');

%         legend('情形1','情形2');

       %% --------------------频率不变波束形成设计--------------------%%
        %-----------------频率信息-------------------------%        

        f_index=linspace(fl,fh,25);
%         f_index=linspace(fl,fh,35);
        f_index_len=length(f_index);                                    %子带数

        %-----------------角度信息--------------------------%
        mainlobe=(beam_point-ang_wide/2):2:(beam_point+ang_wide/2);     %主瓣离散化
       
        mainlobe_num=length(mainlobe(1,:));                             %主瓣离散数

        %-----------------初始化--------------------------%
        beam_d=zeros(1,mainlobe_num);                                   %主瓣期望响应
        w_FIB=zeros(array_num,f_index_len);                             %加权向量

        beam_FIB_SOCP=zeros(f_index_len,length(ang));                        %频率不变波束响应
        beam_temp=zeros(1,length(ang));
                
        %-------------------主瓣期望响应-------------------------------%
        for i=1:length(mainlobe)
               beam_d(i)=conj(exp(-1j*2*pi*fr*(0:(array_num-1))*array_space*sin(mainlobe(i)*pi/180)/C))*w;
        end
        
        beam_d=beam_d/max(abs(beam_d));
        
        %% -------------二阶锥规划-------------------------%%                     
        %-----------------矩阵B-----------------------------------%
        bb=([-1,zeros(1,mainlobe_num),matrix_zero.']).';

        %-----------------二阶锥约束的维数--------------------------%
        Q.l=1;
        Q.q=[4*ones(1,mainlobe_num),3*ones(1,silelobe_num),(2*array_num+1)];  %考虑稳健性
%         Q.q=[4*ones(1,mainlobe_num),3*ones(1,silelobe_num)];  %不考虑稳健性

        %-----------------矩阵C的转置和矩阵A---------------------------%
        cc1=zeros(1,4*mainlobe_num);
        aa1=zeros((array_num*2+mainlobe_num+1),4*mainlobe_num);

        cc2=zeros(1,3*silelobe_num);
        aa2=zeros((array_num*2+mainlobe_num+1),3* silelobe_num);
        
        temp=9.472;       %白噪声增益损失
        coff=(1/array_num)*10^(temp/10);      %这里的coff指的是控制稳健性的gama的值
        cc3=([coff;matrix_zero]).';     %论文中为0.4217（对应的白噪声增益损失为9.472）；  这个值会影响主瓣均方根误差
        aa3=([zeros(1,mainlobe_num+1),matrix_zero.';zeros(2*array_num,mainlobe_num+1),-eye(array_num*2,array_num*2)]).';
       
        cc4=0;        
        aa4=([-1,ones(1,mainlobe_num),zeros(1,2*array_num)]).';
        
        for m=1:f_index_len
                 for i=1:mainlobe_num                          
                      cc1(1,1+4*(i-1):4*i)=([1;2*real(beam_d(i));2*imag(beam_d(i));-1]).';

                      v_i=exp(-1j*2*pi*f_index(m)*(0:(array_num-1))*array_space*sin(mainlobe(i)*pi/180)/C);
                      va_i=[real(v_i),imag(v_i)].';
                      vb_i=[-imag(v_i),real(v_i)].';

                      q_one=zeros(mainlobe_num,1);
                      q_one(i)=1;

                      aa1(1:(array_num*2+mainlobe_num+1),1+4*(i-1):4*i)=([0,-q_one.',matrix_zero.';0,zeros(1,mainlobe_num),2*va_i.';0,zeros(1,mainlobe_num),2*vb_i.';0,-q_one.',zeros(1,2*array_num)]).';
                 end


                 for i=1:silelobe_num
                       cc2(1,1+3*(i-1):3*i)=([10^((sidelobe_level)/20);0;0]).';

                       v_i=exp(-1j*2*pi*f_index(m)*(0:(array_num-1))*array_space*sin(silelobe(i)*pi/180)/C);
                       va_i=[real(v_i),imag(v_i)].';
                       vb_i=[-imag(v_i),real(v_i)].';

                       aa2(1:(array_num*2+mainlobe_num+1),1+3*(i-1):3*i)=([0,zeros(1,mainlobe_num),matrix_zero.';0,zeros(1,mainlobe_num),-va_i.';0,zeros(1,mainlobe_num),-vb_i.']).';
                 end

        %----------------问题的求解-----------------------------%
%         %---------------------考虑稳健性-----------------------------------%
                CC_col=([cc4,cc1,cc2,cc3]).';

                AA_col=([aa4,aa1,aa2,aa3]).';
       %---------------------不考虑稳健性-----------------------------------%
%                 CC_col=([cc4,cc1,cc2]).';
% 
%                 AA_col=([aa4,aa1,aa2]).';

                [x,y,info]=sedumi(AA_col,bb,CC_col,Q);

                w_FIB(:,m)=y(mainlobe_num+2:mainlobe_num+array_num+1)+1j*y((mainlobe_num+array_num+2):(mainlobe_num+array_num*2+1));
                
                ww=w_FIB(:,m);
                
        %---------------波束形成--------------------------------%     
                for i=1:length(ang)
                      beam_temp(i)=conj(exp(-1j*2*pi*f_index(m)*(0:(array_num-1))*array_space*sin(ang(i)*pi/180)/C))*ww;
                end
                
                beam_FIB_SOCP(m,:)=beam_temp;     
               
        end
%         beam_FIB_SOCP=beam_FIB_SOCP.*beam_FIB_SOCP;
        beam_FIB_SOCP_norm=beam_FIB_SOCP/max(max(abs(beam_FIB_SOCP)));
        beam_FIB_SOCP_abs_norm=abs(beam_FIB_SOCP_norm);
        
        %%------------------画图----------------------------------%%             
        figure(2);
        [x,y]=meshgrid(ang,f_index);
        
        mesh(x,y/fh*pi,20*log10(beam_FIB_SOCP_abs_norm));

        axis([-90,90,pi*fl/fh,pi*fh/fh,-80,0]);
    
        title('Frequency Invariant Beamforming Based On SOCP');
        zlabel('Beam Pattern Gain(dB)');
        ylabel('Digital Frequency(pi)');
        xlabel('Azimuth Angle(Degrees)');
        
        figure(3);hold on;box on;       
        
        sl_FIB=-80;
        temp=find(ang<(beam_point-ang_wide/2));
        
        %------------------确认主瓣所在区域----------------------------------%
        errors=zeros(1,f_index_len);
        ang_mainlobe=(beam_point-ang_wide/2):1:(beam_point+ang_wide/2);
        ang_index=zeros(1,length(ang_mainlobe));
        
        for i=1:length(ang_mainlobe)
                ang_index(i)=find(ang==ang_mainlobe(i));
        end
        
                
        for i=1:f_index_len  
                %---------------------确认旁瓣级---------------------%
                temp1=20*log10(max(beam_FIB_SOCP_abs_norm(i,temp)));    %旁瓣级               
                if temp1>sl_FIB
                      sl_FIB=temp1;  
                end
  
                %--------------------求均方误差和---------------------%
                errors(i)=sum(abs(beam_FIB_SOCP_norm(i,ang_index)-beam(ang_index)).^2);                
                
        end
        
        %-----------------不同频率分量对应的主瓣均方根误差------------------------------%
        plot(f_index/fh*pi,10*log10(errors)-10*log10(length(ang_mainlobe))*ones(size(errors)),'k:o');%+,*,。,x,s（方块),d（菱形）,p（五角形）
        
        h = gca;
        set(h,'FontSize',10,'FontName','宋体');
        set(h,'FontName','Times New Roman');

        title('Mean Square Errors of Mainlobe Associated With Different Frequencies');
        ylabel('Mean Square Errors of Mainlobe(dB)');
        xlabel('Digital Frequency(pi)');
%         axis([0.6,3.2,-100,-60]);
%         legend('情形1','情形2');
%         legend('最低数字频率pi/2','最低数字频率pi/3');
%         legend('主瓣宽度10度','主瓣宽度30度','主瓣宽度40度','主瓣宽度46度');
 
