        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                 Dolph-Chebyshev��Ȩ����                   %
        % copyright Liao Fengyi                                      %
        % last modified at 2015.09.22 ��finished��                             %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        clc;
        clear all

        %% -----------------��������-----------------------------------%
        fh=16000;                       %���Ƶ��
%         fl=4000;                        %���Ƶ��
        f0=16000;                          %�ο�Ƶ��

        c=340;                          %����
%         d=(c/fh)*2/3;                       %��Ԫ���Ϊ�ο�Ƶ�ʲ�����һ��
        lamda=c/fh;                     %�ο�Ƶ�ʶ�Ӧ�Ĳ���

        ang=-90:0.1:90;                 %��λ�Ƿ�ΧΪ-90�㵽90��


%         wide_c=[-20,-25,-30,-35,-40];         %�԰꼶ȡֵ
        wide_c=-30;        %�԰꼶ȡֵ

%         wide_c=[15,20,25,30,35,40,45,50,55];  %����������ȡֵ   

%         wide_c=[0,5,10,15,20,25,30,35,40,45,50];        %����ָ��ȡֵ

%         wide_c=[5,10,15,20,25];                       %��Ԫ����ȡֵ


        %% ----------------Dolph-Chebyshev��Ȩ�������խ�������γ���-----------%%
        D=lamda/10:lamda/5:5*lamda/2;
%         D=8*lamda/5;
        ang_temp=zeros(1,length(D));
        sl_FIB=zeros(1,length(D));
        SL_wide=zeros(1,length(D));
        Ang_3dB_SOCP_A=zeros(1,length(D));
    for d_num=1:length(D)
        d=D(d_num);
        for coff_num=1:length(wide_c)   %ϵ��������Ϊcoff_num
                array_num=11;           %��Ԫ����

                %% -------------dolph-Chebyshev��Ȩ--------------------------%%        
                %-------ָ���԰꼶��ָ���԰꼶��������ֻ�ܶ�ѡһ��-----------------------%
%                 sidelobe=-25;                              %�԰꼶�����䲻Ϊ����ʱȡֵΪ-25dB
                sidelobe=wide_c(coff_num);                %�԰꼶Ϊ����
        
                R=10.^(-sidelobe/10); 
                z=cosh(1/(array_num-1)*acosh(sqrt(R)));  %������ʽ2-4��z�Ķ���
            
                %----------------------ָ��������--------------------------------------%
%                 mainlobe_wide=wide_c(coff_num);          %������            
%                 sin_theta_NN=sin(mainlobe_wide/2*pi/180);       
%         
%                 z=cos(pi/(2*(array_num-1)))/cos(pi*d/lamda*sin_theta_NN);

                %---------------------�����б�ѩ���Ȩϵ��------------------------------------%
                wh=zeros(array_num,1);       %array_num���������Ԫ������ʽ2-4��M�ĺ���

                wh(1)=(z^(array_num-1))/2;

                for m=2:floor(array_num/2+1)  %floor������������ȡ������ĳһ��Ԫ�ص�С�ڵ��������������
                        temp=0;               %temp��Ϊһ����ʱ���������ڴ�ż��������Ȩֵ��
                        for i=1:(m-1)
                                temp=temp+0.5*(array_num-1)*factorial_0(m-2)*factorial_0(array_num-i-1)*(z^(array_num-2*m+1))*((z^2-1)^(m-i))/...
                                        factorial_0(m-i)/factorial_0(i-1)/factorial_0(m-i-1)/factorial_0(array_num-m);  
                        end             %factorial_0����������һ�����Ľ׳ˣ���һ��ʵ�ֵ���ʽ2-4��2<=m<(M/2+1)�������Ȩֵ��ȡֵ��

                        wh(m)=temp;
                end

                for m=ceil(array_num/2+1):array_num   %ceil����ȡ���Ǵ��ڵ��ڸ�������С����
                        wh(m)=wh(array_num+1-m);      %wh��Ȩֵ������һ���ԳƵľ�������Ԫ�ظ���Ϊ10����
                                                      %����wh(1)=wh(10),wh(2)=wh(9),......,wh(5)=wh(6),�Ǵ��м�ԳƵģ�
                end                                    

                %-----------------�Ӳ���ָ��ļ�Ȩʸ��---------------------------%
                beam_point=0;        %ָ������ָ��Ƕȣ�����Ƕ�Ҳ����Ϊ������

                v0=(exp(-1j*2*pi*f0*(0:(array_num-1))*d*sin(beam_point*pi/180)/c)).'/array_num;  %.'����˼�Ƿǹ���ת�ã�v0Ϊ����ʸ��

                w=wh.*v0;                           %�Ӳ���ָ��ļ�Ȩʸ��

                %% ----------------�����γ�-----------------------------%
                beam=zeros(1,length(ang));    %beam�ǲ�����Ӧ����

                for i=1:length(ang)
                        v=(exp(-1j*2*pi*f0*(0:(array_num-1))*d*sin(ang(i)*pi/180)/c)).';
                        beam(i)=w'*v;  %�����Ӧ��ÿһ���Ƕ��ϵĲ�����Ӧ�������棩ֵ��
                end

                beam_abs=abs(beam)/max(abs(beam));      %�Բ�����Ӧȡ����ֵ��
                                %% --------------��ͼ-------------------------------%%
                
                figure(1);
%                 grid on;
                hold on;
                
                plot(ang,20*log10(beam_abs));    %��������ͼ
                axis([-90,90,-80,5]);
%         
%                 plot(ang,-6*ones(size(ang)),'m--');              %�����԰꼶ˮƽ��
% %         
%                 plot(45*ones(1,121),-100:20,'r');                     %��������ָ��ǵ�λ��
%         
%                 plot((beam_point-ang_temp(d_num)/2).*ones(1,81),-80:0,'r--');       
%          
%                 plot((beam_point+ang_temp(d_num)/2)*ones(1,81),-80:0,'r--');        %�������겨��������λ�ü�������������λ�ã�
%         
%                 axis([min(ang),max(ang),-80,0]);
%                 
                xlabel('��λ��/(Degrees)');
                ylabel('��������/(dB)');
%                 title('��������');

                %----------------y���ɡ�����ʾ------------------------%
%                 coff_y(coff_num)=ang_temp(d_num);               %coff_y��¼����ÿһ�β�ͬ���԰꼶����Ӧ�Ĳ�������������ˮƽ
%                  
%                 %----------------x���ɡ�����ʾ------------------------%
%                 coff_x(coff_num)=sidelobe_level;         % coff_x��¼����ÿһ�ε��԰꼶ˮƽ

%                 coff_x(coff_num)=array_num;
        
%                 coff_x(coff_num)=beam_point;

                %------------------ȷ�ϲ���3dB����-----------------%

                IndMin=find((abs(beam_abs-((0.5)*beam_abs(901))))<0.02);
                temp=ang(IndMin);
                Ang_3dB=min(abs(temp-beam_point))*2;
                Ang_3dB_SOCP_A(1,d_num)=Ang_3dB;
               

                
                
                %--------------ȷ������������-------------------------%
%                 MainLobe_down=find(diff(sign(diff(beam_abs)))>0)+1;  %diff�����������Ƕ�һ���������΢�֣���֣����㣻
%                                                               %sign�������������ж�һ��ֵ�������������������ֵΪ1����������ֵΪ-1,0�򷵻�0��
%                                                               %find���������������ڷ�������������������Ԫ�ص�����λ��(λ�õ��ж����ھ����У���һ�п�ʼ��
%                                                               %���϶��£�����Ϊ1��2��3...,Ȼ���ٴӵڶ��У�����������������)
%                                                               %�˴��������������ǣ�Ѱ�Ҳ�����Ӧ�İ��ݵ㼴����ͼ�йȵ׵�λ�ö�Ӧ�ĽǶ��±ꣻ
%                 temp=ang(MainLobe_down);                      %�����tempָ�����ǲ���ͼ���ݵ�����Ӧ�ĽǶ�ֵ��
%                 temp=temp-beam_point*ones(1,length(temp));    
%                 ang_temp(d_num)=min(abs(temp))*2;                    %�����ang_temp����Ӧ�ľ��ǲ�������Ŀ�ȼ���������������
%                 
%                 %--------------ȷ���԰꼶-------------------------%
%                 temp=find(ang<(beam_point-ang_temp(d_num)/2));          %����Ѱ�ҵ��ǽǶ�С��-13.5��ĽǶȵ��±�ֵ
%                 sidelobe_level=20*log10(max(beam_abs(temp))) ;   %�����������԰�ˮƽ���԰꼶�����ֵ�ķֱ���ʾֵ��
%                 sl_FIB(d_num)=sidelobe_level;                           %sl_FIBָ�����԰�ˮƽ�����ֵ
% 
%                 %--------------ȷ���԰꼶����-------------------------%
%                 temp=find((20*log10(beam_abs)-sidelobe_level*ones(size(beam_abs)))>0.001);%�����0.001����Ϊ�趨��һ����ֵ��
%                                                                                           %��������������Ѱ�Ҳ�����Ӧֵ���԰꼶��ĽǶȵ��±ꣻ
%                 temp=ang(temp);
%                 SL_wide(d_num)=temp(end)-temp(1)     %SL_wide��Ϊ�԰꼶����temp��ʾ��������ˮƽ�����԰꼶ˮƽ�ĽǶȵķ�Χ��
                


        end
    end
        figure(2);
        plot(D/lamda,Ang_3dB_SOCP_A,'r*:');hold on; 
        xlabel('��Ԫ���/{\lambda}');
        ylabel('��������3dB���/(Degrees)');
%         grid on;
%         plot(D,SL_wide,'bx-.');hold on; grid on;
    
%         figure(2);hold on;box on;
% 
%         plot(coff_x,coff_y,'k-o');
%         grid on;
%         h = gca;
%         set(h,'FontSize',10,'FontName','����');
%         set(h,'FontName','Times New Roman');
% 
%         xlabel('�԰�ˮƽ(dB)');
% %         xlabel('����ָ���(��)');
% %         xlabel('��Ԫ����');
%         
%         ylabel('����������(��)');
%         axis([-40.1,-20,26,44]);
%         title('���������԰꼶�Ĺ�ϵ');
%         % title('����������Ԫ���Ĺ�ϵ');

%         % legend('SOCP','Dolph-chebyshev');
%         % legend('���Ϊ����֮һ����','���Ϊ�ķ�֮һ����');
