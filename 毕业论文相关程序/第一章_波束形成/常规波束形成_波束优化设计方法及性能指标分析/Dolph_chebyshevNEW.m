        
        %************************************************************%
        %             Dolph-Chebyshev��Ȩխ�������γɷ���             %
        %        copyright Liao Fengyi                               %
        %        last modified at 2015.11.24 ��Done��                %
        %************************************************************%

        clc;
        clear all

                     %% ************************��������*************************** %%
%         fh=12000;                       %���Ƶ��
%         fl=4000;                        %���Ƶ��
        fr=4000;                        %�ο�Ƶ��

        c=340;                          %����
        d=c/fr/4;                       %��Ԫ���Ϊ�ο�Ƶ�ʲ�����һ��
        lamda=c/fr;                     %�ο�Ƶ�ʶ�Ӧ�Ĳ���
        ang=-90:0.1:90;                 %��λ�Ƿ�ΧΪ-90�㵽90��

%         wide_c=lamda/8:lamda/8:lamda;
          wide_c=-10:-5:-50;                    %�԰꼶ȡֵ

%         wide_c=[15,20,25,30,35,40,45,50,55];        %����������ȡֵ   

%         wide_c=0:5:50;    %����ָ��ȡֵ

%         wide_c=5:5:30;                     %��Ԫ����ȡֵ
          
          
          N=10;
%           a=length(wide_c);
                     %% ************************Dolph-Chebyshev��Ȩ�������խ�������γ���************************ %%
for  n=1:1:N                     
        for coefficient_number=1:length(wide_c)       %ϵ��������Ϊcoefficient_number
%                 array_number=wide_c(coefficient_number);            %��Ԫ����
                
                array_number=15;                                    %��Ԫ����
                     %% ************************dolph-Chebyshev��Ȩ************************ %%        
                %**************ָ���԰꼶��ָ���԰꼶��������ֻ�ܶ�ѡһ��***************%
%                 sidelobe=-30;                                     %�԰꼶�����䲻Ϊ����ʱȡֵΪ-25dB
                sidelobe=wide_c(coefficient_number);                %�԰꼶ȡֵ
        
                R=10.^(-sidelobe/10); 
                z=cosh(1/(array_number-1)*acosh(sqrt(R)));          %������ʽ2-4��z�Ķ���
            
                %*************ָ��������************%
%                 mainlobe_wide=wide_c(coefficient_number);         %������ȡֵ            
%                 sin_theta_NN=sin(mainlobe_wide/2*pi/180);         
%         
%                 z=cos(pi/(2*(array_number-1)))/cos(pi*d/lamda*sin_theta_NN);

                    %% *************�����б�ѩ���Ȩϵ��************* %%
                wh=zeros(array_number,1);             %array_number���������Ԫ������ʽ2-4��M�ĺ��壬wh��ŵ����б�ѩ���Ȩ��ϵ������������Ԫ������ͬ

                wh(1)=(z^(array_number-1))/2;         %��m=1ʱ�ļ�Ȩֵ�ļ���

                for m=2:floor(array_number/2+1)       %floor������������ȡ������ĳһ��Ԫ�ص�С�ڵ�������������������ָ����2<=m<(M/2+1)�����
                        temp1=0;                      %temp1��Ϊһ����ʱ���������ڴ�ż��������Ȩֵ��
                        for i=1:(m-1)
                                temp1=temp1+0.5*(array_number-1)*factorial_0(m-2)*factorial_0(array_number-i-1)*(z^(array_number-2*m+1))*((z^2-1)^(m-i))/...
                                        factorial_0(m-i)/factorial_0(i-1)/factorial_0(m-i-1)/factorial_0(array_number-m);  
                        end                           %factorial_0����������һ�����Ľ׳ˣ���һ��ʵ�ֵ���ʽ2-4��2<=m<(M/2+1)�������Ȩֵ��ȡֵ��

                        wh(m)=temp1;
                end

                for m=ceil(array_number/2+1):array_number       %ceil����ȡ���Ǵ��ڵ��ڸ�������С����
                        wh(m)=wh(array_number+1-m);             %wh��Ȩֵ������һ���ԳƵľ�������Ԫ�ظ���Ϊ10��������wh(1)=wh(10),wh(2)=wh(9),......,wh(5)=wh(6),�Ǵ��м�ԳƵģ�
                                                               
                end                                    

                   %% *************�Ӳ���ָ��ļ�Ȩʸ��**************** %%
%               beam_point=wide_c(coefficient_number);
                beam_point=0;               %ָ������ָ��Ƕȣ�����Ƕ�Ҳ����Ϊ������
%                 d=wide_c(coefficient_number);

                v0=(exp(-1j*2*pi*fr*(0:(array_number-1))*d*sin(beam_point*pi/180)/c)).'/array_number;  %.'����˼�Ƿǹ���ת�ã�v0Ϊ���еĵ���ʸ��

                w=wh.*v0;                   %�Ӳ���ָ��ļ�Ȩʸ��

                  %% **************�����γ�*************** %%
                beam=zeros(1,length(ang));             %beam�ǲ�����Ӧ����

                for i=1:length(ang)                    %�����Ӧ��ÿһ���Ƕ��ϵĲ�����Ӧ�������棩ֵ��
                        v=(exp(-1j*2*pi*fr*(0:(array_number-1))*d*sin(ang(i)*pi/180)/c)).';
                        beam(i)=w'*v;             
                end

                beam_abs=abs(beam)/max(abs(beam));     %�Բ�����Ӧȡ����ֵ����һ��

                  %% **************ȷ������������************** %%
                MainLobe_down=find(diff(sign(diff(beam_abs)))>0)+1;  %diff�����������Ƕ�һ���������΢�֣���֣����㣻
                                                                     %sign�������������ж�һ��ֵ�������������������ֵΪ1����������ֵΪ-1,0�򷵻�0��
                                                                     %find���������������ڷ�������������������Ԫ�ص�����λ��(λ�õ��ж����ھ����У���һ�п�ʼ��
                                                                     %���϶��£�����Ϊ1��2��3...,Ȼ���ٴӵڶ��У�����������������)
                                                                     %�˴��������������ǣ�Ѱ�Ҳ�����Ӧ�İ��ݵ㼴����ͼ�йȵ׵�λ�ö�Ӧ�ĽǶ��±ꣻ
                temp2=ang(MainLobe_down);                            %�����temp2ָ�����ǲ���ͼ���ݵ�����Ӧ�ĽǶ�ֵ��
                temp3=temp2-beam_point*ones(1,length(temp2));        %temp3��temp2����Ķ��ǲ������ݵ�ĽǶ�ֵ��
                ang_temp=min(abs(temp3))*2;                          %�����ang_temp����Ӧ�ľ��ǲ�������Ŀ�ȼ���������������
                
                 %% **************ȷ���԰꼶************** %%
                temp4=find(ang<(beam_point-ang_temp/2));             %����Ѱ�ҵ��ǽǶ�С�ڲ�������������һ���ȵĽǶȵ��±�ֵ
                sidelobe_level=20*log10(max(beam_abs(temp4))) ;      %�����������԰�ˮƽ���԰꼶�����ֵ�ķֱ���ʾֵ��
                sidelobe_max=sidelobe_level;                         %sidelobe_maxָ�����԰�ˮƽ�����ֵ

                 %% **************ȷ���԰꼶����************** %%
                temp5=find((20*log10(beam_abs)-sidelobe_level*ones(size(beam_abs)))>0.0001);%�����0.001����Ϊ�趨��һ����ֵ����������������Ѱ�Ҳ�����Ӧֵ���԰꼶��ĽǶȵ��±ꣻ
                                                                                          
                temp6=ang(temp5);
                SL_wide=temp6(end)-temp6(1);                          %SL_wide��Ϊ�԰꼶����temp6��ʾ��������ˮƽ�����԰꼶ˮƽ�ĽǶȵķ�Χ��
                
                 %% **************��ͼ************** %%
                
%                 figure(1);grid on;
%                 hold on;
%                 
%                 plot(ang,20*log10(beam_abs));    %��������ͼ
%         
%                 plot(ang,sidelobe_level*ones(size(ang)),'m--');              %�����԰꼶ˮƽ��
%         
%                 plot(beam_point*ones(1,81),-80:0,'r');                     %��������ָ��ǵ�λ��
%         
%                 plot((beam_point-ang_temp/2).*ones(1,81),-80:0,'r--');       
%          
%                 plot((beam_point+ang_temp/2)*ones(1,81),-80:0,'r--');        %�������겨��������λ�ü�������������λ�ã�
%         
%                 axis([min(ang),max(ang),-80,0]);
%                 
%                 xlabel('��λ/��');
%                 ylabel('����/dB');
%                 title('��������');

                %**************y���ɡ�����ʾ**************%
                coff_y(coefficient_number,n)=ang_temp;               %coff_y��¼����ÿһ�β�ͬ���԰꼶����Ӧ�Ĳ�������������ˮƽ
                 
                %**************x����........��ʾ**************%
                                           % coff_x��¼����ÿһ�ε��԰꼶ˮƽ/��Ԫ����/����ָ���
                coff_x(coefficient_number,n)=sidelobe_level;         

%                 coff_x(coefficient_number,n)=array_number;
        
%                 coff_x(coefficient_number,n)=beam_point;
                
%                 coff_x(coefficient_number,n)=d/lamda;

        end
        figure(1);grid on;
                hold on;
                
                plot(ang,20*log10(beam_abs));    %��������ͼ
        
%                 plot(ang,sidelobe_level*ones(size(ang)),'m--');              %�����԰꼶ˮƽ��
        
                plot(beam_point*ones(1,81),-80:0,'r');                     %��������ָ��ǵ�λ��
        
                plot((beam_point-ang_temp/2).*ones(1,81),-80:0,'r--');       
         
                plot((beam_point+ang_temp/2)*ones(1,81),-80:0,'r--');        %�������겨��������λ�ü�������������λ�ã�
        
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
        % legend('���Ϊ����֮һ����','���Ϊ�ķ�֮һ����');
