        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %    ���ڶ�ά����Ҷ��任�Ŀ��Ƶ��-���򲻱�㶨�����ά�����γ�             %
        %    �������Ի���SOCP�Ķ�ά���򲻱�㶨�������γ�����Ƴ��Ĳ���            %
        %    ��Ϊ����������ɢ����Ҷ��任��IFT���ĺ㶨����Ƶ�ʲ��䲨����             %
        %    �ɵĲο���������������ʵ�ֿ������Ƶ�ʲ��䲨���γɡ�
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        clc;
        clear all
        % close all

        %% ************************** �������� ************************** %%
        
        array_numx=15;                  %ƽ����X���ϵ���Ԫ����
        array_numy=15;                  %ƽ����Y���ϵ���Ԫ����

        fh=20000;                       %�ο�Ƶ��
        fl=fh/2;                      %���Ƶ��
        f0=fh;                          %���Ƶ��

        c=340;                          %�����ٶ�
        d=c/f0/2;                       %��Ԫ���
        lamda=c/f0;                     %�ο�Ƶ�ʶ�Ӧ�Ĳ���

        theta=pi/180*(0:0.2:90);          %������ɨ�跶Χ
        phi=pi/180*(-180:0.2:180);        %��λ��ɨ�跶Χ

        num=15;                         %�����������õ���Ԫ��Ϊnum*num
        vec=-(num-1)/2:(num-1)/2;


               %% ************************** Ƶ�ʲ��䲨���γ� ************************** %%
        theta_point=pi/180*60;                       %�����ǲ���ָ��60��
        
        phi_point=pi/180*0;              %��λ�ǲ���ָ��Χ

        multicoff=dlmread('coeffcient_phi0_theta60_all.txt');             %coff_w.txt�ļ��д�ŵ�����steer_invariant_constent_width_2D_BF������
                                                                   %�Ѿ�����õ�����������Ӧ��Ȩֵ��������Ϊ�˼������
                                                                   
               %% ************************** ��ʼ���ο�����ͬƵ�ʲ�����Ӧ�����Լ������Ⱦ��� ************************** %%                                             
        
        OMEGA_range=2*pi*(linspace(fl,fh,20));    %�����Ӵ���Χ

        theta_beam_plus60=zeros(length(OMEGA_range),length(theta));   %�����Ƿ����ϵĲ�����Ӧֵ
        phi_beam_plus60=zeros(length(OMEGA_range),length(phi));       %��λ�Ƿ����ϵĲ�����Ӧֵ
        
        beam_theta_frequency_plus60=zeros(length(theta_point),length(OMEGA_range));
                           %��ʼ��������Ϊ�ο�����ʱ��ͬƵ�ʶ�Ӧ�Ĳ�����Ӧ����Ϊ31x20ά����¼����Ƶ�ʱ仯ʱ�õ��Ĳ���ͼ��ÿһ����Ƶ����λ�ڲο�
                           %�������򣨷�λ��0�㣬������46�㣩λ�õĲ�����Ӧֵ
%         beam_theta_frequency_plus40_mainlobewidth=zeros(length(theta_point),length(OMEGA_range));
                           %��¼beam_phi_frequency��ÿ��������������ֵ�����������ͼ
        

        for angle_num=1:length(theta_point)        
                h1=multicoff(:,angle_num);           %ȡ�þ����ÿһ��
                h1=h1';                              %Ϊ*15x15��*1��������

               % ************************** ��������������Ӧ ************************** %
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
                      %�����γɵĲ����Ѿ�ʵ���˻����б�ѩ���Ȩ��խ����ά�ο���������Լ�ʵ���˸ò�������SOCP�����ķ��򲻱��ԣ�
                      %������Ϊ�������Ļ��ڸ���Ҷ��任��Ƶ�ʲ��䲨���γɵĲο���������ʵ�����ÿһ�������ϵ�Ƶ�ʲ��䲨���γ�
                      %�ļ�Ȩֵ�����
                                                                          

                %% ************************** ֱ������ϵ�Ĳ���ͼ ************************** %
                
                [x,y]=meshgrid(phi*180/pi,theta*180/pi);
                figure(1);           
                mesh(x,y,20*log10(beam_temp));hold on;grid on;

                xlabel('��λ��/��');
                ylabel('������/��');
                zlabel('����/dB');
                title('��������ͼ');

                %% ************************** ������ϵ�Ĳ���ͼ ************************** %
                figure(2);
                Patt3d(beam_temp.',3);hold on;grid on;        

%                 %% ************************** ��λ�Ƿ����ϵĲ���ͼ��������Ϊ����ָ��� ************************** %
%                 figure(3);       
%                 plot(phi*180/pi,20*log10(beam_temp(find(theta==theta_point),:)),'k:');hold on;box on;
%               
%                 axis([-180,180,-80,0]);
%         
%                 h = gca;
%                 set(h,'FontSize',10,'FontName','����');
%                 set(h,'FontName','Times New Roman');
%                 
%                 xlabel('��λ��(��)');
%                 ylabel('����(dB)');
%         %         title('��λ�Ƿ����ϵ���������ͼ');
%         
%         %         legend('(30,0)','(30,120)');

                %% ************************** �����Ƿ����ϵĲ���ͼ����λ��Ϊ����ָ��� ************************** %
                figure(4);
                plot(theta*180/pi,20*log10(beam_temp(:,find(phi==phi_point))),'k:');hold on;box on;
                
                axis([0,90,-80,0]);
                
                h = gca;
                set(h,'FontSize',10,'FontName','����');
                set(h,'FontName','Times New Roman');
                
                xlabel('������(��)');
                ylabel('����(dB)');
        %         title('�����Ƿ����ϵ���������ͼ');
        
        %         legend('(30,0)','(30,120)');


               %% ************************** ���ڶ�ά����Ҷ��任��FIB ************************** %%
               
               
                for OMEGA_num=1:length(OMEGA_range) 

                        Omega=OMEGA_range(OMEGA_num);     

                        M=64;                                         %Ƶ����ɢ����
                        Omega1=-pi*ones(1,M)+2*pi/M*(0:(M-1));        %Ƶ��w1�ķ�Χ
                        Omega2=-pi*ones(1,M)+2*pi/M*(0:(M-1));        %Ƶ��w2�ķ�Χ

                        % ************************** ������Ӧ��ֵ ************************** %
                        P=zeros(length(Omega1),length(Omega2));    %��ʼ�������γɾ���

                        for i=1:length(Omega1)         %�����γ�ѭ��
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

                        % ************************** �Բ�����Ӧ������ɢ����Ҷ�任�õ���Ȩϵ�� ************************** %
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

                        w=reshape(D,array_numx*array_numy,1);  %�õ���Ȩϵ������

                        % ************************** �����γ� ************************** %
                        
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

                        theta_beam_plus60(OMEGA_num,:)=beam_FIB(:,find(phi==phi_point));  %�����Ƿ����ϵĲ�����Ӧֵ
                        phi_beam_plus60(OMEGA_num,:)=beam_FIB(find(theta==theta_point(angle_num)),:);   %��λ�Ƿ����ϵĲ�����Ӧֵ
                        
                        beam_theta_frequency_plus60(angle_num,OMEGA_num)=beam_FIB_original(find(theta==theta_point(angle_num)),find(phi==phi_point));
                                                %���㵱�����ǹ̶�ʱ��ͬ��λ���¸�Ƶ���ϵĲ�����Ӧֵ

                        % ************************** ��ͼ ************************** %                 
                        figure(5); 
                        [x,y]=meshgrid(phi,theta); 
                        mesh(x*180/pi,y*180/pi,20*log10(beam_FIB));
                        hold on;
                        grid on;

                        ylabel('��λ��/��');
                        xlabel('������/��');
                        zlabel('����/dB');

                        title('������ɢ����Ҷ�任��FIB');

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
        
        xlabel('������/��');
        ylabel('����Ƶ��');
        zlabel('����/dB');
        
        title('��ͬƵ�ʶ�Ӧ�ĸ������ϵĲ�����Ӧ');
        
%         figure(8); 
%         hold on; 
%         box on;
%         [x,y]=meshgrid(phi*180/pi,OMEGA_range/2/fh);
%         
%         mesh(x,y,20*log10(phi_beam_plus40));
%         
%         xlabel('��λ��/��');
%         ylabel('����Ƶ��');
%         zlabel('����/dB');
%         
%         title('��ͬƵ�ʶ�Ӧ�ķ�λ���ϵĲ�����Ӧ');