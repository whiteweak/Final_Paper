        
    % ************************** ���ڶ�ά��ɢ����Ҷ�任��FIB **************************%  
 
 
        clc;
        clear all
%         close all

        %% ************************** ����������ã��������ϲ�������һ�£�************************** %%
        array_number=21;           %������Ԫ����Ϊ21��
        fl=10000;                   %���Ƶ�� 
%         fl=6000;                   %���Ƶ�� 
        fh=20000;                  %���Ƶ��
        fr=fl;                     %�ο�Ƶ��

        c=340;                     %����
        d=c/fr/2;                  %��Ԫ���ȡ���Ƶ�ʶ�Ӧ������һ�룬��������ص��������������
        lamda=c/fr;                %�ο�Ƶ�ʶ�Ӧ�Ĳ���

        ang=-90:1:90;              %��λ�Ƿ�ΧΪ-90��- 90��
        
        

       %% **************************Dolph-Chebyshev��Ȩ�������������Ӧ **************************%
        
        array_num=13;                          %�������������Ӧ���õĲο���Ԫ��
        
        vec=-(array_num-1)/2:(array_num-1)/2;                        %����������Ӧʱ�ĵ���ʸ�������м���ԪΪ�ο���������ģ���Ϊ����11������vecȡΪ-5:5
        
       %% *************** ָ���԰꼶ˮƽֵ����������ֵ *************** %%
        
        % *************** ָ���԰꼶 *************** %
        
        sidelobe=-30;                         %�԰꼶ˮƽ����Ϊ-30dB
        R=10.^(-sidelobe/10);            
        z=cosh(1/(array_num-1)*acosh(sqrt(R)));     %����������ʽ2-4�ж����zֵ

        % *************** ָ�������� *************** %
%         mainbeam_wide=24;                     %������
%         sin_theta_NN=sin(mainbeam_wide/2*pi/180);       
%         z=cos(pi/(2*(array_num-1)))/cos(pi*d/lamda*sin_theta_NN);

        % *************** �����б�ѩ���Ȩϵ�� *************** %
        
        wh=zeros(array_num,1);              %��ʼ���б�ѩ���Ȩֵ����Ϊ0����

        wh(1)=(z^(array_num-1))/2;          %���㵱m=1ʱ��Ȩֵ

        for m=2:floor(array_num/2+1)        %floor������������ȡ������ĳһ��Ԫ�ص�С�ڵ�������������������ָ����2<=m<(M/2+1)�����
                                            %temp1��Ϊһ����ʱ���������ڴ�ż��������Ȩֵ��
              temp=0;                       %tempΪ��ʱ����
              for i=1:(m-1)
                   temp=temp+0.5*(array_num-1)*factorial_0(m-2)*factorial_0(array_num-i-1)*(z^(array_num-2*m+1))*((z^2-1)^(m-i))/...
                      factorial_0(m-i)/factorial_0(i-1)/factorial_0(m-i-1)/factorial_0(array_num-m);  
              end                           %factorial_0����������һ�����Ľ׳ˣ���һ��ʵ�ֵ���ʽ2-4��2<=m<(M/2+1)�������Ȩֵ��ȡֵ��
              
              wh(m)=temp;                   %��Ȩֵ����ֵ
        end

        for m=ceil(array_num/2+1):array_num       
              wh(m)=wh(array_num+1-m);      %ceil����ȡ���Ǵ��ڵ��ڸ�������С��������һ��ʵ�ֵ���m>=M/2+1����µ�Ȩֵ�ļ���   
        end


        % *************** �����б�ѩ���Ȩʸ�� *************** %
        
        beam_point=0;            %������λ����Ϊ0��

        v0=(exp(-1j*2*pi*fr*vec*d*sin(beam_point*pi/180)/c)).'/array_num;     %���沨���γɵļ�Ȩʸ��

        wd=wh.*v0;               %�����б�ѩ���Ȩ��ĵ���ʸ��

        % *************** ����������Ӧ *************** %
        
        hh=wd.';                              %�����Ȩ����ʸ���ķǹ���ת�ã�����֮������Ӧ�ļ���

        beam_temp=zeros(1,length(ang));       %�Լ�Ȩ������Ӧ�����ʼ��

        for i=1:length(ang)                   %�����б�ѩ���Ȩ������ÿһ����λ���ϵĲ�����Ӧ
                beam_temp(i)=exp(-1j*vec*pi*sin(ang(i)*pi/180))*(hh');
        end
        
        beam=(beam_temp)/max(abs(beam_temp));          %�Բ�����Ӧ���й�һ��
        beam_temp=abs(beam);                           %��ʱ��beam_tempΪ��������Ӧ���й�һ����ȡ����ֵ֮��Ĳ�����Ӧ����

        % *************** ȷ�������� *************** %
        IndMin=find(diff(sign(diff(beam_temp)))>0)+1;  %diff�����������Ƕ�һ���������΢�֣���֣����㣻
                                                       %sign�������������ж�һ��ֵ�������������������ֵΪ1����������ֵΪ-1,0�򷵻�0��
                                                       %find���������������ڷ�������������������Ԫ�ص�����λ��(λ�õ��ж����ھ����У���һ�п�ʼ��
                                                       %���϶��£�����Ϊ1��2��3...,Ȼ���ٴӵڶ��У�����������������)
                                                       %�˴��������������ǣ�Ѱ�Ҳ�����Ӧ�İ��ݵ㼴����ͼ�йȵ׵�λ�ö�Ӧ�ĽǶ��±ꣻ
                                                       
        temp=ang(IndMin);                              %�����tempָ�����ǲ���ͼ���ݵ�����Ӧ�ĽǶ�ֵ��
        temp=temp-beam_point*ones(1,length(temp));     %temp��temp����Ķ��ǲ������ݵ�ĽǶ�ֵ��   
        
        temp=sort(abs(temp));                          %�Բ������ݵ�Ƕ�ֵ�ľ���ֵ���д�С���������  
        
        mainbeam1=temp(1);                             
        mainbeam2=temp(2);                             %mainbeam1��mainbeam2������ǲ���������ֵ���Ƕȣ���һ��
        
        mainbeam_wide=mainbeam1+mainbeam2;             %mainbeam_wide��ʾ���ǲ��������ȵ�ֵ
        
        % *************** ȷ���԰꼶 *************** %
        
        temp=find(ang<(beam_point-mainbeam_wide/2));   %�����temp��������԰�����ĽǶ�ֵ����Ӧ�Ĳ�����Ӧ�е��±�ֵ
        sidelobe_level=20*log10(max(beam_temp(temp))); %������԰�ˮƽָ�����԰���������ֵ��dBֵ������ȷ���԰�ֵ�Ƿ�Ϊ�趨��-30dB
        
        
     

        %% *************** �������б�ѩ���Ȩ��Ƴ�������������Ӧ�Ĳ���ͼ *************** %%
        
        figure(1);hold on;box on;
        plot(ang,20*log10(beam_temp),'k-');
                
%         plot(beam_point*ones(1,81),-80:0,'k');
%         plot(ang,sidelobe_level*ones(size(ang)),'k--');

        axis([min(ang),max(ang),-80,0]);
        
        h = gca;
        set(h,'FontSize',10,'FontName','Times New Roman');
        set(h,'FontName','Times New Roman');
        
        xlabel('Azimuth Angle(Degrees)');
        ylabel('Beam Pattern Gain(dB)');
        title('Desired Beam Pattern');

%         legend('0��','20��','40��','60��');
%         legend('�������Ƶ��pi/2','�������Ƶ��pi/3');
%         legend('�ο���Ԫ��9','�ο���Ԫ��11');
%         legend('SR ����','SOCP ����','IFT ����');
%         legend('SOCP ����','IFT ����');
%         legend('����1','����2');

        
        %% *************** Ƶ�ʲ��䲨���γ� *************** %%
        Omega_range=2*pi*linspace(fl,fh,25);                    %�Ӵ�����,�õ�w�ķ�ΧΪ15��խ���Ӵ�
%         Omega_range=2*pi*linspace(fl,fh,35); 

%         Omega_range=2*pi*[linspace(fl,6000,4),linspace(6000,fh,12)];     %�Ӵ�����

        M=64;                                           %w1����ɢ��������Ϊ��IDFT�ĵ���
        Omega1_range=-pi*ones(1,M)+2*pi/M*(0:(M-1));    %w1�ķ�ΧΪ-pi����M-1/M��*pi����pi/MΪ���������64��ֵ��

        % *************** ������Ӧ��ֵ *************** %
        P=zeros(length(Omega1_range),length(Omega_range));% P��ŵ�������Ĳ�����Ӧ��ά����w1��w�Ĵ�С����

        for i=1:length(Omega1_range)      %��ѭ����������ͨ�������滻������ÿһ��Ƶ�ʷ�����ÿһ���Ӵ�
                                          %����Ӧ�Ĳ�����Ӧ��ͬʱ����A��w1��Ϊ0�����γɵĲ���ͼ���ӹ⻬
                                          %��Ӧ�������Ϲ�ʽ3-29��
                for k=1:length(Omega_range)      %ѭ����ΧΪ1��w�ĳ���
                        if abs(Omega1_range(i))<(Omega_range(k)*d/c)     
                                temp=(Omega1_range(i)*c/Omega_range(k)/d);  
                                P(i,k)=exp(-1j*vec*pi*temp)*(hh');
                                %����������ʱ����(c*w1)/(w*d)���벨����Ӧ���̼�����˱����滻,�����õ���Ȩ������Ӧ
                        else
                                 P(i,k)=0;     %���������A��Omega1����Ϊ0����ʹ��P������ƽ����
                        end
                end
        end


        % *************** �Բ�����Ӧ����ɢ����Ҷ��任�õ���Ȩϵ�� *************** %
        B=zeros(M,M);   %MΪw1����ɢ������Ϊ��IDFT�ĵ������빫ʽ��N~���Ӧ

        for i=1:M
            B(i,:)=exp(-1j*(0:(M-1))*Omega1_range(i));    
        end

        % temp=abs(sum(diag(B)))/M;
        % B=B+temp/1000*eye(M,M);

        D=B\P;   %��B����ȥ���P���󣬵õ��Ľ����inv��B��*P��P�������B�������

        w=[D(((M-(array_number-3)/2):M),:);D((1:((array_number-1)/2+1)),:)];

        %-----------------�����γ�-----------------------------------%
        beam_FIB=zeros(length(Omega_range),length(ang));

        vec1=(-(array_number-1)/2):(array_number-1)/2;

        for i=1:length(Omega_range)
                for k=1:length(ang)                              
                        beam_FIB(i,k)=exp(-1j*vec1*Omega_range(i)*sin(ang(k)*pi/180)*d/c)*w(:,i);
                end
        end
        
        beam_FIB=beam_FIB/max(max(abs(beam_FIB)));
        beam_FIB_a=abs(beam_FIB);
         
        %---------------ȷ��������������--------------------------%
        mainlobe=(beam_point-mainbeam1):1:(beam_point+mainbeam2); 
        ang_index=zeros(1,length(mainlobe));

        for i=1:length(mainlobe)
               ang_index(i)=find(ang==mainlobe(i));
        end
        %%-------------------��ͼ---------------------------------------%%
        figure(2);grid on;
        [x,y]=meshgrid(ang,Omega_range);

        mesh(x,y/2/fh,20*log10(beam_FIB_a));

        axis([min(ang),max(ang),fl/fh*pi,pi,-80,0]);

        title('Frequency Invariant Beamforming Based On IFT');
        zlabel('Beam Pattern Gain(dB)');
        ylabel('Digital Frequency(pi)');    
        xlabel('Azimuth Angle(Degrees)');
        
        h = gca;
        set(h,'FontSize',10,'FontName','����');
        set(h,'FontName','Times New Roman');

        figure(3);hold on;box on;
        errors=zeros(1,length(Omega_range));  
        sl_FIB=-100;
        
        for i=1:length(Omega_range)
                %----------------ȷ���԰꼶----------------------------%
                temp=find(ang<(beam_point-mainbeam_wide/2-3));
                
                temp1=20*log10(max(beam_FIB_a(i,temp)));                  
                if temp1>sl_FIB
                      sl_FIB=temp1;  
                end
                
                %---------------���������---------------------------%
                errors(i)=sum(abs(beam_FIB(i,ang_index)-beam(1,ang_index)).^2);

                
        end
        
        %-----------------��ͬƵ�ʷ�����Ӧ��������������------------------------------%      
        plot(Omega_range/2/fh,10*log10(errors)-10*log10(length(mainlobe))*ones(size(errors)),'k--o');%+,*,��,x,s������),d�����Σ�,p������Σ�
        
        h = gca;
        set(h,'FontSize',10,'FontName','����');
        set(h,'FontName','Times New Roman');
        
        title('Mean Square Errors of Mainlobe Associated With Different Frequencies');
        ylabel('Mean Square Errors of Mainlobe(dB)');
        xlabel('Digital Frequency(pi)');
        axis([1.4,3.2,-90,-60]);
%         axis([0.6,3.2,-90,-50]);
%         legend('SR ����','SOCP ����','IFT ����');
%         legend('����1','����2');
%         legend('0��','20��','40��');
%         legend('�������Ƶ��pi/2','�������Ƶ��pi/3');
%         legend('�ο���Ԫ��9','�ο���Ԫ��11');
%         legend('SOCP ����','IFT ����');
     

