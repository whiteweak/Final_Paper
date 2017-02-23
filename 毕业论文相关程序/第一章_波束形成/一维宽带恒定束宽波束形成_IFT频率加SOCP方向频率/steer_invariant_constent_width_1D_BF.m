       

        %% ****************************************************** %%
         %                 ���򲻱�㶨����խ��һά�����γ�  
         %  �˳����ǻ��ڶ���׶�滮��SOCP���ķ��򲻱�㶨����խ��һά�����γɣ������õķ�
         %  ��������26������ģ����Ļ���˼���ǣ�
         %    �ڱ�֤����ָ��������Ʋ�������������Ӧ�Ĳο�����������֮��ľ��������С��
         %    ����£�ʹ������ƵĲ������԰�ˮƽ��ͣ��Ӷ�ʵ���ڸ��������Ͼ��к㶨����
         %    ���������Ҫ���խ�������γ����Ĳ��������Ȼ��沨��ָ��ǵ�������𽥱��
         %    �����⣬�ܹ�ʹ��һ������ָ��Ƿ�Χ�ڵĲ��������ȱ��ֺ㶨���۽���һά����
         %    �γɼ�ָ���Ϊ��λ�ǡ�
         %      
        %% ****************************************************** %%         

        clc;
        clear;
%         close all

        %% ******************* �������� ******************* %%
        array_num=13;                   %������Ԫ��

        f0=20000;                       %�������Ƶ�ʼ��ο�Ƶ�ʣ���������ֻ��Է��򣬹ʲ�����Ƶ�ʲ�����
        C=340;                          %����Ϊ340m/s
        array_space=C/2/f0;             %��Ԫ���ȡΪ�ο�Ƶ�ʶ�Ӧ������һ�룬��֤�����������
        
        sidelobe_level=-35;             %�ο������԰꼶

        ang_reference=0;               %���òο�����ָ���
        mainbeam_width=25;              %���òο����������ȣ������õ����԰꼶����

        sidelobe=[-90:1:(ang_reference-mainbeam_width/2),(ang_reference+mainbeam_width/2):1:90];  
                                        %ȷ���԰�����ȡֵ��Χ��������ɢ��

        silelobe_num=length(sidelobe);  %�԰���ɢ��֮��ĸ���

        ang=-90:0.1:90;                   %���÷�λ�ǵķ�Χ

        beam=zeros(1,length(ang));      %��ʼ��������Ӧ����

        %% ******************* ����SOCP��Ʋο�������Ӧ ******************* %%
        
        matrix_zero=zeros(array_num*2,1);           %����һ����ʼ26*1ά�����

        % ******************* ���þ���b ******************* %
        
        b=([-1,matrix_zero.']).';        %���þ���b��Ϊһ��27*1ά�ľ���

        % ******************* ���ö���׶Լ����ά�� ******************* %
        K.f=2;                          
        K.q=[3*ones(1,silelobe_num),(2*array_num+1)];

        % ******************* ����C��ת�� ******************* %
        
        c1=([1;0]).';                  %c1Ϊ1*2ά 
        c2=zeros(1,3*silelobe_num);    %c2Ϊ1*456ά����ʼ��Ϊ0���� 
        for i=1:silelobe_num
                c2(1,1+3*(i-1):3*i)=([10^(sidelobe_level/20);0;0]).';
        end                            %��һ���c2������и�ֵ������ÿ����ֵ���ظ�

        c3=([0;matrix_zero]).';        %c3Ϊ1*27�������

        % ******************* �������A ******************* %
        
        a3=([-1,matrix_zero.';matrix_zero,-eye(array_num*2,array_num*2)]).'; %a3Ϊ27*27�ĸ���λ����
                                            
        v_0=exp(-1j*2*pi*f0*(0:(array_num-1))*array_space*sin(ang_reference*pi/180)/C);
                                             %V_0��ʾ�ڲο�Ƶ���²ο������ϵ����еĵ���ʸ����Ϊ1*13ά
        va_0=[real(v_0),imag(v_0)].';        %va_0ΪV_0��ʵ�����鲿�����ɵľ���Ϊ26*1ά
        
        vb_0=[-imag(v_0),real(v_0)].';       %vb_0ΪV_0��ʵ�����鲿�����ɵľ���Ϊ26*1ά����va_0��Щ΢���� 

        a1=([0,va_0.';0,-vb_0.']).';         %a1��ά��Ϊ27*2��a1��Ӧ�ڲο������ϵĲ������γ� 

        a2=zeros((array_num*2+1),3*silelobe_num);     %a2Ϊ27*456ά��0���� 

        for i=1:silelobe_num
                v_i=exp(-1j*2*pi*f0*(0:(array_num-1))*array_space*sin(sidelobe(i)*pi/180)/C);
                va_i=[real(v_i),imag(v_i)].';
                vb_i=[-imag(v_i),real(v_i)].';

                a2(1:(array_num*2+1),1+3*(i-1):3*i)=([0,matrix_zero.';0,-va_i.';0,-vb_i.']).';
        end                                         %a2Ϊ��Ӧ���԰���ɢ������ÿһ����λ�ǶȵĲ������γ� 

        % ******************* �������� ******************* %
                      % ����ǰ�渳ֵ�õ���ص�a,b��c�����Լ�sedumi����������sedumi�����䣨����������Ȩֵ����� 
        
        C_col=([c1,c2,c3]).';       %����C_col����

        A_col=([a1,a2,a3]).';       %����A_col����
        
        [x,y,info]=sedumi(A_col,b,C_col,K);     %����sedumi�����䣨��������������Ķ���׶�滮��͹�Ż������⣬
                                                %�ο������ļ�Ȩֵ����y�У������yΪ27*1ά

        wr=y(2:array_num+1)+1j*y((array_num+2):(array_num*2+1));    %�õ��ο������ϵĲ�����Ȩֵ���Ǵ�y��ȡ����
        
        % ******************* �����γ� ******************* %
        
        for i=1:length(ang)
                beam(i)=conj(exp(-1j*2*pi*f0*(0:(array_num-1))*array_space*sin(ang(i)*pi/180)/C))*wr;
        end                                     %������ÿһ����ɢ�Ƕ��ϵĲ�����Ӧֵ

        beam_abs=abs(beam)/max(abs(beam));      %����õĲ�����Ӧ���������ֵ�������������ֵ���й�һ�����Է��㻭ͼ

        % ******************* ��ͼ ******************* %   
        
        figure(1);hold on;box on;
        
        plot(ang,20*log10(beam_abs),'r:');

        plot(ang_reference*ones(1,81),-80:0,'g');
        plot(ang,sidelobe_level*ones(size(ang)),'k--');
        
        axis([min(ang),max(ang),-80,0]);

        xlabel('��λ��(��)');
        ylabel('����(dB)');
 %         title('��������');   
 
        h = gca;
        set(h,'FontSize',10,'FontName','����');
        set(h,'FontName','Times New Roman');
        
%         legend('40��','60��')

       %% ******************* ���򲻱䲨���γ���� ******************* %%
        
        % ******************* ��ʼ����Ӧ��25������ָ��ǵļ�Ȩʸ���Ͳ�����Ӧ����
        ang_range=-90:5:90;                                %��������ָ��ǵķ�Χ��ÿ5��ȡһ��ֵ
        multi_w=zeros(array_num,length(ang_range));        %��ͬ����ָ��Ƕ�Ӧ�ļ�Ȩʸ������Ϊ13*25ά
        beam_multi=zeros(length(ang_range),length(ang));   %��ͬ����ָ��Ƕ�Ӧ�Ĳ�����Ӧ����Ϊ25*181ά
        MainWidth=zeros(1,length(ang_range));
        % ******************* �Ƕ���Ϣ ******************* %
        
        for ang_num=1:length(ang_range)           %ѭ����1��25,��ѭ��ʵ�ֵ��Ƕ���ѡȡ��-60��60�ȷ�Χ��
                                                  %ÿһ������ָ����ϵĲ����γɣ�������Ȼ��SOCP
                beam_point=ang_range(ang_num);    %����ָ��Ƿֱ�ȡ-60��5��60�㣬��25������ָ���
                
                mainlobe=(beam_point-mainbeam_width/2):1:(beam_point+mainbeam_width/2);  
                                                  %��������1�Ƚ�����ɢ�����Ա��ں��潫ƽ�ƺ�Ĳ����������һһ�Ա�                
                sidelobe=[-90:1:(beam_point-mainbeam_width/2),(beam_point+mainbeam_width/2):1:90];  
                                                  %��ÿһ����������Ӧ���԰���1�Ƚ�����ɢ��

                mainlobe_num=length(mainlobe(1,:));   %ȷ��������ɢ��
                silelobe_num=length(sidelobe);        %ȷ���԰���ɢ��

                % ******************* ����������Ӧ ******************* %
                beam_d=zeros(1,mainlobe_num);               %����������Ӧ
               
                for i=1:length(mainlobe)   %���ѭ���������������ڸò���ָ���ÿһ��������ɢ�Ƕ��ϵĲ�����Ӧֵ
                                           %������Ϊ֮��Լ����������д����
                       beam_d(i)=conj(exp(-1j*2*pi*f0*(0:(array_num-1))*array_space*sin((mainlobe(i)-(beam_point-ang_reference))*pi/180)/C))*wr;
                end
                beam_d=beam_d/max(abs(beam_d));    %�Ըò���ָ����ϵĲ�����Ӧ���й�һ��

        %% ******************* ����׶�滮������Ըò���ָ��ǵļ�Ȩֵ ******************* %%
                % ******************* ���Ʋ���ָ����ϵĲ�����ӦΪ1 ******************* %
                cc1=([1;0]).';                      %ȷ��cc1
                 
                v_0=exp(-1j*2*pi*f0*(0:(array_num-1))*array_space*sin(beam_point*pi/180)/C);
                                                    %v_0��ʾ�ò���ָ��Ƕ��ϵĵ���ʸ��  
                va_0=[real(v_0),imag(v_0)].';
                vb_0=[-imag(v_0),real(v_0)].';      %va_0��vb_0����v_0��ʵ�����鲿����Ϲ���
                
                aa1=([0,zeros(1,mainlobe_num),va_0.';0,zeros(1,mainlobe_num),-vb_0.']).';   %����aa1����
                
                % ******************* �������������С��ĳ���� ******************* %
                cc2=10;       
                aa2=([0,ones(1,mainlobe_num),zeros(1,2*array_num)]).';       %����aa2����
               
                % ******************* �������� ******************* %
                cc3=zeros(1,4*mainlobe_num);         
                aa3=zeros((array_num*2+mainlobe_num+1),4*mainlobe_num);

                % ******************* �����԰� ******************* %
                cc4=zeros(1,3*silelobe_num);
                aa4=zeros((array_num*2+mainlobe_num+1),3* silelobe_num);
                
                % ******************* �����Ƚ��� ******************* %
                cc5=([0.4;matrix_zero]).';
                aa5=([zeros(1,mainlobe_num+1),matrix_zero.';zeros(2*array_num,mainlobe_num+1),-eye(array_num*2,array_num*2)]).';   
                        
                % ******************* ����B ******************* %
                bb=([-1,zeros(1,mainlobe_num),matrix_zero.']).';

                % ******************* ����׶Լ����ά�� ******************* %
                Q.f=2;
                Q.l=1;
                Q.q=[4*ones(1,mainlobe_num),3*ones(1,silelobe_num)];

                           
                % ******************* ����C��ת�ú;���A ******************* %

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

                % ******************* �������� ******************* %
                CC_col=([cc1,cc2,cc3,cc4]).';

                AA_col=([aa1,aa2,aa3,aa4]).';

                [x,y,info]=sedumi(AA_col,bb,CC_col,Q);

                ww=y(mainlobe_num+2:mainlobe_num+array_num+1)+1j*y((mainlobe_num+array_num+2):(mainlobe_num+array_num*2+1));
                      
                multi_w(:,ang_num)=ww;
                
                % ******************* �����γ� ******************* %     
                for i=1:length(ang)
                      beam(i)=conj(exp(-1j*2*pi*f0*(0:(array_num-1))*array_space*sin(ang(i)*pi/180)/C))*ww;
                end
                
                beam_abs=abs(beam)/max(abs(beam));  %��һ��
        
                beam_multi(ang_num,:)=beam_abs;
               
                %-----------------ȷ�����������----------------------------%
                IndMin=find(diff(sign(diff(beam_abs)))>0)+1;
                temp=ang(IndMin);
                temp=temp-beam_point*ones(1,length(temp));    
                mainlobe_wide=min(abs(temp))*2;
                MainWidth(1,ang_num)=mainlobe_wide;
                
                %----------------------ȷ���԰꼶----------------------------%
%                 temp=find(ang>(beam_point+mainlobe_wide/2));
                temp=find(ang<(beam_point-mainlobe_wide/2));
                sidelobe=20*log10(max(beam_abs(temp)));   
     
               
        end
        %% -------------------��ͼ--------------------------------%   
        figure(2);hold on;box on;

        plot(ang,20*log10(beam_abs),'k:');

        plot(beam_point*ones(1,81),-80:0,'k--');
        plot(ang,sidelobe*ones(size(ang)),'k--');

        axis([min(ang),max(ang),-80,0]);
        
        xlabel('��λ��(��)');
        ylabel('����(dB)');
        title('�㶨������');
        
        h = gca;
        set(h,'FontSize',10,'FontName','����');
        set(h,'FontName','Times New Roman');
        
%         legend('40��','60��')
        

        figure(3);
        [x,y]=meshgrid(ang,ang_range);
        mesh(y,x,20*log10(beam_multi));
        
        xlabel('����ָ���(��)');
        ylabel('��λ��(��)');
        zlabel('����(dB)');
        
        h = gca;
        set(h,'FontSize',10,'FontName','����');
        set(h,'FontName','Times New Roman');   
        
        figure(4);
        plot(ang_range,MainWidth,'r*-');
        axis([-65,65,20,40]);
        
        dlmwrite('multicoff.txt',multi_w,';');   %�豣���Ȩʸ��ʱ�Ų���
         
