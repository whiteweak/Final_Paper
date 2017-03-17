
       %*************************************************************%
       %          ���ڶ���׶�滮��SOPC���Ĳ����Ż����               %
       %          ���������ԭ�����Ƚ��԰���Ƹ����沨�����         %
       %   copyright Liao Fengyi                                     %
       %   last modified at 2015.11.24 (done)                        %

        clc;
        clear;

                     %% ****************��������**************** %%
%         fl=4000;                                %���Ƶ��
%         fh=12000;                               %���Ƶ��
        fr=20000;                                %�ο�Ƶ��
        C=340;                                  %����

%         wide_c=0:5:20;                                          %����ָ��Ƕ�ȡֵ
%         wide_c=8:1:15;                                          %��Ԫ����ȡֵ
%         wide_c=-14:-1:-32;               %�԰�ˮƽȡֵ
        wide_c=-40;               %�԰�ˮƽȡֵ

%         mainlobe_c=20.2:0.8:34.6;                      %�벨��ʱ��Ӧ���԰꼶����
        mainlobe_c=25;                      %�벨��ʱ��Ӧ���԰꼶����
        
%         mainlobe_c=37:2.8:65;                     %�ķ�֮һ����ʱ��Ӧ���԰꼶����
                     %% ****************�����Ż����**************** %% 
        for coefficient_number=1:length(wide_c)               %coefficient_number����ֵΪ1��5��Ϊ�԰�ˮƽ��ȡֵ��
            sidelobe_level=wide_c(coefficient_number);        %�԰꼶ˮƽȡֵΪ����
%             sidelobe_level=-25;
            
            array_number=15;                                  %��Ԫ�� 
%             array_space=C/4/fr; 
            array_space=C/2/fr;                               %��Ԫ���Ϊ�ķ�֮һ����������ΪC/fr

            beam_point=0;
%             beam_point=wide_c(coefficient_number);                                    %����ָ��Ƕ�����Ϊ0��

%             mainbeam_wide=29;
            mainbeam_wide=mainlobe_c(coefficient_number);     %�԰꼶����ȡֵ��������

            silelobe=[-90:1:(beam_point-mainbeam_wide/2),(beam_point+mainbeam_wide/2):1:90];  %��1�ȵļ����԰�������нǶ��ϵ���ɢ��(��Է�λ��)
 
            silelobe_number=length(silelobe);        %�԰���ɢ�������Ŀ�������sidelobe_num�������SOCP�е�Nsl

            ang=-90:0.1:90;                          %�趨��λ�ǵı仯��ȡֵ����Χ

            beam=zeros(1,length(ang));               %�����������󣬾���ά��Ϊ1x361�����������ÿһ����ɢ�Ƕȣ�0.5�ȣ���˵��

            matrix_zero=zeros(array_number*2,1);     %����һ����2����Ԫ������x1��0���󣬶�Ӧ���ǹ�ʽ�е������
            
            
            
                     %% ***************������SeDuMi������������SOCP��Ȩֵ�������ز�����Ϊ�����㷨�ĺ��Ĳ���***************** %%
                     
            %*************����b*************%
            b=([-1,matrix_zero.']).';                %����һ����2����Ԫ����+1��x1�ľ���������������Ϊ��ʽ�е�b����

            %*************����׶Լ����ά��*************%
            K.f=2;                                   % K.fָ����sedumi�����У�K�ṹ�������Ķ�Ӧ�ڡ���Լ���ġ���ά��������
                                                     % ��ֵ��ʾ��������Լ���������м�����ʽԼ��

            K.q=[3*ones(1,silelobe_number),(2*array_number+1)];      % K.qȷ������sedumi�����У�������Լ����ά�ȣ����ֵ���빫ʽ��

            %*************����C��ת��*************%
            c1=([1;0]).';                            %��.'�������������þ���ķǹ���ת��,c1��ά��Ϊ1x2��

            c2=zeros(1,3*silelobe_number);           %����һ��c2���󣬳���ʼԪ��Ϊ0
            for i=1:silelobe_number                  %��ѭ���������Ƕ�c2������Ӧ�ĸ�ֵ����
                    c2(1,1+3*(i-1):3*i)=([10^(sidelobe_level/20);0;0]).';   %��c2����ֵ,c2����ָ�����ǹ�ʽ��C2+i��i=0,1,2,...Nsl������      
            end
                
          
            c3=([0;matrix_zero]).';                 %��c3����ֵ

            %*************����A*************%
            a3=([-1,matrix_zero.';matrix_zero,-eye(array_number*2,array_number*2)]).';       % eye����һ��22x22�ĵ�λ�󣬱����ִ���������һ����׼�ĸ��ĵ�λ��a3��

            v_0=exp(-1j*2*pi*fr*(0:(array_number-1))*array_space*sin(beam_point*pi/180)/C);  %array_spaceΪ��Ԫ��࣬�൱��d��v_0Ϊ����ָ���Ϊ0��ʱ�ĵ���ʸ������
                                                                                             %�˴�Ϊһ1x11��ȫ1����
            va_0=[real(v_0),imag(v_0)].';      
            vb_0=[-imag(v_0),real(v_0)].';                      %�ֱ����ɵ���ʸ����ʵ�����鲿���ɵľ���va_0��vb_0

            a1=([0,va_0.';0,-vb_0.']).';                        %����a1����

            a2=zeros((array_number*2+1),3*silelobe_number);     %����a2����

            for i=1:silelobe_number      % i��1��116
                    v_i=exp(-1j*2*pi*fr*(0:(array_number-1))*array_space*sin(silelobe(i)*pi/180)/C);          %���㾭����ɢ��֮����԰�����
                                                                                                              %���Ƕ��ϵĵ���ʸ��v_i��
                    va_i=[real(v_i),imag(v_i)].';
                    vb_i=[-imag(v_i),real(v_i)].';              %�ֱ����ɵ���ʸ����ʵ�����鲿���ɵľ���va_i��vb_i��

                    a2(1:(array_number*2+1),1+3*(i-1):3*i)=([0,matrix_zero.';0,-va_i.';0,-vb_i.']).';         %����a2����
            end

            %*************����׶�滮�������⣨����matlab�Ĺ����䣩*************%
            C_col=([c1,c2,c3]).';

            A_col=([a1,a2,a3]).';

            [x,y,info]=sedumi(A_col,b,C_col,K);    %����sedumi����������׶�滮�������ز���

            w=y(2:array_number+1)+1j*y((array_number+2):(array_number*2+1));   %w�Ǽ��������Ȩֵ��������ȡֵ��y�����2��Nsl+1��Ԫ��
            
            
            
                     %% ****************����������Ƴ�����Ȩֵ��Ȩ�Ĳ����γ�***************       
            %*************�����γ�*************%
            for i=1:length(ang)
                    beam(i)=conj(exp(-1j*2*pi*fr*(0:(array_number-1))*array_space*sin(ang(i)*pi/180)/C))*w;           %beamΪ��Ȩ������Ӧ����
            end
%             beam1=max(abs(beam));
            beam_abs=abs(beam)/max(abs(beam));              %�Բ�����Ӧ���й�һ��
            beam_abs_1=beam_abs;
            %*************ȷ������������*************%
            MainLobe_down=find(diff(sign(diff(beam_abs)))>0)+1;    %diff�����������Ƕ�һ���������΢�֣���֣����㣻
                                                            %sign�������������ж�һ��ֵ�������������������ֵΪ1����������ֵΪ-1,0�򷵻�0��
                                                            %find���������������ڷ�������������������Ԫ�ص�����λ��(λ�õ��ж����ھ����У���һ�п�ʼ��
                                                            %���϶��£�����Ϊ1��2��3...,Ȼ���ٴӵڶ��У�����������������)
                                                            %�˴��������������ǣ�Ѱ�Ҳ�����Ӧ�İ��ݵ㼴����ͼ�йȵ׵�λ�ö�Ӧ�ĽǶȣ�
            temp1=ang(MainLobe_down);
            temp2=temp1-beam_point*ones(1,length(temp1));    
            ang_temp=min(abs(temp2))*2;                      %�����ang_temp����Ӧ�ľ��ǲ�������Ŀ�ȼ���������������
            
            %*************�԰꼶�����ָ��ıƽ��̶�*************%
            temp3=find(ang<(beam_point-ang_temp/2));
            sl_FIB=20*log10(max(beam_abs(temp3)));  
%             sl_FIB-sidelobe_level
                   
            %% ****************��ͼ**************** %%
            
            figure(1);hold on;grid on;
            
%             plot((beam_point-ang_temp/2).*ones(1,81),-80:0,'r--');
%             plot((beam_point+ang_temp/2)*ones(1,81),-80:0,'r--');

            plot(ang,20*log10(beam_abs));

%             plot(beam_point*ones(1,81),-80:0,'m--');
%             plot(ang,sl_FIB*ones(size(ang)),'m--');
            
            axis([min(ang),max(ang),-100,0]);

            xlabel('��λ/��');
            ylabel('����/dB');
            title('��������');
         
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
%         xlabel('�԰꼶(dB)');
% %         ylabel('�����ȣ��ȣ�');
% %         title('���������԰꼶�Ĺ�ϵ');
%         ylabel('��������(dB)');
%         title('��������������԰꼶�Ĺ�ϵ');
% %         axis([min(wide_c),max(wide_c),-10,10]);