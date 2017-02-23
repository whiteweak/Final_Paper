        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                 ���򲻱�㶨����խ����ά�����γ�             %
        %        ����ʵ�ֻ��ھ���ƽ����Ķ�ά���򲻱�㶨�������γ�   %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        clc;
        clear all;
        % close all
        %% ****************** �������� ****************** %%
        array_numx=15;  %ƽ����X���ϵ���Ԫ����
        array_numy=15;  %ƽ����Y���ϵ���Ԫ����

        f0=20000;       %�ο�Ƶ��

        c=340;          %�����ٶ�
        d=c/f0/2;       %��Ԫ���
        lamda=c/f0;     %��С����

        theta=(0:0.2:90);          %�����Ƿ�Χ
        phi=(-180:0.2:180);        %��λ�Ƿ�Χ
        %% ****************** ����SOCP�Ż���Ȩ����Ʋο�������Ӧ ****************** %%
        
        num=15;                  %��Ʋο�������Ӧ��Ӧ����Ԫ��
        vec=-(num-1)/2:(num-1)/2;
        
        wh=dlmread('coff_w_socp.txt');   %��������SOCP������Ƴ�����Ȩֵ
        wh_x=wh.';                       %x�����ϵļ�Ȩֵ��
        wh_y=wh.';                       %y�����ϵļ�Ȩֵ�� 
        
        % ****************** �ο�������Ӧ ******************%
        theta_p=46;  %�����ǲο�����ָ��
        phi_p=0;     %��λ�ǲο�����ָ��

        v_x=conj(wh_x.*exp(-1j*(vec)*pi*sin(theta_p*pi/180)*cos(phi_p*pi/180)));  %�ο���������x����ĵ���ʸ����
        v_y=conj(wh_y.*exp(-1j*(vec)*pi*sin(theta_p*pi/180)*sin(phi_p*pi/180)));  %�ο���������y����ĵ���ʸ����

        h1=kron(v_y,v_x);     %������ָ��ļ�Ȩʸ����Ϊ���������ϵ���ʸ���Ŀ����ڿ˻�  
     
        beam_temp=zeros(length(theta),length(phi));   %��ʼ���ο�������Ӧ����

        for k=1:length(theta)         %��ѭ���������ǶԸ����Ǻͷ�λ�ǹ۲�Ƕȷ�Χ�ڵ����нǶȽ��в����γɣ�������ɨ��
                for i=1:length(phi)
                        v_x=(exp(-1j*(vec)*pi*sin(theta(k)*pi/180)*cos(phi(i)*pi/180))).';
                        v_y=(exp(-1j*(vec)*pi*sin(theta(k)*pi/180)*sin(phi(i)*pi/180))).';
                        v=kron(v_y,v_x);

                        beam_temp(k,i)=(h1*v);     %���㲨����Ӧ
                end
        end
        beam_temp=abs(beam_temp)/max(max(abs(beam_temp)));      %�Բ�����Ӧ����ȡ����ֵ�͹�һ������

        %% ****************** ��������ϵ�µĲ���ͼ ****************** %%
        
         figure(1);
         Patt3d(beam_temp.',3);     %�����������µĲ���ͼ
         hold on;
         grid on;   
        
         
        %% ****************** ��ֱ������ϵ�µĲ���ͼ ****************** %%
        
        figure(2);    
        [x,y]=meshgrid(phi,theta);       %�Բ������зָ��γ����񣬱��ڻ�ͼ         
        mesh(x,y,20*log10(beam_temp));   %����ֱ������ϵ�µĲ���ͼ
        
        xlabel('��λ��/��');
        ylabel('������/��');
        zlabel('����/dB');
        h = gca;
        set(h,'FontSize',10.5,'FontName','����');
        set(h,'FontName','Times New Roman');  
        title('�ο�����ͼ');
         
        %% ****************** ��λ�Ƿ����ϵĲ���ͼ��������Ϊ����ָ��� ****************** %%
        
        figure(3);
        subplot(2,1,1);
        plot(phi,20*log10(beam_temp(find(theta==theta_p),:)),'r');   %������λ�Ƿ����ϣ��ο������Ƿ���Ĳ���ͼ
        hold on;
        grid on;
        
        beam_abs=beam_temp(find(theta==theta_p),:);    %Ѱ�ҵ�������ȷ��ʱ�ĸ�����λ���ϵĲ�����Ӧֵ
        
        % ****************** ȷ������������ ****************** %
        
        IndMin=find(diff(sign(diff(beam_abs)))>0)+1;
        temp=phi(IndMin);
        temp=temp-phi_p*ones(1,length(temp));    
        phi_wide1=min(abs(temp))*2;                  %ȷ����λ�Ƿ���������������

        % ****************** ȷ���԰꼶 ****************** %
        
        temp=find(phi<(phi_p-phi_wide1/2));
        sidelobe1=20*log10(max(beam_abs(temp)));     %ȷ����λ�Ƿ���������԰�ˮƽ���԰꼶
        % sidelobe1=sidelobe;
        % sl_phi=-25;

        % ****************** ȷ���԰꼶���� ****************** %
        
        temp=find((20*log10(beam_abs)-sidelobe1*ones(size(beam_abs)))>0.1); 
        SL_phiwide1=phi(temp(end));
        SL_phiwide2=phi(temp(1));
        SL_phiwide3=SL_phiwide1-SL_phiwide2;   %���㷽λ�Ƿ������԰꼶����      

        xlabel('��λ��(��)');
        ylabel('����(dB)');
        % title('�ο�����ͼ����λ�ǣ�');
        % legend('(30,0)','(50,0)')

        %% ****************** �����Ƿ����ϵĲ���ͼ����λ��Ϊ����ָ��� ****************** %
        
        subplot(2,1,2);
        plot(theta,20*log10(beam_temp(:,find(phi==phi_p))),'b');  %���������Ƿ����ϣ��ο���λ�Ƿ���Ĳ���ͼ
        hold on;
        grid on;

        beam_abs=beam_temp(:,find(phi==phi_p));   %Ѱ�ҵ���λ��ȷ��ʱ�ĸ����������ϵĲ�����Ӧֵ
        
        % ****************** ȷ������������ ****************** %
        
        IndMin=find(diff(sign(diff(beam_abs)))>0)+1;
        temp=theta(IndMin);
        temp=temp-theta_p*ones(1,length(temp));    
        theta_wide1=min(abs(temp))*2;                  %ȷ�������Ƿ���������������

         % ****************** ȷ���԰꼶 ****************** %
        temp=find(theta>(theta_p+theta_wide1/2));
        sidelobe2=20*log10(max(beam_abs(temp)));       %ȷ�������Ƿ���������԰�ˮƽ���԰꼶
        % sidelobe2=sidelobe;
        % sl_theta=-25;

         %------------ȷ���԰꼶����----------------------------%
        temp=find((20*log10(beam_abs)-sidelobe2*ones(size(beam_abs)))>0.2); 
        SL_thetawide1=theta(temp(end));
        SL_thetawide2=theta(temp(1));
        SL_thetawide3=SL_thetawide1-SL_thetawide2      %���㸩���Ƿ������԰꼶����

        xlabel('������(��)');
        ylabel('����(dB)');
        % title('�ο�����ͼ�������ǣ�');

        sidelobe=max(sidelobe1,sidelobe2);

        %% ****************** ����SOCP��Ʋ�ͬ����ָ���µļ�Ȩʸ�� ****************** %%
        
        phi_range=0;                        %��λ��ɨ�跶Χ��������Ϊ�̶�ֵ45�㣬���ټ��������ڻ�ͼ
        theta_range=20;                 %������ɨ�跶Χ
        w_coff=zeros(num*num,length(phi_range)*length(theta_range));      %��ʼ��Ȩֵ����
        
        theta_beam_mainlobewidth=zeros(1,length(theta_range));
        beam_theta_plus60=zeros(length(theta_range),length(theta));
                 %�˳��������ڸ����ǹ̶����ʴ洢���Ǹ����ǹ̶�Ϊ46��ʱ����λ�Ƿ����ϲ�ͬ�ǶȵĲ�����Ӧֵ���ڻ�ͼ�ͼ���������
                 
        for beamnum_phi=1:length(phi_range)
                for beamnum_theta=1:length(theta_range)
                        phi_point=phi_range(beamnum_phi);     %��λ�ǲ���ָ��
                        theta_point=theta_range(beamnum_theta);                   %�����ǲ���ָ��

                        temp1=SL_phiwide3;        %��λ�Ƿ����ϵ��԰꼶����
                        temp2=SL_thetawide3;      %�����Ƿ����ϵ��԰꼶����

                        mainlobe_theta=(theta_point-temp2/2):2:(theta_point+temp2/2);  %������ɢ��
                        mainlobe_phi=(phi_point-temp1/2):5:(phi_point+temp1/2);  %������ɢ��

                        silelobe_theta=[0:2:(theta_point-temp2/2),(theta_point+temp2/2):2:90];  %�԰���ɢ��
                        silelobe_phi=[-180:10:(phi_point-temp1/2),(phi_point+temp1/2):10:180];  %�԰���ɢ��

                        slphi_num=length(silelobe_phi);  
                        sltheta_num=length(silelobe_theta);
                        silelobe_num=slphi_num*sltheta_num;   %�԰���ɢ��

                        mlphi_num=length(mainlobe_phi);  
                        mltheta_num=length(mainlobe_theta);
                        mainlobe_num=mlphi_num*mltheta_num;  %�԰���ɢ��

                % ****************** ��ʼ�� ****************** %
                        matrix_zero=zeros(num*num*2,1);
                 
                        beam_d=zeros(mltheta_num,mlphi_num);               %�ο�������Ӧ
                
                        phi_temp=mainlobe_phi-(phi_point-phi_p)*ones(size(mainlobe_phi));         %�ο�������Ӧ�ķ�λ�Ƿ�Χ
                        theta_temp=mainlobe_theta-(theta_point-theta_p)*ones(size(mainlobe_theta));   %�ο�������Ӧ�ĸ����Ƿ�Χ

                % ****************** ����������Ӧ ****************** %
                        for k=1:mltheta_num
                                for i=1:mlphi_num
                                        v_x=(exp(-1j*(vec)*pi*sin(theta_temp(k)*pi/180)*cos(phi_temp(i)*pi/180))).';
                                        v_y=(exp(-1j*(vec)*pi*sin(theta_temp(k)*pi/180)*sin(phi_temp(i)*pi/180))).';
                                        v=kron(v_y,v_x);  

                                        beam_d(k,i)=(h1*v);                          
                                end
                        end
                
                        beam_d=beam_d/max(max(abs(beam_d)));  %��һ��


              %% ****************** ����׶�滮��Ʒ��򲻱�㶨����խ����ά������Ӧ ****************** %%
                
                % ****************** ���Ʋ���ָ�����ӦΪ1 ****************** %
                        cc1=([1;0]).';

                        v_x=(exp(-1j*(vec)*pi*sin(theta_point*pi/180)*cos(phi_point*pi/180))).';
                        v_y=(exp(-1j*(vec)*pi*sin(theta_point*pi/180)*sin(phi_point*pi/180))).';
                        v_0=kron(v_y,v_x);  

                        va_0=[real(v_0);imag(v_0)];
                        vb_0=[-imag(v_0);real(v_0)];
                        aa1=([0,zeros(1,mainlobe_num),va_0.';0,zeros(1,mainlobe_num),-vb_0.']).';

                % ****************** �������������С��ĳ���� ****************** %
                        cc2=1;  %0.4
                        aa2=([0,ones(1,mainlobe_num),zeros(1,2*num*num)]).';

                % ****************** ���������� ****************** %
                        cc3=zeros(1,4*mainlobe_num);
                        aa3=zeros((num*num*2+mainlobe_num+1),4*mainlobe_num);

                % ****************** �����԰�ˮƽ ****************** %
                        cc4=zeros(1,3*silelobe_num);
                        aa4=zeros((num*num*2+mainlobe_num+1),3* silelobe_num);

                % ****************** ���Ʋ����Ƚ��� ****************** %
                        cc5=([0.1;matrix_zero]).';
                        aa5=([zeros(1,mainlobe_num+1),matrix_zero.';zeros(2*num*num,mainlobe_num+1),-eye(num*num*2,num*num*2)]).';  
                
                % ****************** ��������B ****************** %
                        bb=([-1,zeros(1,mainlobe_num),matrix_zero.']).';

                % ****************** ����׶Լ����ά�� ****************** %
                        Q.f=2;
                        Q.l=1;
                        Q.q=[4*ones(1,mainlobe_num),3*ones(1,silelobe_num),num*num*2+1];

                % ****************** ����C��ת�ú;���A ****************** %
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

                % ****************** �������� ****************** %
                        CC_col=([cc1,cc2,cc3,cc4,cc5]).';

                        AA_col=([aa1,aa2,aa3,aa4,aa5]).';

                        [x,y,info]=sedumi(AA_col,bb,CC_col,Q);
                        y(1)

                        w=y(mainlobe_num+2:mainlobe_num+num*num+1)+1j*y((mainlobe_num+num*num+2):(mainlobe_num+num*num*2+1));
                        w=w/max(abs(w));

                % ****************** ���򲻱�㶨�������γ� ****************** %               
                        for k=1:length(theta)
                                for i=1:length(phi)
                                        v_x=(exp(-1j*(vec)*pi*sin(theta(k)*pi/180)*cos(phi(i)*pi/180))).';
                                        v_y=(exp(-1j*(vec)*pi*sin(theta(k)*pi/180)*sin(phi(i)*pi/180))).';
                                        v=kron(v_y,v_x);


                                        beam_temp(k,i)=v'*w;
                                end
                        end

                        beam_temp=abs(beam_temp)/max(max(abs(beam_temp)));      %��һ��

                        w_coff(:,(beamnum_phi-1)*length(theta_range)+beamnum_theta)=w;
                        
                        beam_theta_plus60(beamnum_theta,:)=beam_temp(:,find(phi==phi_point));
                        
                        %% ****************** ������λ�ǹ̶�ʱ�����������温���Ǳ仯���� ****************** %%
%                           
%                 % ****************** ȷ������������ ****************** %
%                 
%                         IndMin=find(diff(sign(diff(beam_theta_plus50(beamnum_theta,:))))>0)+1;
%                         temp=theta(IndMin);
%                         temp=temp-theta_range(beamnum_theta)*ones(1,length(temp));    
%                         theta_wide1=min(abs(temp))*2;                  %ȷ����λ�Ƿ���������������
% 
%                 % ****************** ȷ���԰꼶 ****************** %
%         
%                         temp=find(theta<(theta_range(beamnum_theta)-theta_wide1/2));
%                         sidelobe1=20*log10(max(beam_theta_plus50(temp)));     %ȷ����λ�Ƿ���������԰�ˮƽ���԰꼶
%                         % sidelobe1=sidelobe;
%                         % sl_phi=-25;
% 
%                 % ****************** ȷ���԰꼶���� ****************** %
%         
%                         temp=find((20*log10(beam_theta_plus50(beamnum_theta,:))-sidelobe1*ones(1,length(theta)))>0.1); 
%                         SL_thetawide1=theta(temp(end));
%                         SL_thetawide2=theta(temp(1));
%                         SL_thetawide3=SL_thetawide1-SL_thetawide2;   
%                                     %���㷽λ�Ƿ������԰꼶����,SL_phiwide3��Ϊ�ڹ̶�������ʱ��������λ���ϵĲ��������԰꼶����ֵ  
%                                     
%                           %% ****************** �洢��������Ĳ���������ֵ���԰꼶������߲������������� ****************** %%            
% %                         phi_beam_mainlobewidth(1,beamnum_phi)=SL_phiwide3;  
%                                     %��ʱphi_beam_sidelobewidth����洢���Ƿ�λ�Ƿ����������԰꼶����ֵ�����ں�����ͼ
%                         theta_beam_mainlobewidth(1,beamnum_theta)=theta_wide1;  
%                                     %��ʱphi_beam_sidelobewidth����洢���Ƿ�λ�Ƿ��������겨������������ֵ�����ں�����ͼ


                %% ****************** ��ͼ ****************** %
                
%                         figure(4);   %����������ϵ�µĲ���ͼ
%                         Patt3d(beam_temp.',3);

                        figure(5);    %����ֱ������ϵ�µĲ���ͼ
                        [x,y]=meshgrid(phi,theta);
                        mesh(x,y,20*log10(beam_temp));

                        xlabel('��λ��(��)');
                        ylabel('������(��)');
                        zlabel('����(dB��');
                
                        h = gca;
                        set(h,'FontSize',10,'FontName','����');
                        set(h,'FontName','Times New Roman'); 
                
                        title('����SOCP�ķ��򲻱䲨���γ�');
              
                        figure(6);
                % ****************** ��λ�Ƿ����ϵĲ���ͼ��������Ϊ�ο�����ָ������ ****************** %
%                         subplot(2,1,1);
                        plot(phi,20*log10(beam_temp(find(theta==theta_point),:)),'k');hold on;grid on;

                        h = gca;
                        set(h,'FontSize',10,'FontName','����');
                        set(h,'FontName','Times New Roman');  

                        xlabel('��λ��(��)');
                        ylabel('����(dB��');
                        axis([-180,180,-80,0]);
        
                % ****************** �����Ƿ����ϵĲ���ͼ����λ��Ϊ�ο�����ָ��λ�� ****************** %
                        figure(7);
                        plot(theta,20*log10(beam_temp(:,find(phi==phi_point))),'k');hold on;grid on;

                        h = gca;
                        set(h,'FontSize',10,'FontName','����');
                        set(h,'FontName','Times New Roman'); 
                
                        xlabel('������(��)');
                        ylabel('����(dB��');
%                         axis([0,90,-50,0]);
               

                end
        end

%                 dlmwrite('coeffcient_phi0_theta60_all.txt',w_coff,';');       %��Ȩֵ���Ϊ.txt�ļ�����������ĵ������޸�

%                 figure(8);
%                 theta=(30:2:66);
%                 plot(theta,theta_beam_mainlobewidth(1,:),'r*');grid on;
