

%      clc;
%      clear;
    %% ********** »­Í¼ ********** %%
    
    theta=(0:0.25:90); 
    beam_temp=beam_theta_plus30to68;
    figure(1);
    plot(theta,20*log10(beam_temp),'k');hold on;grid on;
    h = gca;
    set(h,'FontSize',10,'FontName','ËÎÌå');
    set(h,'FontName','Times New Roman'); 
                
    xlabel('¸©Ñö½Ç(¶È)');
    ylabel('²¨Êø(dB£©');
%     axis([0,90,-40,0]);
    
    %% ********** ¼ÆËã¸÷²¨ÊøÖ÷°ê¿í¶È²¢»­Í¼ ********** %%    
    beam_temp_max=max(beam_temp');
    theta1=30:2:68;
    beam_3db_width=zeros(1,length(theta1));
    for i=1:20
            indmin=find(0.5<((beam_temp_max(1,i))-beam_temp(i,:))<0.6);
            indmin=max(indmin);
%             temp=beam_temp(i,indmin)
            beam_3db_width(1,i)=(abs(theta(indmin)-theta1(i)))*2;
%             beam_3db_width(1,i)=
    end
    
%     figure(2);
    
%     plot(theta,beam_3db_width,'k');
%     hold on;grid on;
    
    