function [xc,yc,zc]=CrossSecT(thetac,phic,ratioc,sclc,rminc,beam_temp)
        [THETAC,PHIC]=meshgrid(thetac,phic);
%         rc=f(THETAC,PHIC);
%         rc=beam_f(thetac,phic);
        for i=1:length(thetac)
                rc(:,i)=beam_temp;
        end
        rc = Scale(rc,ratioc,sclc,rminc);
        % Mesh the cross-section plane
        [IC,JC]=size(rc);
        for i=1:IC
            rc(i,JC)=0;
        end
        [xc,yc,zc]=Sph2rec(THETAC,PHIC,rc);

     