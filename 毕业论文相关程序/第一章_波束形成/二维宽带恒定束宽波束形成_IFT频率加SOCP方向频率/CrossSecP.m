   % Calulate Cross-section of the Cut - Phi
        function [xc,yc,zc]=CrossSecP(thetac,phic,ratioc,sclc,rminc,beam_temp)
        [THETAC,PHIC]=meshgrid(thetac,phic);
%         rc=f(THETAC,PHIC);
%         rc=beam_f(thetac,phic);
        for i=1;length(phic)
                rc(i,:)=beam_temp;
        end
        rc = Scale(rc,ratioc,sclc,rminc);
        % Mesh the cross-section plane
        [IC,JC]=size(rc);
        for j=1:JC
            rc(IC,j)=0;
        end
        [xc,yc,zc]=Sph2rec(THETAC,PHIC,rc);