function Bdyplot(thetab,phib,ratiob,sclb,rminb)
        rb= Scale(f(thetab,phib),ratiob,sclb,rminb);
        [xb,yb,zb]=Sph2rec(thetab,phib,rb);
        plot3([0,xb],[0,yb],[0,zb],'k');