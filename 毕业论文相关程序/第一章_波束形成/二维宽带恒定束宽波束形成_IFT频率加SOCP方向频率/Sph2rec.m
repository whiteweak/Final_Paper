function [xt,yt,zt]=Sph2rec(thetat,phit,rhot)
        xt=rhot.*sin(thetat).*cos(phit);
        yt=rhot.*sin(thetat).*sin(phit);
        zt=rhot.*cos(thetat);