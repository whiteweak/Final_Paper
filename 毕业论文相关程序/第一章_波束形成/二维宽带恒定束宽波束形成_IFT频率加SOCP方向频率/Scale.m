function fs=Scale(rs,ratios,scls,rmins)
        if scls < 3
           fs=(rs/ratios).^scls;
        else
           rs=(rs/ratios);
           fs=20*log10(rs);
           idxs=find(fs>rmins);
           fs(idxs)=fs(idxs)-rmins;
           idxs=find(fs<=rmins);
           fs(idxs)=0;
        end