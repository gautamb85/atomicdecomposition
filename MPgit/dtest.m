function C = dtest(res,window,scale,fftlen)
      
    hop = scale/2;
    
    l = length(res);
    %no. of hops
    nhop = ceil((l-scale)/hop);
 
    
    for ptr = 0:nhop,
        
        
        samp_ind = 1 + ptr*hop;
        corr_ind = 1+ ptr*fftlen;
        
        %for fftlen = scale
        %fcorr = fft(window.*ip(samp_ind ; samp_ind + scale),scale);
        
        fcorr = fft(window.*res(samp_ind : samp_ind + scale-1),fftlen);
        
        C(corr_ind:corr_ind+fftlen-1,1) = fcorr;
        
    end
  
 end    