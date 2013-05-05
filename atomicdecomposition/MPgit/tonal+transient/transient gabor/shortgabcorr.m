% function to calculate correlations for short length(64)
% gabor atoms, time shift = 1 sample
% M = input signal length without zero padding

function X = shortgabcorr(x,window,hop,fftlen,M) 

 %window parameters
 x = x(:);
 winlen = length(window);
 
 nstride = ceil((M-winlen)/hop);
 
 pad = nstride*hop + winlen - M;
 x = [x;zeros(pad,1)];
 %window = hanning(winlen);
 %window = window/norm(window);

 %pad if needed
 
 X = zeros(((nstride+1)*fftlen),1);
 
 for ptr = 0:nstride,
         
     samp_ind = 1 + ptr*hop;
     fftind = 1 +ptr*fftlen;
     
     chunk = window .* x(samp_ind:samp_ind+winlen-1);
    
     X(fftind:fftind+fftlen-1) = fft(chunk,fftlen);
     
     
 end
                

end
 