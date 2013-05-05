% function to calculate correlations for short length(64)
% gabor atoms, time shift = 1 sample

function X = damptrans1(x,N) 

 %window parameters
 l = 64;
 winlen = N;
 win = hanning(l);
 fftlen = N;
 X = zeros(N*(N-l+1),1);
 
 for ptr = 0:winlen-l,
      
     prepad = zeros(ptr,1);
 
     postpad = zeros(winlen - length(prepad)-l,1);
    
     window = [prepad;win;postpad];
     window = window/norm(window);
     
     fftind = 1+ptr*fftlen;
     
     chunk = fft(window .* x,N);
     X(fftind:fftind + fftlen -1) = chunk;
    
     
 end
                

end
 