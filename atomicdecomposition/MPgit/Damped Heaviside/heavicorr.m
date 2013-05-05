function X = heavicorr(x,N) 

 %window parameters
 T = 0.15;
 fftlen = N;
 winlen = N;
 nstrides = winlen-1;
 X =zeros(N*nstrides,1);
 
 
 for ptr = 0:winlen-1,
     
     l = winlen - ptr; 
     a0 = log(T)/l;
    
     win = exp(a0*[0:l-1]');
 
     prepad = zeros(ptr,1);
     %one sided exponentially damped window
     dampwin = [prepad;win];
     dampwin = dampwin/norm(dampwin);
     
     fftind = 1 +ptr*fftlen;

     chunk = dampwin .* x;
     X(fftind:fftind+fftlen-1) = fft(chunk,N); 
     
 end
                

