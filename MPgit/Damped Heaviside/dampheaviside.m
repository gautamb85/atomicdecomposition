function [alpha index] = dampheaviside(x,N) 

 %window parameters
 T = 0.15;
 winlen = N;
 nstrides = winlen-1;
 
 Corr = zeros(winlen,1);
 ind = zeros(winlen,1);
 
 for ptr = 0:winlen-1,
     
     l = winlen - ptr; 
     a0 = log(T)/l;
    
     win = exp(a0*[0:l-1]');
 
     prepad = zeros(ptr,1);
     %one sided exponentially damped window
     dampwin = [prepad;win];
     dampwin = dampwin/norm(dampwin);
     
     chunk = fft(dampwin .* x,N);
     
     [~,i] = max(abs(chunk));
     Corr(ptr+1) = chunk(i);
     ind(ptr+1) = i + ptr*winlen;
     
 end
                
[~,i1] = max(abs(Corr));
  alpha = Corr(i1);
  index = ind(i1);
end
