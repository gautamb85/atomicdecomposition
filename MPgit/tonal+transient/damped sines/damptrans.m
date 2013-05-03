function [X] = damptrans(x,winlen,hop,fftlen) 

 %window parameters
 T = 0.15;
 %winlen = 64;
 %winlen = N;
 a0 = log(T)/winlen;
 X =[];
 %Corr = zeros(fftlen,1);
 %ind = zeros(fftlen,1);
 
 nhop = ceil((length(x)-winlen)/hop) ;

 for ptr = 0:nhop,  

 %for ptr = 0:winlen-l,
     
     win = exp(a0*[0:winlen-1]');
 
     prepad = zeros(hop*ptr,1);
 
     postpad = zeros(fftlen - length(prepad)-winlen,1);
    
     %exponentially damped window
     dampwin = [prepad;win;postpad];
     dampwin = dampwin/norm(dampwin);
        
     chunk = fft(dampwin .* x,fftlen);
     fft_ind = 1 + ptr*fftlen;
     X(fft_ind:fft_ind+fftlen-1) = chunk;   

%      [~,i] = max(abs(chunk));
%      Corr(ptr+1) = chunk(i);
%      ind(ptr+1) = i + ptr*winlen;

     
 end
                
end 