% function to calculate correlations for short length(64)
% gabor atoms, time shift = 1 sample

function X = transgabcorr(x,winlen,hop,fftlen) 

 % Analysis window 
 win = hanning(winlen);
 %window = window/norm(window);

 
 %hop size
 
%  %preprocessing signal
%  %padding with a hop size of zeros
%  pad = zeros(winlen/2,1);
%  %postpad so that signal is a multiple of hop length
%  l = winlen/2 + length(x);
%  post = ceil(l/hop)*hop - l;
%  postpad = zeros(post,1);
%  
%  ip = [pad;x;postpad];
%  
 %build in a hop size

 nhop = ceil((length(x)-winlen)/hop) ;
 
 %Corr = zeros(fftlen,1);
 %ind = zeros(fftlen,1);
 
 %all correlations
 X = [];
      
       
   for ptr = 0:nhop,  
     
     prepad = zeros(hop*ptr,1);    
     postpad = zeros(fftlen - length(prepad)-winlen,1);
     window = [prepad;win;postpad];
     window = window/norm(window);

        %samp_ind = 1+ptr*hop;
        fft_ind = 1 + ptr*fftlen;
        %tmp = ip(samp_ind:samp_ind +winlen-1);
        chunk = fft((window .*x),fftlen);
     
        X(fft_ind:fft_ind+fftlen-1) = chunk;   
%         [~,i] = max(abs(chunk));
%         Corr(ptr+1) = chunk(i);
%         ind(ptr+1) = i + ptr*fftlen;   
        
   end
 
%   [~,i1] = max(abs(Corr));
%   alpha = Corr(i1);
%   index = ind(i1);

end
 