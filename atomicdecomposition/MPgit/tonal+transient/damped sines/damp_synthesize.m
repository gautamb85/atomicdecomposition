function op = damp_synthesize(book,winlen,hop,fftlen)        

 T = 0.15;
 %winlen = 64;
 a0 = log(T)/winlen;
     
 %fftlen = 2048;
            
  book = book(:); % make row vector
  op = 0; % init. output buffer
  nhop = ceil((fftlen-winlen)/hop);
          
            for ptr =0:nhop
                
                prepad = zeros(hop*ptr,1);
                postpad = zeros(fftlen - length(prepad)-winlen,1);
                
                window = [prepad;exp(a0*[0:winlen-1]');postpad];
                window = window/norm(window);
                
                fft_ind = 1 + ptr*fftlen;
                fft_sup = fft_ind:fft_ind+fftlen-1;
                block = book(fft_sup);
                
                chunk = ifft(block, fftlen) * fftlen;
           
                op = op + window .* chunk; % mult synth win
                
            end
                                      
end
            
      