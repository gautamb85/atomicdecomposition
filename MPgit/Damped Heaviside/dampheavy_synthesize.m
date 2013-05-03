function op = dampheavy_synthesize(book)        

 T = 0.15;

 winlen = 2048;    
 fftlen = 2048;
            
  book = book(:); % make row vector
  blocklen = floor(length(book)/fftlen)-1;
  op = 0; % init. output buffer
            
            for ptr =0:blocklen
                
                l = winlen-ptr;
                a0 = log(T)/l;
                prepad = zeros(ptr,1);
                %postpad = zeros(fftlen - length(prepad)-l,1);
                
                window = [prepad;exp(a0*[0:l-1]')];
                window = window/norm(window);
                
                fft_ind = 1 + ptr*fftlen;
                fft_sup = fft_ind:fft_ind+fftlen-1;
                block = book(fft_sup);
                
                chunk = ifft(block, fftlen) * fftlen;
           
                op = op + window .* chunk; % mult synth win
                
            end
                                      
end
            
      