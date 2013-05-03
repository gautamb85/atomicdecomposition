function op = gabtrans_synthesize(book,window,hop,fftlen,M)        

        
       %currently trying to synthesize o/p using gabor atoms only
            %l = 64;
            winlen = length(window);
            %window = window/norm(window);
            %fftlen = 2048;
            
            book = book(:); % make row vector
            blocklen = floor(length(book)/fftlen)-1;
            %op = zeros(blocklen*hop+winlen,1);  
            %nhop = ceil((M-winlen)/hop);
            op = zeros(blocklen*hop+winlen,1); % init. output buffer

            for ptr =0:blocklen,
                  
%                window = hanning(winlen);
%                window = window/norm(window);
                
                samp_ind = 1 + ptr*hop;
                samp_sup = samp_ind : samp_ind+winlen-1;
                
                fft_ind = 1 + ptr*fftlen;
                fft_sup = fft_ind:fft_ind+fftlen-1;
                
                block = book(fft_sup); % conv w/ synth win? --> do in time-domain instead
                chunk = ifft(block, fftlen) * fftlen;
                op(samp_sup) = op(samp_sup) + window.* chunk(1:winlen);
                    
            end
            
            %op = op(1:M);
            
end
            
      
                                      

            
      