function op = synthesize(book,choice)        
       %currently trying to synthesize o/p using gabor atoms only
            window = hanning(2048);
            window = window/norm(window);
            fftlen = 2048;
            
            book = book(:); % make row vector
            blocklen = floor(length(book)/fftlen)-1;
            op = 0; % init. output buffer
            
            for ptr =0:blocklen
                
                fft_ind = 1 + ptr*fftlen;
                fft_sup = fft_ind:fft_ind+fftlen-1;
                block = book(fft_sup);
                
                chunk = ifft(block, fftlen) * fftlen;
                
                if(choice == 't')
                    op = op + window .* chunk; % time domain synthesis
                elseif(choice == 'f')
                    op = op + fft(window .* chunk,fftlen); % time domain synthesis
                end
                
            end
            
end
            
      