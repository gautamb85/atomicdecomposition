function y = gabsynthesize(X,window,stride,fftlen,M)        
        % synthesis matrix
        % assumptions:  1) X is aligned with fft block boundary
        %               2) the length of X is a multiple of the fftlen
            winlen = length(window);
            X = X(:); % make row vector
            blocklen = floor(length(X)/fftlen)-1;
            
            y = zeros(blocklen*stride+winlen,1); % init. output buffer
            
            for ptr = 0:blocklen,
                
                samp_ind = 1 + ptr*stride;
                samp_sup = samp_ind : samp_ind+winlen-1;
                fft_ind = 1 + ptr*fftlen;
                fft_sup = fft_ind:fft_ind+fftlen-1;
                
                block = X(fft_sup); % conv w/ synth win? --> do in time-domain instead
                chunk = ifft(block, fftlen) * fftlen;
                y(samp_sup) = y(samp_sup) + window .* chunk(1:winlen); % mult synth win
                
            end
            
            y = y(1:M);
end