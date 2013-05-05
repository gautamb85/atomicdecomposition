function X = gabanalyze(y,window,stride,fftlen)         
        % analysis matrix
        % assumptions:  1) y is aligned with a fft block boundary
        %
        % (caller may need to pad y with extra blocks to get desired
        % boundary conditions right)
            M = length(y);
            winlen = length(window);
            nstrides = ceil((M-winlen)/stride);
            pad = nstrides*stride + winlen - M;
            N = (nstrides+1)*fftlen;
            % make row vector, pad if necessary
            y = y(:);                         
            y = [y; zeros(pad,1)];
            
            X = zeros(N,1); % init. output buffer
            
            % fft analysis
            for ptr = 0:nstrides,
                
                samp_ind = 1 + ptr*stride;
                fft_ind = 1 + ptr*fftlen;
                
                chunk = window .* y(samp_ind : samp_ind+winlen-1);
                X(fft_ind : fft_ind+fftlen-1) = fft(chunk, fftlen);                
            end                 
        end         