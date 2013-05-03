function X = gabrep_kern(atom, st,ed,window,stride,fftlen,M)
        % fast computation of atom-to-dictionary inner products
        % (similar to analyze function except with limited range)
            
        
            winlen = length(window);
            nstrides = ceil((M-winlen)/stride);
            N = (nstrides+1)*fftlen;
                     
            st_block = max(0, floor(st/stride) - floor(winlen/stride-1));
            ed_block = min(floor(ed/stride), nstrides);
            blocklen = ed_block-st_block;
            
            pre_pad = zeros(st - st_block*stride, 1);
            post_pad = zeros(ed_block*stride + winlen - ed - 1, 1);                        
            y = [pre_pad; atom; post_pad];            
            
            X = zeros(N,1); % init. output buffer
            
            for ptr = 0:blocklen,
                
                samp_ind = (0 + ptr)*stride + 1;
                fft_ind = (st_block + ptr)*fftlen + 1;
                
                chunk = window .* y(samp_ind : samp_ind+winlen-1);
                X(fft_ind : fft_ind+fftlen-1) = fft(chunk, fftlen);

            end
        end