% calculate atom-atom correlation
% based on the atom selected, the aim is to calculate the 
% correlations with all overlapping atoms
% input is an atom index

% function for atom-atom or window cross correlation

function X = crosscorr(atom,st,ed,window,hop,fftlen,M)
        % fast computation of atom-to-dictionary inner products
        % (similar to analyze function except with limited range)
            
        
            winlen = length(window);
            %fftlen = 2048;
            nstride = ceil((M-winlen)/hop);
            N = (nstride+1)*fftlen;    
            
            st_block = max(0, floor(st/hop) - floor(winlen/hop-1));
            ed_block = min(floor(ed/hop),nstride);
            
            blocklen = ed_block-st_block;
            
            pre_pad = zeros(st - st_block*hop, 1);
            post_pad = zeros(ed_block*hop + winlen - ed - 1, 1);
            
            y = [pre_pad; atom; post_pad];            
            
            X = zeros(N,1); % init. output buffer
            
            for ptr = 0:blocklen,
                
                samp_ind = (0+ptr)*hop + 1;
                fft_ind = (st_block + ptr)*fftlen + 1;
                
                chunk = window .* y(samp_ind : samp_ind+ winlen-1);
                X(fft_ind : fft_ind+fftlen-1) = fft(chunk,fftlen);

           end   
       
end