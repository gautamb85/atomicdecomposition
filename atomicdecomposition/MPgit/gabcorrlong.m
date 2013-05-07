%lets try to calculate correlations quickly
% function returns the coefficient value - alpha
% and the index of the atom selected by MP

function [alpha index] = gabcorrlong(x,N)

winlen = N;
fftlen = 8192;

win = hanning(N);
win = win/norm(win);
X=[];

for ptr = 0:winlen-1,
                
                fft_ind = 1+ ptr*fftlen;
                
                chunk = fft(win .* x,fftlen);

                X(fft_ind : fft_ind +fftlen-1) = chunk;
                
                
end
                [~,i1] = max(abs(X));
                alpha = X(i1);
                index = i1;
end





