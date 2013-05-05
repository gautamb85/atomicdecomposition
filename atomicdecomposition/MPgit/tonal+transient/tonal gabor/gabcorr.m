%lets try to calculate correlations quickly
% function returns the coefficient value - alpha
% and the index of the atom selected by MP

function [alpha index] = gabcorr(x,N)

winlen = N;
fftlen = 8192;
Corr = zeros(2048,1);
ind = zeros(2048,1);

win = hanning(N);
win = win/norm(win);

for ptr = 0:winlen-1,
                
                
                
                
                
                chunk = fft(win .* x,fftlen);

                [~,i] = max(abs(chunk));
                Corr(ptr+1) = chunk(i);
                ind(ptr+1) = i + ptr*N;
                
end
                [~,i1] = max(abs(Corr));
                alpha = Corr(i1);
                index = ind(i1);
end





