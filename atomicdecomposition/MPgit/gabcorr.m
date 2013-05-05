%lets try to calculate correlations quickly
% function returns the coefficient value - alpha
% and the index of the atom selected by MP

function [alpha index] = gabcorr(x,N)

winlen = N;
fftlen = N;


win = hanning(N);
win = win/norm(win);


                
chunk = fft(win .* x,fftlen);

[~,i] = max(abs(chunk));

alpha = chunk(i);
index = i;

end





