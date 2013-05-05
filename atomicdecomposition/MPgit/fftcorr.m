%lets try to calculate correlations quickly
% function returns the coefficient value - alpha
% and the index of the atom selected by MP

function [alpha index] = fftcorr(x,N)

winlen = N;
Corr = zeros(2048,1);
ind = zeros(2048,1);

for ptr = 0:winlen-1,
                T = 0.01;
                l = winlen - ptr;
                a0 = log(T)/l;
                normal = sqrt((1-a0^2)/(1-(a0^2)^l));
                
                pad = zeros(ptr,1);
                
                win = [pad;exp(a0*[0:l-1]')];
                win = win/norm(win);
                
                chunk = fft(win .* x,N);

                [~,i] = max(abs(chunk));
                Corr(ptr+1) = chunk(i);
                ind(ptr+1) = i + ptr*N;
                
end
                [~,i1] = max(abs(Corr));
                alpha = Corr(i1);
                index = ind(i1);
end





