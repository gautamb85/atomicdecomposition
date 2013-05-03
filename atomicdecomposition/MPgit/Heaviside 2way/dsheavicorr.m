function [alpha index] = dsheavicorr(x,fftlen)

Tf = 0.15;
Tb = 0.15;
winlen = 2048;

Corr = zeros(winlen,1);
ind = zeros(winlen,1);

window = zeros(winlen,1);

%for ptr = 0:winlen-1
    
    % forward & backward length
    bl = winlen/2;
    fl = bl;

    af = log(Tf)/fl;
    ab = log(Tb)/bl;

    % damped exponential window - forward & backward
    fwin = exp(af.*[0:fl-1]');
    fwin = fwin/norm(fwin);
    fpad = zeros(1024,1);
    fwin = [fpad;fwin];


    
    bwin = exp(-1*ab*[0:bl-1]');
    bwin = bwin/norm(bwin);
    bpad = zeros(bl,1);
    bwin = [bwin;bpad];
    % combined window
    window = fwin + bwin;
    window = window/norm(window);

for ptr = 0:winlen-1,
    
    chunk = fft(window .* x,fftlen);
     
     [~,i] = max(abs(chunk));
     Corr(ptr+1) = chunk(i);
     ind(ptr+1) = i + ptr*winlen;
end
                
[~,i1] = max(abs(Corr));
alpha = Corr(i1);
index = ind(i1);
end
