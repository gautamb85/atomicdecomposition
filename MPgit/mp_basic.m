%basic matching pursuit for damped sine dictionary
%seems to work - needs more testing
%need to implement conjugate subspace pursuit
% need to make a transient block of short (32) length
%damped sines

clear all;
clc;
% need to test with a dictionary atom as input
y = cos(2*pi*440*[0:2047]');
%y = y/norm(y);
res = y;
its =0;
book = zeros(2048*2047,1);

tic;
while (its <100)
    
    [alpha i] = fftcorr(res,2048);
    book(i) = alpha;
    atom = atomgen(i);
    
    res = res - real(alpha.*atom);
    its = its+1;
end
    mp_toc = toc;
    