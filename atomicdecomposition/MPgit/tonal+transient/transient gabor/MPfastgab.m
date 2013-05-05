function [book output] = MPfastgab(a,fftlen,maxiter)
%add variable time shifts
%stopping condition based on number of atoms
%modify to genetate either spectral or temporal reconstruction
winlen = 64;
%window = window/norm(window);
%booklen = N*(N - l +1) ;
%booklen = 2048*(2048-63);
booklen = fftlen*(fftlen - winlen +1);
%initialize book
book = zeros(booklen,1);
%initialize correlation vector
%X = zeros(booklen,1);
%initialize residual
res = a;
pow = 10*log10(res'*res);
its = 1;
 
X = shortgabcorr(res,fftlen);
 
while(its<=maxiter)
  
    [~,ind] = max(abs(X));
    alpha = X(ind);
        
    book(ind) = alpha;
    
    [atom asup] = gabtrans_atomgen(ind);
    
    dx = crosscorr(real(alpha*atom),asup(1),asup(end));
    
    X = X - dx;
    
    yhat = gabtrans_synthesize(book);
    res = a - yhat;
    srr = pow - 10*log10(res'*res);
    fprintf('atom selected = %d , SRR = %d \n',ind,srr);         

    
    its = its+1;
   
end

output = real(gabtrans_synthesize(book));
end