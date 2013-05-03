function [book output] = MPfastheavy(a,fftlen,maxiter)
%add variable time shifts
%stopping condition based on number of atoms
%modify to genetate either spectral or temporal reconstruction
winlen = 2048;
%window = window/norm(window);
%booklen = N*(N - l +1) ;
%booklen = 2048*(2048-63);
booklen = fftlen*(winlen - 1);
%initialize book
book = zeros(booklen,1);
%initialize correlation vector
%X = zeros(booklen,1);
%initialize residual
res = a;

its = 1;
 
X = heavicorr(res,fftlen);
 
while(its<=maxiter)
  
    [~,ind] = max(abs(X));
    alpha = X(ind);
    
    fprintf('atom selected = %d \n',ind);         
    
    book(ind) = alpha;
    
    [atom asup] = dampheavy_atomgen(ind);
    
    dx = crosscorr(real(alpha*atom),asup(1),asup(end));
    
    X = X - dx;
    its = its+1;
   
end

output = real(gabtrans_synthesize(book));
end