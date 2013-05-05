function [book srr output] = MPdamped(y,winlen,hop,fftlen,maxiter)

nstride = ceil((fftlen-winlen)/hop);

booklen = fftlen*(nstride +1) ;
 
 book = zeros(booklen,1);
 
 res = y;
 its = 1;
 pow = 10*log10(res'*res);
 
 while(its<=maxiter)
  
 
   X = damptrans(res,winlen,hop,fftlen);
   [~,ind] = max(abs(X));
   alpha = X(ind);
   book(ind) = alpha;
   
   f = damp_synthesize(book,winlen,hop,fftlen);
   res = y-f;
   srr = pow - 10*log10(res'*res);
   fprintf('atom selected = %d , SRR = %d \n',ind,srr);         

   
   its = its+1;
   
 end
   output = real(damp_synthesize(book,winlen,hop,fftlen));
end