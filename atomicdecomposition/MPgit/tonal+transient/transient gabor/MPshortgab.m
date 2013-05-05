function [book srr output] = MPshortgab(y,winlen,hop,fftlen,maxiter)

nstride = ceil((fftlen-winlen)/hop);

booklen = fftlen*(nstride +1) ;
 
 book = zeros(booklen,1);
 
 res = y;
 pow = 10*log10(res'*res);

 its = 1;
 
 while(its<=maxiter)
  
 
   X = transgabcorr(res,winlen,hop,fftlen);
   [~,ind] = max(abs(X));
   
   alpha = X(ind);
   book(ind)=alpha;
   
   f = gabtrans_synthesize(book,winlen,hop,fftlen);
   res = y-f;
   srr = pow - 10*log10(res'*res);

   fprintf('atom selected = %d , SRR = %d \n',ind,srr);         

   its = its+1;
   
 end
  output = real(gabtrans_synthesize(book,winlen,hop,fftlen));
end