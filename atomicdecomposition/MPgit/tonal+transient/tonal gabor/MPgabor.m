function [book srr] = MPgabor(y,N,maxiter)
 
booklen = N;
 
 book = zeros(booklen,1);
 
 res = y;
 pow = 10*log10(res'*res);
 
 its = 1;
 
 while(its<=maxiter)
  
 
   [alpha ind] = gabcorr(res,N);
   
   book(ind) = alpha;
   f = synthesize(book,'t');
   
   res = y-f;
   
   srr = pow - 10*log10(res'*res);
   fprintf('atom selected = %d , SRR = %d \n',ind,srr);
   
   its = its+1;
   
 end
 
end