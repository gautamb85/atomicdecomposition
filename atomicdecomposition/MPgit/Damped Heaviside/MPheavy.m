function [book subset  corr srr output] = MPheavy(a,N,maxiter)

winlen = 2048;
booklen = N*(winlen-1);
 
 book = zeros(booklen,1);
 
 res = a;
 pow = 10*log10(res'*res);
 corr =[];
 subset =[];
 its = 1;
 
 while(its<=maxiter)
  
 
   [alpha ind] = dampheaviside(res,N);
   
   book(ind) = alpha;
   corr(its) = alpha;
   subset(its) = ind;
   
   f = dampheavy_synthesize(book);
   res = a-f;
   
   srr = pow - 10*log10(res'*res);
   fprintf('atom selected = %d , SRR = %d \n',ind,srr);         

   its = its+1;
   
 end
 output = real(dampheavy_synthesize(book));
end