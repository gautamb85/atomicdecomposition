function [book subset  corr srr output] = MP2cheavy(a,N,maxiter)

winlen = 2048;
booklen = N;
 
 book = zeros(booklen,1);
 
 res = a;
 pow = 10*log10(res'*res);
 corr =[];
 subset =[];
 its = 1;
 
 while(its<=maxiter)
  
 
   [alpha ind] = dsheavicorr(res,N);
   
   book(ind - 1) = alpha;
   corr(its) = alpha;
   subset(its) = ind;
   
   f = dsheavy_synth(book);
   res = a-real(f);
   
   srr = pow - 10*log10(res'*res);
   fprintf('atom selected = %d , SRR = %d \n',ind,srr);         

   its = its+1;
   
 end
 output = real(dsheavy_synth(book));
end