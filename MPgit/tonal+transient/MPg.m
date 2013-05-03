%% Matching Pursuit for a union of frames

% Gabor block uses a fftlen = 8192 (make it user defined)
% transient block uses a window of length 64 (make it user defined);
% incorporate SRR calucation in this fuction and make a stopping condition
% based on SRR
% tonal block preceeds the transient block 
% support of tonal block : tonal_sup = 1:8192;
% support for transient block : trans_sup = 8192+1:8192 + (2048*1985) -1;
 


function [book] = MPcomb(y,N,maxiter)
 
booklen = 2048 + 2048*1985;
 
 book = zeros(booklen,1);
 
 res = y;
 its = 1;
 
 while(its<=maxiter)
  
 
   [alpha1 ind1] = gabcorr(res,N);
    % book(ind) = alpha1;
  [alpha2 ind2] = damptrans(res,N);
 
% combined : transient + tonal 
 if (abs(alpha1)>abs(alpha2))
     alpha = alpha1;
     atom = gab_atomgen(ind1);
     book(ind1) = alpha1;
 else
     alpha = alpha2;
     atom = trans_atomgen(ind2);
     book(8192 +ind2) = alpha2;
 end
 
        f = synthesize(book);
        res = y - f;
        its = its+1;
 end
 
end

