%function to generate the final book from 
%the combination of the tonal + transient blocks

% input signal length = 2048
% tonal block = 2048x8192
% transient block = 2048x4065280
 x = cos(2*pi*440*[0:2047]'/11000);
 N = 2048;
 [alpha1 ind1] = gabcorr(x,N);
 [alpha2 ind2] = damptrans(x,N);
 
 if (abs(alpha1)>abs(alpha2))
     index = ind1;
     
 else
     index = ind2;
 end




