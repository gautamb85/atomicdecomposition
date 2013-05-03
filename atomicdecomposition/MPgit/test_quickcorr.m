   T = 0.01;
   l = 400;
   a0 = log(T)/l;
   normal = sqrt((1-a0^2)/(1-(a0^2)^l));
   
   s = exp(a0*[0:l-1]').*cos(2*pi*440*[0:l-1]'/11000);
   pad = zeros((2048-400),1);
               
   
   sf = [pad;s];
   
   [C index X] = fftcorr(sf,2048);
   
   [m iu] = max(abs(X));
   
   C1 = X(iu);