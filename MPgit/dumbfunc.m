% function to chose the dictionary block with the highest correlation
% and generates the appropriate atom

function  X = dumbfunc(x,N) 

 %window parameters
 T = 0.01;
 l = 64;
 winlen = N;
 a0 = log(T)/l;
 
 X =[];

 
 win = exp(a0*[0:l-1]');
 
 for ptr = 0:winlen-l,
 
     prepad = zeros(ptr,1);
 
     postpad = zeros(winlen - length(prepad)-l,1);
     
     fft_ind1 = 1+ptr*N;
     
     %exponentially damped window
     dampwin = [prepad;win;postpad];
     dampwin = dampwin/norm(dampwin);
        
        
     chunk = fft(dampwin .* x,N);
     
     X(fft_ind1: fft_ind1 + N-1) = chunk;
     
     [~,i] = max(abs(chunk));
    
 end
                
end
 