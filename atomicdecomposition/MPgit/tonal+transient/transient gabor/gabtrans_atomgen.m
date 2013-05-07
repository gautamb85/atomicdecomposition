% GIVEN AN INDEX NUMBER
%Generate a Gabor atom of any size

function  [atom atomsup]  = gabtrans_atomgen(ind,window,hop,fftlen)


winlen = length(window);
%l = 64;
%window = hanning(winlen);
%window = window/norm(window);

bin = mod(ind-1,2048);                        
st = floor(ind/fftlen)*hop;

bl = floor(ind/2048);

%prepad = zeros(bl-1,1);
%postpad = zeros(winlen - length(prepad)-l,1);
     
                
sup = (0:winlen-1)';
atomsup = st+sup;               
            
atom = window .* exp(1j*2*pi/fftlen * bin * sup);
%atom = [prepad;atom;postpad];
end


