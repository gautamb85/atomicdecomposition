function [atom st bin] = genscaleatom(ind,selscale,fftlen)
            
shop = selscale/2;

%time block in which the atom lies

%amount of zero padding needed
st = floor(ind/fftlen)*shop;

bin = mod(ind-1,fftlen);

window = hanning(selscale);
window = window/norm(window);

sup = [0:selscale-1]';
atom = window.*exp(1j*2*pi/fftlen*bin*sup);
atomsup = st+sup;
end
        