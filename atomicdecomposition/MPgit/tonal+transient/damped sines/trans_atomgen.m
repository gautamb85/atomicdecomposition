%FUNCTION TO GENERATE A DAMPED SINE ATOM
% GIVEN AN INDEX NUMBER
function  atom  = trans_atomgen(ind)

T = 0.15;
winlen = 2048;
bin = mod(ind-1,2048);                        
bl = ceil(ind/2048);

l = 64;

a = log(T)/l;

prepad = zeros(bl-1,1);
postpad = zeros(winlen - length(prepad)-l,1);
                
dampwin = exp(a*[0:l-1]');
dampwin = dampwin/norm(dampwin);
                
sup = (0:l-1)';
               
            
atom = dampwin .* exp(1j*2*pi/2048 * bin * sup);
atom = [prepad;atom;postpad];
end


