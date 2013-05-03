%FUNCTION TO GENERATE A DAMPED SINE ATOM
% GIVEN AN INDEX NUMBER
function  atom  = atomgen(ind)

T = 0.01;
winlen = 2048;
bin = mod(ind-1,2048);                        
bl = ceil(ind/2048);
l = 2048-bl+1;
a = log(T)/l;
normal = sqrt((1-a^2)/(1-(a^2)^l));
pad = zeros(bl-1,1);
                
win = [pad;exp(a*[0:l-1]')];
win = win/norm(win);
                
sup = (0:winlen-1)';
               
            
atom = win .* exp(1j*2*pi/2048 * bin * sup);
end


