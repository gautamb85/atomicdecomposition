%FUNCTION TO GENERATE A DAMPED SINE ATOM
% GIVEN AN INDEX NUMBER
function  atom  = gab_atomgen(ind)
N = 2048;
winlen = 2048;
bin = mod(ind-1,2048);                        
                
win = hanning(winlen);
win = win/norm(win);
                
sup = (0:winlen-1)';
               
            
atom = win .* exp(1j*2*pi/N * bin * sup);

end


