function op = dsheavy_synth(book)        

 Tf = 0.15;
 Tb = 0.15;
 winlen = 2048;
 fl = winlen/2;
 bl = fl;
 af = log(Tf)/fl;
 ab = log(Tb)/bl;
 
 winlen = 2048;    
 fftlen = 2048;
            
  book = book(:); % make row vector
  blocklen = floor(length(book)/fftlen)-1;
  op = 0; % init. output buffer
  
  for ptr =0:blocklen              
                
                %postpad = zeros(fftlen - length(prepad)-l,1);
                
                fwin = af.*exp(af*[0:fl-1]');
                fwin = fwin/norm(fwin);
                fpad = zeros(1024,1);
                fwin = [fpad;fwin];
                
                bwin = ab.*exp(-ab*[0:bl-1]');
                bwin = bwin/norm(bwin);
                bpad = zeros(bl,1);
                bwin = [bwin;bpad];
                
                window = fwin + bwin;
                window = window/norm(window);
                
                fft_ind = 1 + ptr*fftlen;
                fft_sup = fft_ind:fft_ind+fftlen-1;
                block = book(fft_sup);
                
                chunk = ifft(block, fftlen) * fftlen;
           
                op = op + window .* chunk; % mult synth win
                
            end