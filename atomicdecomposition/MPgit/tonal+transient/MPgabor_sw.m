function [book1] = MPgabor_sw(a,winlen,fftlen,maxiter)
%add variable time shifts
%stopping condition based on number of atoms
%modify to genetate either spectral or temporal reconstruction

%window = window/norm(window);
%booklen = N*(N - l +1) ;
%booklen = 2048*(2048-63);
%booklen = fftlen*(fftlen - winlen +1);
%initialize book
book1 = zeros(4065280,1);
%initialize correlation vector
%X = zeros(booklen,1);
%initialize residual
res = a;

its = 1;
 
X = fftcorra(res,64,2048);
 
while(its<=maxiter)
  
    [~,ind] = max(abs(X));
    alpha = X(ind);
    
    fprintf('atom selected = %d \n',ind);         
    
    book1(ind) = alpha;
    
    [atom asup] = atomgen(ind,64);
    
    dx = wincorra(real(alpha*atom),asup(1),asup(end),64,2048);
    
    X = X - dx;
    its = its+1;
   
end

%output = synthesize(book,winlen,fftlen,'t');
end


%FUNCTION TO GENERATE GABOR ATOM OF DESIRED LENGTH
function  [atom atomsup]  = atomgen(ind,winlen4)
%add functionality to input any window
bin = mod(ind-1,2048);                        
bl = floor(ind/2048);     
window = hanning(winlen4) ;
sup = (0:winlen4-1)';
atomsup = bl+sup;               
            
atom = window .* exp(1j*2*pi/2048 * bin * sup);
end

%calculate signal-dictionary correlations
%using FFT
function Z = fftcorra(x,winlen3,fftlen3) 

 %window parameters
 window = hanning(winlen3);
 window = window/norm(window);
 nstrides = fftlen3 - winlen3+1;
 
 Z = zeros(nstrides*fftlen3,1);
 
 for ptr = 0:nstrides-1,
         
     samp_ind = 1 + ptr;
     fftind = 1 +ptr*fftlen3;
     
     chunk = window .* x(samp_ind:samp_ind+winlen3-1);
    
     Z(fftind:fftind+fftlen3-1) = fft(chunk,fftlen3);
          
 end
                
end

%function to calculate at correlation changes
function C = wincorra(atom,st,ed,winlen2,fftlen2)
            
            window = hanning(winlen2);
            window = window/norm(window);
            
            st_block = max(0, st - winlen2);
            ed_block = min(ed,(fftlen2-winlen2));
            blocklen = ed_block-st_block;
       
            pre_pad = zeros(st - st_block, 1);
            post_pad = zeros(ed_block - length(pre_pad) + winlen2-1, 1);                        
            y = [pre_pad; atom; post_pad];            
            
             C = zeros(2048*(2048-63),1); % init. output buffer
             
             
             %C = zeros(booklen2,1); % init. output buffer

            for ptr = 0:blocklen-1,
                
                samp_ind = ptr + 1;
                fft_ind = (st_block + ptr)*fftlen2 + 1;
                
                chunk = window .* y(samp_ind : samp_ind+ winlen2-1);
                C(fft_ind : fft_ind+fftlen2-1) = fft(chunk,fftlen2);

           end   
       
end

function op = synthesize(book1,winlen1,fftlen1,ch)        

        
       %currently trying to synthesize o/p using gabor atoms only
       %experiment with changing the synthesis window    
            window = hanning(winlen1);           
            book1 = book1(:); % make row vector
            blocklen = floor(length(book1)/fftlen1)-1;
            
            op = 0; % init. output buffer
            
            for ptr =0:blocklen
                
                prepad = zeros(ptr,1);
                postpad = zeros(fftlen - length(prepad)-l,1);
                at = [prepad;window;postpad];
                at = at/norm(at);
                
                fft_ind = 1 + ptr*fftlen;
                fft_sup = fft_ind:fft_ind+fftlen1-1;
                block = book1(fft_sup);
                
                chunk = ifft(block, fftlen1) * fftlen1;
                
                if(ch=='t')
                    op = op + at .* chunk;
                elseif(ch=='f')
                    op = op + fft(at .* chunk,fftlen); 
                end
                
            end
end
            

