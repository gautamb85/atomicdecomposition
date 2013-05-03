function [book subset] = MPshortfastgab(y,winlen,hop,fftlen,maxiter)

M = length(y);
y = y(:);
nstride = ceil((M-winlen)/hop);
disp = 10;
booklen = fftlen*(nstride +1) ;

book = zeros(booklen,1);
X = zeros(booklen,1);
srr = 0;

res = y;
pow = 10*log10(res'*res);

%Analysis Window
window = hanning(winlen);
window = window/norm(window);

subset = [];
% Conjugate subspace setup

scal = [1/sqrt(2); ones(fftlen/2-1,1); 1/sqrt(2); zeros(fftlen/2-1,1)];
scalar = repmat(scal,nstride+1,1);
            

wincorr = fft(window.^2,fftlen);
c = [0; wincorr(3:2:end-1); 0; wincorr(3:2:end-1)];            
gamma = repmat(c,nstride+1,1); 



 %res = res/norm(res);
 its = 0;
 
 X = shortgabcorr(res,window,hop,fftlen,M);

 
 while(its<maxiter)
  
    a = (X - gamma.*conj(X))./(1 - abs(gamma').^2)';
    Theta = 2*real(a .* conj(X));
    [~, ind] = max(Theta.*scalar);
     
    conj_ind = conj_index(ind,fftlen);
    alpha = a(ind);
    
    [~,ind] = max(abs(X));
    alpha = X(ind);
        
    book(ind) = book(ind) + alpha;
    
    [atom asup] = gabtrans_atomgen(ind,window,hop,fftlen);
    
%     
    if ind == conj_ind % dc or nyquist atom
            
            % record projection
            book(ind) = book(ind) + alpha;
            subset = [subset; ind]; % record index
            %corr = [corr; alpha]; % record weig
            dx = gabrep_kern(real(alpha*atom),asup(1),asup(end),window,hop,fftlen,M);
    
    else % conj atoms
            
            conj_alpha = conj(alpha);
            
            % record projection
            book(ind) = book(ind) + alpha;
            book(conj_ind) = book(conj_ind) + conj_alpha;
            subset = [subset; ind; conj_ind]; % record index

             dx = crosscorr(2*real(alpha*atom),asup(1),asup(end),window,hop,fftlen,M);
    end
            %update inner products
            X = X - dx;
    
    
             its = its+1;
    
%     y_hat = real(gabtrans_synthesize(book,winlen,hop,fftlen,M));
%     res = y - y_hat;
%     srr = pow - 10*log10(res'*res);
    
    %if(mod(its,disp)==0)
        fprintf('atom selected = %d \n',ind);         
end
    
end
 
