function conj_ind = conj_index(ind,fftlen)
           
bin = mod(ind-1,fftlen);
conj_bin = mod(fftlen-bin, fftlen);
conj_ind = conj_bin + floor(ind/fftlen)*fftlen + 1;

end