%Matching Pursuit
%Algorithm :
% 1. Calculte correlations btw sound and dictionary

% 2. Select the abs max value, create a seach space of atoms of same
%    frequency, at all time steps
% 3. Build a molecule on this basis.
% 4. Figure out how to stop
% 5. Orthogonalization of coefficients
% 
% Scales considered - 8192,4096,2048,1024,512 complex gabor atoms
% hop for all scales is scale/2

%adding a special comment!!!!





function mbook = MolecularMP(sc,ip,fftlen,window)

    res = ip;
    res = res/norm(res);
    lip = length(ip);
    % find initial seed atom
        
    % power
      pow = 10*log10(res'*res);
    scale = sc;
    

    window = window/norm(window);
    
    hop = scale/2;
    
    %no. of hops
    nhop = ceil((lip-scale)/hop);
    
    %initialize book (check gablab)
    book = zeros(fftlen*(nhop+1),1);
    mbook = zeros(fftlen*(nhop+1),1);
    
    % find correlation between signal and entire gabor dictionary
    C = seedcorr(res,window,scale,fftlen);
     
    %maximum absoulute value of correlation
    % index of seed atom
    [mms,ind] = max(abs(C));
      
     alpha = C(ind);
     
     % store seed atom in book. This should work such that i can synthesize
     % the first molecule
     
     book(ind) = book(ind) + alpha;
    
    
    %stblk is the amount of zero padding
    % ind1 = 8192*5+30;
     [atom stblk bin] = genscaleatom(ind,scale,fftlen);
    
     % pad this atom appropriately and subtract from
     %residual
     prepad = zeros(stblk,1);
     po = length(ip)-length(prepad) - length(atom);
     popad = zeros(po,1);    
     %pad atom upto input length
     fatom = [prepad;atom;popad];
     fatom = fatom/norm(fatom);
     
     
     %subtract this atom from the residual
     res = res - real(alpha.*fatom);
     
     % calculate srr (based on seed atom)
      pres = 10*log10(res'*res);
      srr = pow -pres;
      fprintf('SRR acheived by seed atom = %d \n',srr);   

     
     % determine analysis window to consider for the molecule
     dicblk = overlap(stblk,hop,lip);
     
     %generate explicit dictionary
     D = dictgen(lip,scale,dicblk,window,bin);
     
     %start building a molecule with this seed atom
      mbook = buildmolecule(D,res,window,book,bin,dicblk,scale,fftlen,mms,0,pres);
end
  
  %function to calculate the correlations for seed atom determination
 function C = seedcorr(res,window,scale,fftlen)
      
    hop = scale/2;
    
    l = length(res);
    %no. of hops
    nhop = ceil((l-scale)/hop);
    nblocks = ceil(l/hop);
    C =[];
    
    for ptr = 0:nhop,
        
        
        samp_ind = 1 + ptr*hop;
       
        corr_ind = 1+ ptr*fftlen;
        
        %for fftlen = scale
        %fcorr = fft(window.*ip(samp_ind ; samp_ind + scale),scale);
        
        fcorr = fft(window.*res(samp_ind : samp_ind + scale-1),fftlen);
        
        C(corr_ind:corr_ind+fftlen-1,1) = fcorr;
        
    end
  
 end    
     
     
     
% modify the dictionay

% function to find the number of ovelapping blocks
 
function dicblk = overlap(stblk,hop,lip)
   
     % the selected analysis window block in time  
     ablok = floor(stblk/hop)+1;
     %block 1 = no zero pad, block 2 = 1hop zize zero pad
     
     
     % total number of analysis blocks
     totblk = ceil(lip/hop)-1;
     
     fblk = 0;
     sblk = 0;
     
     %fwd and bwd atoms, make it a parameter
     fbat = 5;
     
     for j=1:fbat,
        
         if(ablok + j > totblk )
             
             fblk = totblk;
             sblk = ablok-j;
         
         elseif(ablok-j<1)
             
             fblk = ablok + j;
             sblk = 1;
             
         else
             sblk = ablok - j;
             fblk = ablok + j;
         end
     end
     
     dicblk = sblk:fblk; 
     
end
     
    
% correlation with all ovelapping blocks at all frequencies

%      
%      dicblk = sblk:fblk;
%         
%      s = length(dicblk);
%      Corr1 = [];
% 
%     for j=1:s-1,
%         
%        samp_ind = 1+(dicblk(j)-1)*hop;
%        fft_ind = 1+j*fftlen;
%        
%        chunk1 = window.*res(samp_ind:samp_ind + fftlen-1);
%        
%        Corr1(fft_ind:fft_ind+fftlen-1,1) = fft(chunk1,fftlen);
%     end
    
%%
 %%[~,ind1] = max(abs(Corr1));
 

 
 %function to generate an atom using block no. and
 %bin - right now we are only considering the same frequency throughout  
 
 function D = dictgen(lip,scale,dicblk,window,bin)

%generate explicit dictionary
        
     % inputs is a vector of blocks overlapping in time
     % with the seed atom
 
     k = length(dicblk);
     D = [];
    
     %remember to change this
     hop = scale/2; 
     
     for g = 1:k,

         pre = (dicblk(g)-1)*hop;
         prepad = zeros(pre,1);

         pos = lip - length(prepad) - scale;
         postpad = zeros(pos,1);

         atom = window.*exp(1j*2*pi*bin/scale*[0:scale-1]');
         ratom = [prepad;atom;postpad];

         D(:,g) = ratom/norm(ratom);
     end
 
end
 %function to perform matching pursuit with this explicit dictionary
 
 %analysis
 %%
 
 function [book done] = buildmolecule(D,res,window,book,bin,dicblk,scale,fftlen,mms,done,pres)
 
 r =size(D,2);
 %book2 = zeros(fftlen*r,1);
 
 Cs = D'*res;
 lipr = length(res);
 hop = scale/2;
 
 
 [mma,in3] = max(abs(Cs));
  
 %counter for the number of atoms in a molecule
 %apart from the seed atom
 c = 0;
 
 %stop the molecule and search for a new seed atom
 % based on the correlation of each atom selected in a 
 %perticular molecule
 %done =0;
 
 while (done==0)
     
     if(mma <0.1*mms)
         done =1;
        

     else

         alpha2 = Cs(in3);

         %book2(in3) = book2(in3) + alpha;

        % function molatomgen(in3)

        %analysis block this window belongs to
         st1 = dicblk(in3);

         % block relative to the original signal

          %stblk1 = (act_blk)*hop;
         if (st1==1)
             actindex = 1 + bin;
         % index relative to main book
         else
            actindex = (st1-1)*fftlen + bin;
         end

          % generate the atom

    %       pre1 = (st1-1)*hop;
    %       pppad = zeros(pre1,1);
    % 
    %       poi = length(ip) - scale - length(pppad);
    %       poipad = zeros(poi,1);
    % 
    %       at = window.*exp(1j*2*pi*bin/scale*[0:scale-1]');
    %       tom = [pppad;at;poipad];

         [gatom stblkg bing] = genscaleatom(actindex,scale,fftlen);

         %pad the atom an subtract it from the residual
         pp = zeros(stblkg,1);
         poe = zeros((lipr - length(pp) - scale),1);

         gatom = [pp;gatom;poe];
         %normailize all atoms and ip
         gatom = gatom/norm(gatom);
         
         res = res - real(alpha2*gatom);
         presm = 10*log10(res'*res);
         % SRR
         srr = pres -  presm;
          fprintf('SRR = %d \n',srr);   
          
          if (srr >= 1)
              done=1;
          end

          %update main book 
          book(actindex) = book(actindex) + alpha2;

         % determine backward and forward overlapping atoms based on the
         % currently selected atom and repeat the MP process 
          dicblk2 = overlap(stblkg,hop,lipr);
         %generate explicit dictionary based on this block of atoms
          Dm = dictgen(lipr,scale,dicblk2,window,bing);
         %run matching pursut and repeat until the correlation of the selected
         % atom is too low or 5 atoms are selected 
          [book done] = buildmolecule(Dm,res,window,book,bin,dicblk2,scale,fftlen,mms,done,presm);
         
          
%           c = c +1;
%           if (c ==5)
%               done =1;
%           end     

     end
    
 end
 
 end
 
 
 
     
     
 
        







    
     
     
     
     
     
     
     
     
     
     
     
         
         
     
     
    
    
    
    
    
   