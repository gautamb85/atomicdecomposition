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

clear all;
clc;

a = 0.666*cos(2*pi*666*[0:8192*4-1]'/11000);
l = 8192*4;
a1 = cos(2*pi*333*[0:2047]'/11000);
p1 = zeros(l-2048,1);
a1 = [a1;p1];

a2 = 0.533*cos(2*pi*145*[0:4095]'/11000);
p2 = zeros(l-4096,1);
a2  = [a2;p2];


ip = a1 + a2 + a;
%do signal postpadding

%initialize residual
res = ip;

lip = length(ip);

sc = 8192;

fftlen = 8192;


C = [];



    % find initial seed atom
        
    scale = sc;
    
    window = hanning(scale);
    window = window/norm(window);
    
    hop = scale/2;
    
    %no. of hops
    nhop = ceil((lip-scale)/hop);
 
 %function to calculate the correlations for seed atom determination
 function C = seedcorr(res,lip,scale)
      
    hop = scale/2;
    
    l = length(res);
    %no. of hops
    nhop = ceil((l-scale)/hop);
 
    
    for ptr = 0:nhop,
        
        
        samp_ind = 1 + ptr*hop;
        corr_ind = 1+ ptr*fftlen;
        
        %for fftlen = scale
        %fcorr = fft(window.*ip(samp_ind ; samp_ind + scale),scale);
        
        fcorr = fft(window.*res(samp_ind : samp_ind + scale-1),fftlen);
        
        C(corr_ind:corr_ind+fftlen-1,1) = fcorr;
        
    end
  
 end
    
    %maximum absoulute value of correlation
    % index of seed atom
    [mms,ind] = max(abs(C));
      
     alpha = C(ind);
    
    % store values in a book
    
    %stblk is the amount of zero padding
    % ind1 = 8192*5+30;
     [atom stblk bin] = genscaleatom(ind1,scale,fftlen);
    
     % pad this atom appropriately and subtract from
     %residual
     prepad = zeros(stblk,1);
     po = length(ip)-length(prepad) - length(atom);
     popad = zeros(po,1);    
     %pad atom upto input length
     fatom = [prepad;atom;popad];
     
     
     %subtract this atom from the residual
     res = res - real(alpha.*fatom);
     
     % determine analysis window to consider for the molecule
     dicblk = overlap(stblk,hop,lip);
     
     %generate explicit dictionary
     D = dictgen(ip,dicblk,bin);
     
     %start building a molecule with this seed atom
      %moleculeMP(D,res)
     
     
     
     
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
     fbat = 3;
     
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
 
 function D = dictgen(lip,blks,bin)

%generate explicit dictionary

     % inputs is a vector of blocks overlapping in time
     % with the seed atom
 
     k = length(dicblk);
     D = [];

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
 
 function moleculeMP(D,res,dicblk)
 
 r =size(D,2);
 book2 = zeros(fftlen*r,1);
 
 Cs = D'*res;
 lipr = length(res);
 
 [mma,in3] = max(abs(Cs));
  
 %stop the molecule and search for a new seed atom
 % based on the correlation of each atom selected in a 
 %perticular molecule
 
 if(mma >0.3*mms)
     
     alpha2 = Cs(in3);

     %book2(in3) = book2(in3) + alpha;

    % function molatomgen(in3)
    
    %analysis block this window belongs to
     st1 = dicblk(in3);

     % block relative to the original signal
      
      %stblk1 = (act_blk)*hop;
     
     % index relative to main book
     actindex = (st1-1)*fftlen + bin;
     
     
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
     poe = zeros((length(ip) - length(pp) - scale),1);
      
     gatom = [pp;gatom;poe];
     res = res - real(alpha2*gatom);
     
      %update main book 
      book(actindex) = alpha2;
      
     % determine backward and forward overlapping atoms based on the
     % currently selected atom and repeat the MP process 
      dicblk2 = overlap(stblkg,hop,lipr);
     %generate explicit dictionary based on this block of atoms
      Dm = dictgen(lipr,dicblk2,bing);

      moleculeMP(D,res,dicblk);

 else
     
     break;

 end
 
 
 
     
     
 
        







    
     
     
     
     
     
     
     
     
     
     
     
         
         
     
     
    
    
    
    
    
   