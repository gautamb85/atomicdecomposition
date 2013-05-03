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

lip = length(ip);

sc = [8192,4096,2048,1024,512];

fftlen = 8192;

%initialize correlation vector
c = [];
C = [];

for i = 1:5,

    % find initial seed atom
        
    scale = sc(i);
    
    window = hanning(scale);
    window = window/norm(window);
    
    hop = scale/2;
    
    nhop(i) = ceil((lip-scale)/hop);
    
    
    
    for ptr = 0:nhop(i),
        
        
        samp_ind = 1 + ptr*hop;
        corr_ind = 1+ ptr*fftlen;
        
        %for fftlen = scale
        %fcorr = fft(window.*ip(samp_ind ; samp_ind + scale),scale);
        
        fcorr = fft(window.*ip(samp_ind : samp_ind + scale-1),fftlen);
        
        c(corr_ind:corr_ind+fftlen-1,1) = fcorr;
        
    end
    C = [C;c];
    
    % block divisions
    bldiv(i) = length(c);
    
end

%%  find seed atom = working

[~,ind] = max(abs(C));
k = length(bldiv);

%make a matrix of correlations

    
     C1 = C(1:bldiv(1));
     
     C2 = C(bldiv(1)+1:bldiv(1) + bldiv(2));
     
     C3 = C(bldiv(2)+1:bldiv(2)+bldiv(3));

     C4 = C(bldiv(3)+1:bldiv(3)+bldiv(4));

     C5 = C(bldiv(4)+1:bldiv(4)+bldiv(5));

     
     if (ind<=bldiv(1))
        selscale = sc(1);
        [atom stblk asup] = genscaleatom(ind,selscale,fftlen);

     elseif(ind<=bldiv(1)+bldiv(2))
        selscale = sc(2);
        in = ind - bldiv(1);
        [atom stblk asup] = genscaleatom(in,selscale,fftlen);

     elseif(ind<=bldiv(1)+bldiv(2)+bldiv(3))
        selscale = sc(3);
        in = ind - bldiv(1)-bldiv(2);
        [atom stblk asup] = genscaleatom(in,selscale,fftlen);

     elseif(ind<=bldiv(1)+bldiv(2)+bldiv(3)+bldiv(4))
        selscale = sc(4);
        in = ind - bldiv(1)-bldiv(2) - bldiv(3);
        [atom stblk asup] = genscaleatom(in,selscale,fftlen);


     elseif(ind<=bldiv(5))
        selscale = sc(5);
        in = ind - bldiv(1)-bldiv(2) - bldiv(3)-bldiv(4);
        [atom stblk asup] = genscaleatom(in,selscale,fftlen);
        
     end

%% create a new subspace in order to track partial


        

        
        
        

