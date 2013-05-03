%% want to do a structured decomposition

%start with stft analysis
% need to keep track of this stuff
% we assume that we are modelling the sustain part
%of a sound
% 
% fs = 11000;
% fftlen = 4096;
% winlen = 4096;
% hop = winlen/8;
% 
% 
% ip = cos(2*pi*440*[0:9999]'/fs);

function X = gtbSTFT(ip,winlen,hop,fftlen)

    window = blackman(winlen);
    window = fftshift(window);
    window = window/norm(window);

    l = length(ip);
    % pad the input to a multiple of hop/window
    p = ceil(l/hop)*hop-l;
    pad = zeros(p,1);
    ipp = [ip;pad];

    %number of hops
    nhop = ceil((length(ipp)-winlen)/hop) ;

    % STFT analysis
    mid = winlen/2;
    pf = mid/hop;

    % initialize STFT buffer
    X = zeros(fftlen,nhop);

    for i = 0:nhop-1,
        chunk = fft(window.*ipp(mid - pf*hop +1:mid + pf*hop),fftlen);
        % normalized spectrum
        X(:,i+1) = chunk/norm(chunk);
    end

end

