% clear all;
% 
% %i/p signal
% x = rand(33563,1);
% 
%  l = 64;
%  winlen = 2048;
%  win = hanning(l);
%  
%  %hop size
%  hop = l/2;
%  
%  %preprocessing signal
%  %padding with a hop size of zeros
%  pad = zeros(hop,1);
%  %postpad so that signal is a multiple of hop length
%  l = hop + length(x);
%  post = ceil(l/hop)*hop - l;
%  postpad = zeros(post,1);
%  
%  ip = [pad;x;postpad];
%  
%  nhop = length(ip)/hop -1;
%  
%  X = [];
%  mid = hop;
%  
%  for i = 0:nhop-1
%      
%      X(:,i+1) = fft(win.*ip(mid-hop+1:mid+hop),2048);
%      mid = mid+hop;
%  end
%  
%  %% 
%  [s fs] = wavread('GLOCK_A5.wav');
% 
% s = s(:,1);
% 
% 
% samp = s(1:2048);
%      
% %%
% T = 0.15;
% a = log(T)/64;
% win = exp(a.*[0:63]');
% 
% d = win.*exp(1j*2*pi*50/2048*[0:63]');
% d = [zeros(1024,1);d;zeros(2048-1088,1)];
% 
% 
% [x] = damptrans(d,64,32,2048);
% [~,ior] = max(abs(x));
% 
% %%
% Tf = 0.01;
% Tb = 0.01;
% winlen = 2048;
% pos = 666;
% atom = [];
% % forward & backward lengths
% bl = pos;
% fl = winlen - bl;
% 
% 
% af = log(Tf)/1024;
% ab = log(Tb)/1024;
% 
% % damped exponential window - forward & backward
% fwin = exp(af.*[0:1024-1]');
% %fwin = fwin/norm(fwin);
% 
% fpad = zeros(1024,1);
% fwin = [fpad;fwin];
% fwin = fwin/norm(fwin);
% %plot(real(fwin));
% %figure;
% bwin = exp(-1*ab*[0:1024-1]');
% bwin = bwin/norm(bwin);
% bpad = zeros(1024,1);
% bwin = [bwin;bpad];
% 
% win =fwin+bwin;
% 
% 
% 
% %forward and backward atoms
% atom = win.*exp(1j*2*pi*50/2048*[0:winlen-1]');
% %plot(real(atom));
% %batom = bwin.*exp(1j*2*pi*50/2048*[0:bl-1]');
% 
% 
 %% 
clear;

[x fs] = wavread('GLOCK_A5.wav');
x = x(:,1);
samp = x(1:12*2048);

%specgram2(l, 1024, fs, 240, 230,'bark');
   


% modes for bark spectrogram

%mode = 'e''m'or 'b' - 'f' for in hz


%[g gsup]= gabtrans_atomgen(7700,64);


wina = hanning(2048);
D = GaborBlock(wina,512,2048,length(samp));
book = MP(D,samp,100,10,10);
%%
oop = D*book;
plot(real(oop));
figure;
%%
 we = hanning(2048);
 we = we/norm(we);
 opi = gabsynthesize(book,we,512,2048);
 plot(real(opi));
% %%
% [at asup] = D.get_atom(7700);
% s = D.nstrides;
% 
% 
% l1 = D'*samp;
% l2 = shortgabcorr(samp,2048,512,2048,length(samp));
% [at asup] = D.get_atom(379);
% [gt gsup] = gabtrans_atomgen(379,2048);
% vv = asup - gsup;
% v = at - gt;
% %%
% XC = D.rep_kern(at,asup(1),asup(end));
% X = crosscorr(gt, asup(1), asup(end),2048,512,2048,length(samp));
% 
% f = length(find(X-XC));
% 
% f1 = length(find(l1-l2));
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 








 