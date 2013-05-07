clear all;

%i/p signal
x = rand(33563,1);

 l = 64;
 winlen = 2048;
 win = hanning(l);
 
 %hop size
 hop = l/2;
 
 %preprocessing signal
 %padding with a hop size of zeros
 pad = zeros(hop,1);
 %postpad so that signal is a multiple of hop length
 l = hop + length(x);
 post = ceil(l/hop)*hop - l;
 postpad = zeros(post,1);
 
 ip = [pad;x;postpad];
 
 nhop = length(ip)/hop -1;
 
 X = [];
 mid = hop;
 
 for i = 0:nhop-1
     
     X(:,i+1) = fft(win.*ip(mid-hop+1:mid+hop),2048);
     mid = mid+hop;
 end
 
 %% 
 [s fs] = wavread('GLOCK_A5.wav');

s = s(:,1);


samp = s(1:2048);
     
%%
T = 0.15;
a = log(T)/64;
win = exp(a.*[0:63]');

d = win.*exp(1j*2*pi*50/2048*[0:63]');
d = [zeros(1024,1);d;zeros(2048-1088,1)];


[x] = damptrans(d,64,32,2048);
[~,ior] = max(abs(x));

%%
T = 0.01;
winlen = 2048;
pos = 666;
atom = [];
% forward & backward lengths
bl = pos;
fl = winlen - bl;


af = log(T)/fl;
ab = log(T)/bl;

% damped exponential window - forward & backward
fwin = exp(af.*[0:fl-1]');
fwin = fwin/norm(fwin);

bwin = exp(-1*ab*[0:bl-1]');
bwin = bwin/norm(bwin);

%forward and backward atoms
fatom = fwin.*exp(1j*2*pi*50/2048*[0:fl-1]');

batom = bwin.*exp(1j*2*pi*50/2048*[0:bl-1]');

atom(1:bl) = batom;

atom(bl+1:2048) = fatom;

























 