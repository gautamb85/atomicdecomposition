clear all;
clc;

a = cos(2*pi*4000/8192*[0:8192*20-1]'/11000);
l = 8192*20;
a1 = cos(2*pi*1200*[0:2047]'/11000);
p1 = zeros(l-2048,1);
a1 = [a1;p1];

a2 = 0.533*cos(2*pi*145/8192*[0:8192*20-1]'/11000);



ip = a2 + a+a1;
%do signal postpadding

window = hanning(8192);
C = dtest(ip,window,8192,8192);

D = GaborBlock(window, 4096, 8192, length(ip));

[book, snr, lpn, subset, corr, l0, cpu] = MP(D,ip,5,20,10);
%%
lip = length(ip);

window = hanning(8192);

book =  MolecularMP(8192,ip,8192,window);


%%

wina = window/norm(window);

h = find(book);
op = zeros(length(ip),1);

for i = 1:length(h),
    
    b = mod(h(i)-1,8192);
    sh = floor(h(i)/8192)*4096;


    pre = zeros(sh,1);
    po = length(ip) - length(pre) -8192;
    popad = zeros(po,1);
    
    atom = wina.*exp(1j*2*pi*b/8192*[0:8191]');
    fatom = [pre;atom;popad];
    b=0;
    sh=0;
    
    op = op + fatom;
end
%%
subplot(211);plot(real(ip));subplot(212);plot(real(op));
    
    