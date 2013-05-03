clear all;
clc;

[s fs] = wavread('GLOCK_A5.wav');

s = s(:,1);


samp = s(1:2048*12);
    

% %try with and without window
% win = hanning(2048);
% win = win/norm(win);
% 
% %define a dictionary for this-time frequency space
% sn = samp(1:2048*12);
% 
% % spectrum
% S = fft(sn,2048);

%%

%gabor window analysis
M = length(samp);
[book1 subset] = MPshortfastgab(samp,2048,512,2048,50);
%%
% wina = hanning(2048);
% D = GaborBlock(wina,512,2048,length(samp));
% [book snr lpn sub]= MP(D,samp,50,10,10);
% %%
win = hanning(2048);
win = win/norm(win);
op = real(gabsynthesize(book1,win,512,2048,M));
 plot(op);
%%
% %gabor residual
% gres = sn - op1;
% 
% % transient modelling with short windows 
% [book1r op1r] = MPshortgab(gres,64,16,2048,100);
% gnoise = gres - op1r;
% 
% subplot(511);plot(sn);title('Original sound');
% subplot(512);plot(op1);title('Tonal sound');
% subplot(513);plot(gres);title('Residual sound');
% subplot(514);plot(op1r);title('Transient sound');
% subplot(515);plot(gnoise);title('Residual Noise');xlabel('Gabor Analysis');
% 
% 
% % Damped Sinusoid Analysis
% [book2 srr2 op2] = MPdamped(sn,2048,1024,2048,100);
% dsres = sn - op2;
% 
% [book2r srr2r op2r] = MPdamped(dsres,64,16,2048,100);
% snoise = dsres - op2r;
% 
% figure;
% subplot(511);plot(sn);title('Original sound');
% subplot(512);plot(op2);title('Tonal sound');
% subplot(513);plot(dsres);title('Residual sound');
% subplot(514);plot(op2r);title('Transient sound');
% subplot(515);plot(snoise);title('Residual Noise');xlabel('Damped Sinusoid Analysis');
% 
% 
% 
% % Damped Heaviside analysis
% 
% [book3 subset3 corr3 srr3 op3] = MPheavy(sn,2048,50);
% 
% resh = sn - op3;
% % transient modelling
% [book3r subset3r corr3r srr3r op3r] = MPheavy(resh,2048,50);
% dsnoise = resh - op3r;
% 
% figure;
% subplot(511);plot(sn);title('Original sound');
% subplot(512);plot(op3);title('First reconstructed sound');
% subplot(513);plot(resh);title('Residual sound');
% subplot(514);plot(op3r);title(' Reconstructed Residual');
% subplot(515);plot(dsnoise);title('Residual Noise');xlabel('Damped Heaviside Analysis');
% 
% 
% 
% 
% 
% %%
% % op1 = synthesize(book,'f');
% % op2 = synthesize(book,'t');
% % 
% % %residual
% % res = sn - real(op2);
% % R = fft(res);
% % %spectral
% % sres = S-op1;
% % %frequency axis
% % f = -fs/2:fs/2048:fs/2-fs/2048;
% % 
% % subplot(211);plot(real(op2));title('tonal part');
% % subplot(212);plot(res);title('residual');
% % 
% % figure;
% % subplot(411);plot(f,abs(S));title('Spectrum');
% % subplot(412);plot(f,abs(op1));title('Reconstucted tonal spectrum');
% % subplot(413);plot(f,abs(sres));title('Difference');
% % subplot(414);plot(f,abs(R));title('Residual spectrum');
% % 
% % 
% % %%
% % figure;
% % [book2 sd co srre  op12] =  MPheavy(sn,2048,100);
% % plot(op12)
% % %% decompose residual using damped sinusoid transient block
% % 
% % %
% % 
% % plot(op);
% % %%
% % [booka subset corr srr output] = MPheavy(res,2048,50);
% % %%
% % [book3 op2] = MPdamped(res,64,32,2048,50);
% % 
% % %%
% % subplot(311);plot(op);title('Gabor atoms (64)'); subplot(312);plot(op1);title('Heaviside damped atoms');subplot(313);plot((op2));title('Damped sinusoid atoms (64)');xlabel('No. of atoms = 50, hop = 32 samples');
% % %%
% % saveas(gcf,'\Users\gautambhattacharya\Desktop\transient2.jpg','jpg')
% % 
% % %%
% % 
% % l = length(find(book));
% % l1 = length(find(book1));
% % l2 = length(find(book3));
% % 
% % 
% % 
% % 
% %  