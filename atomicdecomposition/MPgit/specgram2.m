function [] = specgram2(x, nfft, fs, wind, noverlap, scale)
%SPECGRAM2(X, NFFT, FS, WIND, NOVERLAP, SCALE, LIMITS) produces a
%dynamic-limited (linear or Bark) spectrogram of the signal X.
%
%   X is the input signal at sampling frequency FS. (Necessary args.)
%
%   The window is given as a length of a hanning window (scalar) or
%   as a window vector. Overlap defaults to (50%), but is the number
%   of samples to overlap.
%  
%   NFFT is the number of discrete-frequency points. SCALE indicates
%   'linear' analysis or 'Bark' (increased low-frequency resolution)
%   analysis. Can be left empty as in [], which defaults to linear
%   scale.
%   
%   LIMITS is an length-2 vector holding the dynamic range limits. It
%   defaults to -20dB as upper limit and -80dB as lower limit, 
%   [20 80].
%   
%   In order to ease comparison of several plots, spurious peaks
%   artifacts are removed by using an absolute signal maximum
%   (approx. 24dB) relative to a variance-scaled signal (with a
%   std. deviation of .1).
%   
%   If you need to use these plots in a subplot, do e.g.:
%   $ figure; subplot(211);
%   $ specgram2(x1, 1024, fs);
%   $ subplot(212);
%   $ specgram2(x2, 1024, fs);
%   
%   An example of advanced use:
%   $ specgram2(x, 1024, fs, 240, 230, [], [30 90]);
%   
%   Resolution issues. Reasonable results (fast) can be obtained by
%   Hanning tappering 160 samples with 4:1 overlap. To enhance frequency
%   resolution 240 samples can be used while time-compensating using 230
%   samples overlap.
   
% Author: Martin Hammer, 2006-12-05, 2007-01-15
%
% $Id: specgram2.m 481 2007-04-12 15:23:17Z mhje03 $

%    msg = nargchk(4, 7, nargin);
%    error(msg);
% 
% 
%    if (nargin < 6)
%       scale = []; % default linear scale
%    end
%    if (nargin < 7),
%       maxdB = -20; % max. and min. spectral value
%       mindB = -80; 
%    else
%       maxdB = -limits(1);
%       mindB = -limits(2);
%    end
%       
%    if (nargin < 5),  noverlap = [];  end % default 50%
% 
%    if (isempty(scale)), scale = 'linear'; end; % default linear scale
%       
%    [scale, msg] = get_method(scale, {'linear', 'bark'});
%    error(msg);
   
    
   maxdB = -20; % max. and min. spectral value
   mindB = -80;
   
   nwind = length(wind);
   if ( 1 == nwind ),
      wind = hanning(wind); % use hanning window
      nwind = length(wind);
   end % otherwise wind is the window
   if( isempty(noverlap) ),
      noverlap = .5*nwind; % defaults to 50% overlap
   end
   
   x = x(:); % make a column vector for later ease
   wind = wind(:); % consistent window
   nx = length(x); % signal length

   % Variance scaling and absolute maksimum signal value in order to
   % reduce spurious peak artifacts. Signal is scaled to a
   % std. deviation of .1 and max. of approx. 24dB
   x = x*.1/std(x); % variance scaling (std. deviation = .1)
   % y_max = 20*log10(16); % set maximum to approx. 24dB (relative to the variance scaled signal). (See farther down)
  

   %%% Form signals, vectors and filter (DFT)
   segs = floor((nx-noverlap)/(nwind-noverlap)); % no. of segments to process
   index = 1:(nwind-noverlap):nx; % index for first element in each segement
   if ( nx < index(end)-1+nwind ),
      x(index(end)-1+nwind) = 0; % zero pad, if necessary
   end
   
   y = zeros(nwind, segs); % initiate output matrix
  
   if (strcmp(scale, 'bark')), % do segmentation of signal
      lambda = barkwarp(fs); % compute warping factor
      for k = 1:segs,
         y(:,k) = wind.*longchain(x(index(k):index(k)-1+nwind), nwind-1, lambda)';
         % y(:,k) = wind.*warp_impres(x(index(k):index(k)-1+nwind), lambda, nwind-1)';
      end
   else % scale == 'linear'
      for k = 1:segs,
         y(:,k) = wind.*x(index(k):index(k)-1+nwind);
      end
   end
   
   y = fft(y, nfft); % compute DFT of all segements
   if ( rem(nfft,2) ), % odd DFT
      y = y(1:(nfft+1)/2, :);
   else % even DFT
      y = y(1:(nfft/2+1), :);
   end
   [nfreq, segs] = size(y); % segs is not changed, retreive number of freq points

   
   %%% Limit scale to enhance signal information (human perception)
   y = 20*log10(abs(y)+eps);
   % y_max = max(max(y)) % find maximum value over all time and frequency segements (changed to absolute maximum relative variance-scaled signal).
   y_max = 20*log10(16); % absolute max, approx 24dB.
   maxval = y_max + maxdB; % maximum allowed spectrogram magnitude 
   minval = y_max + mindB;
   too_large = find(y > maxval); % extract values larger than maxval
   too_small = find (y < minval);
   y(too_large) = maxval; % set equal to maxval (make dynamic area smaller)
   y(too_small) = minval;

   
   %%% Plotting
   newplot;
   t = (index-1)/fs; % time index
  % if (strcmp(scale, 'bark')),
      nfreq_pts = 8; % use eight freqency points on y-axis
      F = freqz([-lambda 1], [1 -lambda], nfreq_pts); % extract warping curve
      a = -angle(F)/pi; % positive angle, 0 <=a <1
      freq_axis = (0:fs/2000/(nfreq_pts-1):fs/2000).*a'; % form freq. axis values
      freq_axis = round(freq_axis*100)/100;
      imagesc(t, 1:nfreq, y); % plot spectrogram
      set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', freq_axis, 'YTick', 1:nfreq/nfreq_pts:nfreq);
      title('Spectrogram (Bark scale)');
   %else % scale == 'linear'
    figure;
   imagesc(t, linspace(0, fs/2000, nfreq), y); % plot spectrogram
      title('Spectrogram (linear scale)')
   %end
   axis xy; % mirror frequency axis (Matlab's stupid DFT notation)
   colormap('jet');
   xlabel('Time [s]'); ylabel('Frequency [kHz]');

