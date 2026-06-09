function [cleanECGchan, ecgVec] = pickCleanECG(data)
% pickCleanECG   Pick the cleaner ECG channel via QRS‐band SNR
% Only for ECG1 and ECG2.
% Check function pickCleanECG2 if you want to include biopac as well

% [cleanECGchan, ecgVec] = pickCleanECG(data)
%   data       : your ft-style struct (with .label, .trial, .fsample)
%   cleanECGchan   : name of the chosen channel (string)
%   ecgVec         : N×1 column vector of its timecourse

  % find the two ECG channels
  ecgLabels = data.label(contains(data.label,'ECG','IgnoreCase',true));
  if numel(ecgLabels) < 2
    error('Expected at least two ECG channels, found %d', numel(ecgLabels));
  end

  fs     = data.fsample;
  rawEEG = data.trial{1};  % [nCh × nSamples]

  % QRS band parameters
  bpFreq  = [5 40];   % Hz
  bpOrder = 4;

  snrScores = zeros(1,2);
  for i = 1:2
    % pull out channel
    idx = find(strcmp(data.label, ecgLabels{i}),1);
    sig = rawEEG(idx, :)';

    % band‑pass filter around QRS
    [b,a] = butter(bpOrder, bpFreq/(fs/2), 'bandpass');
    sigBP = filtfilt(b, a, sig);

    % noise = raw minus band‑passed
    noise = sig - sigBP;

    % compute SNR
    snrScores(i) = var(sigBP) ./ var(noise);
  end

  % pick the cleaner channel
  [~, best]      = max(snrScores);
  cleanECGchan   = ecgLabels{best};
  fprintf('→ picked ECG channel “%s” (SNR=%.1f vs %.1f)\n', ...
          cleanECGchan, snrScores(best), snrScores(3-best));

  % return the timecourse
  idx    = find(strcmp(data.label, cleanECGchan),1);
  ecgVec = data.trial{1}(idx, :)';  % column N×1
end
