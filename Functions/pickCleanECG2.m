function [cleanECGchan, ecgVec, diagTable] = pickCleanECG2(data)
% pickCleanECG   Pick the cleanest ECG channel via Pan–Tompkins QC + SNR
%
% [cleanECGchan, ecgVec, diagTable] = pickCleanECG(data)
%   data          : ft-style struct with optional .ecg (struct) and/or .ecg_biopac
%   cleanECGchan  : name of the chosen channel (string)
%   ecgVec        : N×1 column vector (with NaNs preserved at start/end)
%   diagTable     : table with diagnostics for all candidates
%
% Requires pan_tompkin.m on the path.

fs = data.fsample;

% --- build candidate list directly from data.label
cands = struct('label',{},'sig',{});

ecgIdx = find(contains(data.label,'ECG','IgnoreCase',true));
if isempty(ecgIdx)
    error('No ECG channels found in this dataset.');
end

raw = double(data.trial{1});  % [nCh x nSamples]
for ii = 1:numel(ecgIdx)
    cands(ii).label = data.label{ecgIdx(ii)};
    cands(ii).sig   = raw(ecgIdx(ii),:);
end

% --- QRS band filter
[b,a] = butter(4,[5 40]/(fs/2),'bandpass');

bestIdx   = 1;
bestScore = -inf;
diagRows  = [];

for ii = 1:numel(cands)
    sig = cands(ii).sig(:)';

    % remove leading/trailing NaNs
    finiteMask = isfinite(sig);
    if nnz(finiteMask) < 0.5*numel(sig)
        snr_qrs   = -inf;
        pct_valid = 0;
        nR        = 0;
        score     = -inf;
    else
        firstValid = find(finiteMask,1,'first');
        lastValid  = find(finiteMask,1,'last');
        sigValid   = sig(firstValid:lastValid);

        % median-detrend and fill small gaps inside
        sigW = fillmissing(sigValid - median(sigValid,'omitnan'), ...
            'linear','EndValues','nearest');

        % SNR
        sigBP = filtfilt(b,a,sigW);
        noise = sigW - sigBP;
        snr_qrs = var(sigBP,'omitnan') / max(var(noise,'omitnan'), eps);

        % Pan–Tompkins detection
        [~, qrs_i_raw, delay] = pan_tompkin(sigW, fs, 0);
        locs = qrs_i_raw(:)';
        if ~isempty(delay) && isnumeric(delay)
            locs = locs - round(delay);
        end
        locs = locs(locs>0 & locs<=numel(sigW));
        if ~isempty(locs)
            keep = [true, diff(locs) > round(0.30*fs)];
            locs = locs(keep);
        end

        RR_s = diff(locs)/fs;
        valid = (RR_s >= 0.30 & RR_s <= 2.00);
        pct_valid = 100 * nnz(valid) / max(1,numel(RR_s));

        % Composite score
        score = pct_valid + 5*log10(1+max(snr_qrs,0));

        nR = numel(locs);
    end

    cands(ii).snr_qrs   = snr_qrs;
    cands(ii).pct_valid = pct_valid;
    cands(ii).nR        = nR;
    cands(ii).score     = score;

    diagRows = [diagRows; {cands(ii).label, snr_qrs, pct_valid, nR, score}]; %#ok<AGROW>
end

% --- find the best channel
allScores = [cands.score];
[bestScore, bestIdx] = max(allScores);

cleanECGchan = cands(bestIdx).label;
ecgVec       = cands(bestIdx).sig(:); % return original signal (NaNs preserved)

fprintf('→ Picked ECG channel "%s" (score=%.2f)\n', ...
    cleanECGchan, bestScore);

% --- diagnostics table
diagTable = cell2table(diagRows, ...
    'VariableNames', {'Channel','SNR_QRS','PctValidRR','nR','Score'});
end
