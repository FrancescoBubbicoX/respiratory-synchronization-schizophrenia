%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                 %%%
%%%          BODY-BRAIN TIME SERIES PROCESSING      %%%
%%%                                                 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extracts heart, respiration, and EEG-derived time series from
% preprocessed data.
%
% Step 1 loads preprocessed signals and generates derived signals.
% Step 2 loads these derived signals, resamples them to a common 2-Hz time grid,
% and saves the final time series for network physiology analysis.

% Written by: Francesco Bubbico
% Last updated: June 2026

%% Step 1: Body-brain derived signals extraction

clear; clc;

addpath('C:\Users\francescob\Desktop\Open Data');
addpath('C:\Users\francescob\Desktop\Open Data\Code\Functions');

% Define directories
mainDir   = 'C:\Users\francescob\Desktop\Open Data\Data';
dataDir = fullfile(mainDir, 'Processed data', 'PreprocessedSignals');
outputDir = fullfile(mainDir, 'Processed data', 'BodyBrain_DerivedSignals');

validConds = {'lp','hp','pat','hea'};
selectedCond = lower(strtrim(input('Enter condition to process (lp, hp, pat, hea): ', 's')));
if ~ismember(selectedCond, validConds)
    error('Invalid condition selection.');
end
dataFolder = fullfile(dataDir, selectedCond);
outputFolder = fullfile(outputDir, selectedCond);

% Get list of input files
files = dir(fullfile(dataFolder, 'ReRef_*.mat'));
if isempty(files)
    fprintf('No ReRef files found in %s\n', dataFolder);  
    return;
end
fprintf('Found %d ReRef file(s) in %s\n', numel(files), dataFolder);

% Loop over files
for iFile = 1:numel(files)
    inPath = fullfile(files(iFile).folder, files(iFile).name);
    [~, baseName] = fileparts(inPath);
    subjectID = erase(baseName, 'ReRef_');
    bodyBrainFile      = fullfile(outputFolder, sprintf('bodyBrain_%s.mat', subjectID));
    if ~exist(bodyBrainFile, 'file')
        fprintf('Preparing body brain data for subject %s...\n', subjectID);
        S = load(inPath);                    % FieldTrip-type struct
        data = S.data;
        fs = data.fsample;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 1.1: RR interval extraction
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % -------------------------------------------------------------------------
        % Input:
        %   All available ECG recordings (EEG-derived ECG channels and BIOPAC ECG).
        %
        % ECG selection:
        %   The ECG channel providing the most reliable R-peak detection is
        %   automatically selected. Candidate channels are evaluated using:
        %   (1) the percentage of physiologically plausible RR intervals
        %       (0.30–2.00 s),
        %   (2) a QRS-band signal-to-noise ratio (5–40 Hz),
        %   (3) the proportion of missing data.
        %   When channels have comparable scores, the BIOPAC ECG is preferred.
        %
        % R-peak detection:
        %   R-peaks are detected using the Pan–Tompkins algorithm and refined to
        %   local ECG maxima. Detections occurring within 300 ms of a previous
        %   peak are discarded.
        %
        % Output:
        %   - R-peak times,
        %   - RR intervals,
        %   - RR validity mask.
        % -------------------------------------------------------------------------

        % --- ECG selection ---

        % Build candidates for best ECG
        cands = struct('label',{},'sig',{});
        if isfield(data,'ecg') && isfield(data.ecg,'label') && ~isempty(data.ecg.trial)
            rawECG = double(data.ecg.trial{1});  % [nECG x nSamples]
            for ii = 1:numel(data.ecg.label)
                cands(end+1).label = data.ecg.label{ii};
                cands(end).sig     = rawECG(ii,:); %#ok<SAGROW>
            end
        end
        if isfield(data,'ecg_biopac') && ~isempty(data.ecg_biopac)
            cands(end+1).label = 'ecg_biopac';
            cands(end).sig     = double(data.ecg_biopac(:))';
        end
        if isempty(cands), error('No ECG candidates found in %s', baseName); end

        % Detection-only filter for ECG quality scoring
        [b,a] = butter(4,[5 40]/(fs/2),'bandpass');
        bestIdx = 1; bestScore = -inf;

        minRRsec   = 0.30;
        minTrimSec = 10;                        % minimum 10s of finite data
        minTrimSamp = round(minTrimSec*fs);

        % Maximum score penalty for start/end missing data
        missingPenaltyW = 20;

        for ii = 1:numel(cands)
            sig = cands(ii).sig(:)';

            % Missingness on full signal
            finiteMask = isfinite(sig);
            validFrac  = nnz(finiteMask) / numel(sig);        % overall finite proportion (0..1)

            % Skip only if basically unusable
            if validFrac < 0.5
                cands(ii).snr_qrs     = -inf;
                cands(ii).pct_valid   = 0;
                cands(ii).score       = -inf;
                cands(ii).nR          = 0;
                cands(ii).validFrac   = validFrac;
                cands(ii).missingFrac = NaN;
                continue
            end

            % Trim start/end non-finite ONLY
            f1 = find(finiteMask, 1, 'first');
            f2 = find(finiteMask, 1, 'last');
            sigTrim = sig(f1:f2);

            trimFrac = numel(sigTrim) / numel(sig);           % proportion kept after trimming (0..1)
            missingFrac_edges = 1 - trimFrac;                 % penalize edge missingness

            % If signal is too short, force candidate to lose
            if numel(sigTrim) < minTrimSamp
                cands(ii).snr_qrs     = -inf;
                cands(ii).pct_valid   = 0;
                cands(ii).score       = -inf;
                cands(ii).nR          = 0;
                cands(ii).validFrac   = validFrac;
                cands(ii).missingFrac = missingFrac_edges;
                continue
            end

            % Working copy for detection only
            sigW = sigTrim - median(sigTrim,'omitnan');
            sigW = fillmissing(sigW, 'linear', 'EndValues', 'nearest');

            % QRS-band SNR on trimmed signal
            try
                sigBP = filtfilt(b,a,sigW);
                noise = sigW - sigBP;
                denom = max(var(noise,'omitnan'), eps);
                snr_qrs = var(sigBP,'omitnan') / denom;
            catch
                snr_qrs = 0;  % conservative fallback
            end

            % Pan–Tompkins detection on trimmed signal
            try
                [~, qrs_i_raw, ~] = pan_tompkin(sigW, fs, 0);
            catch
                try
                    [~, qrs_i_raw] = pan_tompkin(sigW, fs, 0);
                catch
                    qrs_i_raw = [];
                end
            end

            locs = qrs_i_raw(:)';
            locs = locs(locs>0 & locs<=numel(sigW));
            locs = unique(locs,'stable');

            % Remove detections occurring within the refractory period
            if ~isempty(locs)
                keep = [true, diff(locs) > round(minRRsec*fs)];
                locs = locs(keep);
            end

            % If too few peaks to define RR reliably, force candidate to lose
            if numel(locs) < 3
                cands(ii).snr_qrs     = -inf;
                cands(ii).pct_valid   = 0;
                cands(ii).score       = -inf;
                cands(ii).nR          = 0;
                cands(ii).validFrac   = validFrac;
                cands(ii).missingFrac = missingFrac_edges;
                continue
            end

            % RR validity (trimmed)
            RR_s = diff(locs)/fs;
            valid = (RR_s >= 0.30 & RR_s <= 2.00);
            pct_valid = 100 * nnz(valid) / max(1,numel(RR_s));

            % Composite ECG quality score
            baseScore = pct_valid + 5*log10(1+max(snr_qrs,0));

            % Penalize edge missingness
            penalty = missingPenaltyW * missingFrac_edges;

            score = baseScore - penalty;

            % Store diagnostics
            cands(ii).snr_qrs     = snr_qrs;
            cands(ii).pct_valid   = pct_valid;
            cands(ii).nR          = numel(locs);
            cands(ii).score       = score;
            cands(ii).validFrac   = validFrac;
            cands(ii).missingFrac = missingFrac_edges;

            % Prefer BIOPAC when scores are nearly identical
            if score > bestScore + 1 || (abs(score - bestScore) <= 1 && strcmp(cands(ii).label,'ecg_biopac'))
                bestScore = score; bestIdx = ii;
            end
        end

        ECG         = cands(bestIdx).sig(:)';     % chosen ECG
        selECGlabel = cands(bestIdx).label;

        % Print diagnostic
        diagStr = strings(1,numel(cands));
        for k = 1:numel(cands)
            mf = NaN; pv = NaN; sn = NaN;
            if isfield(cands(k),'missingFrac') && ~isempty(cands(k).missingFrac), mf = 100*cands(k).missingFrac; end
            if isfield(cands(k),'pct_valid')   && ~isempty(cands(k).pct_valid),   pv = cands(k).pct_valid; end
            if isfield(cands(k),'snr_qrs')     && ~isempty(cands(k).snr_qrs),     sn = cands(k).snr_qrs; end
            diagStr(k) = sprintf('%s: validRR=%.1f%%, SNR=%.2f, miss=%.1f%%', cands(k).label, pv, sn, mf);
        end

        fprintf('   ECG selected: %s | validRR=%.1f%% | QRS-SNR=%.2f | missing=%.1f%%\n  (%s)\n', ...
            selECGlabel, cands(bestIdx).pct_valid, cands(bestIdx).snr_qrs, 100*cands(bestIdx).missingFrac, ...
            strjoin(cellstr(diagStr), ' | '));


        % --- Final R-peak detection on selected ECG ---

        ecg = ECG(:)';  
        nanMask = ~isfinite(ecg);
        nNaN_ECG = nnz(nanMask);
        if nNaN_ECG > 0
            fprintf('   ECG contains %d non-finite samples (%.2f%%).\n\n', ...
                nNaN_ECG, 100*nNaN_ECG/numel(ecg));
        end

        % Trim start/end non-finite ONLY
        finiteMask = ~nanMask;
        if ~any(finiteMask)
            warning('ECG is all non-finite in %s. Skipping.', baseName);
            continue
        end
        f1 = find(finiteMask, 1, 'first');
        f2 = find(finiteMask, 1, 'last');
        ecgTrim = ecg(f1:f2);
        if numel(ecgTrim) < minTrimSamp
            warning('Too little finite ECG after trimming in %s (%.2f s). Skipping.', baseName, numel(ecgTrim)/fs);
            continue
        end

        % Copy ECG for R-peak detection
        ecgPT = ecgTrim - median(ecgTrim, 'omitnan');
        ecgPT = fillmissing(ecgPT, 'linear', 'EndValues', 'nearest');

        % Pan–Tompkins
        try
            [~, qrs_i_raw, ~] = pan_tompkin(ecgPT, fs, 0);
        catch
            try
                [~, qrs_i_raw] = pan_tompkin(ecgPT, fs, 0);
            catch
                warning('pan_tompkin failed in %s. Skipping.', baseName);
                continue
            end
        end
        locsR = qrs_i_raw(:)';
        locsR = locsR(locsR > 0 & locsR <= numel(ecgPT));
        locsR = unique(locsR,'stable');

        % Refine detections to local ECG maxima within ±50 ms
        winSec  = 0.05;                         
        winSamp = max(1, round(winSec * fs));
        locsRef = locsR;
        for k = 1:numel(locsR)
            i0 = locsR(k);
            idx = max(1, i0-winSamp) : min(numel(ecgTrim), i0+winSamp);
            [~, imax] = max(ecgTrim(idx));
            locsRef(k) = idx(imax);
        end
        locsR = unique(locsRef,'stable');

        % Remove detections within refractory period
        if ~isempty(locsR)
            keep = [true, diff(locsR) > round(minRRsec*fs)];
            locsR = locsR(keep);
        end

        % Shift back to original sample indices
        locsR = locsR + (f1 - 1);

        % Safety check: discard detections located on originally non-finite samples
        locsR = locsR(~nanMask(locsR));

        if numel(locsR) < 3
            warning('Few R-peaks detected in %s (n=%d). Skipping.', baseName, numel(locsR));
            continue
        end

        % RR intervals and validity mask
        RR_s    = diff(locsR) ./ fs;
        validRR = (RR_s >= 0.30) & (RR_s <= 2.00);

        % Attach output
        heart = struct();
        heart.rpeaks_time_s = locsR(:)' ./ fs;
        heart.rr_interval_s = RR_s(:)';
        heart.rr_valid      = validRR(:)';

        data.heart = heart;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 1.2: Respiratory phase extraction
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % -------------------------------------------------------------------------
        % Respiratory phase extraction
        %
        % Input:
        %   BIOPAC respiration trace.
        %
        % Processing:
        %   - Detect respiratory landmarks using trough→peak→trough segmentation.
        %   - Estimate continuous respiratory phase from the detected cycles.
        %   - Optionally mask portions of the signal not retained after quality checks.
        %
        % Main outputs:
        %   - resp: continuous respiration signal.
        %   - resp_masked: respiration signal after quality masking.
        %   - phase: continuous respiratory phase.
        %   - phase_masked: respiratory phase after quality masking.
        %
        % Reproducibility note:
        %   Default parameters are provided below. Some recordings required minor
        %   data-specific parameter adjustments after visual inspection of signal quality.
        % -------------------------------------------------------------------------

        plotDir = 'C:\Users\francescob\Desktop\Open Data\Plots\RespCycles';
        if ~isfolder(plotDir)
            mkdir(plotDir);
        end
        
        if isfield(data,'resp_biopac') && ~isempty(data.resp_biopac)
            RESP_row = double(data.resp_biopac(:)).';   % 1xL row
        else
            error('resp_biopac missing in %s.', baseName);
        end

        fprintf('   Respiratory cycle analysis ...\n');
        resp = RESP_row;
        L    = numel(resp);

        % Trim start/end non-finite ONLY
        finiteMask = isfinite(resp);
        if ~any(finiteMask)
            error('Respiration is all non-finite in %s.', baseName);
        end
        f1 = find(finiteMask, 1, 'first');
        f2 = find(finiteMask, 1, 'last');
        respTrim = resp(f1:f2);
        Ltrim    = numel(respTrim);
        if Ltrim < round(10*fs)
            warning('Too little finite respiration after trimming in %s (%.2f s).', baseName, Ltrim/fs);
        end

        % --- Optional start-amplitude correction parameters ---
        % Used only for a few selected recordings with clear initial
        % stabilization artifacts. The parameters were adjusted
        % after visual inspection of the affected initial segment.
        respCorrOpts = struct();
        % Detection and correction parameters
        respCorrOpts.start_s  = 4;      % initial window evaluated for abnormal amplitude
        respCorrOpts.thr_hi   = 3;      % correction threshold: start/reference amplitude ratio
        respCorrOpts.ramp_s   = 2;      % fade-out duration (s)
        respCorrOpts.cap_lo   = 0.33;   % lower bound for amplitude scaling
        respCorrOpts.cap_hi   = 1.66;   % upper bound for amplitude scaling
        respCorrOpts.make_plots = true;
        % Apply only when needed:
        % [respTrim, qc] = resp_start_envelope_correction(respTrim, fs, subjectID, respCorrOpts);


        % --- Respiratory phase estimation -----------------------------------------
        % Respiratory cycles are identified using trough–peak–trough segmentation.
        % Default parameters are used for most recordings. For noisier traces,
        % landmark-detection parameters, especially frac_iqr, may be adjusted after
        % visual inspection to improve cycle identification.

        params = struct( ...
            'interactive_split', false, ...  % set true for manual correction of long cycles
            'make_plots',        true, ...
            'browser_plot',      false, ...
            'frac_iqr',          0.3, ...   % prominence threshold = frac_iqr * IQR
            'min_dist_pk_s',     1, ...     % minimum distance between peaks (s)
            'min_dist_tr_s',     1, ...     % minimum distance between troughs (s)
            'min_width_s',       0.5 ...    % minimum width (s)
            );

        [phaseTrim, Ttrim, Ptrim, half_cycles] = resp_cycles_phase(respTrim, fs, params, plotDir, subjectID);

        % Map trimmed respiration and phase back to the original signal length
        resp = nan(1, L);
        resp(f1:f2) = respTrim;
        phase = nan(1, L);
        phase(f1:f2) = phaseTrim;
        T = Ttrim + (f1 - 1);
        P = Ptrim + (f1 - 1);
        half_cycles(:,1:2) = half_cycles(:,1:2) + (f1 - 1);

        % Define valid respiratory cycles from consecutive troughs containing one peak
        nTT = max(0, numel(T)-1);
        respiration_samp = nan(nTT,4);   % [inh_on, inh_off, exh_on, exh_off] in samples

        for k = 1:nTT
            i1 = T(k);
            i2 = T(k+1) - 1;                 
            pk = P(P > i1 & P < T(k+1));      
            if numel(pk) == 1
                respiration_samp(k,:) = [i1, pk-1, pk, i2];
            else
                respiration_samp(k,:) = [i1, NaN, NaN, i2];  
            end
        end

        % Initial validity mask: retain only complete trough–peak–trough cycles
        mask_keep = false(1,L);   
        for k = 1:size(respiration_samp,1)
            rowk = respiration_samp(k,:);
            if all(isfinite(rowk)) && rowk(1) < rowk(2) && rowk(2) < rowk(3) && rowk(3) < rowk(4)
                s = max(1, rowk(1)); e = min(L, rowk(4));
                mask_keep(s:e) = true;
            end
        end
        if ~isempty(T)
            mask_keep(1:max(1,T(1)-1)) = false;   
            mask_keep(min(L,T(end)):L)  = false;
        end
        nan_mask_base = ~mask_keep;               

        if f1 > 1, nan_mask_base(1:f1-1) = true; end
        if f2 < L, nan_mask_base(f2+1:end) = true; end

        % Convert cycle boundaries to seconds and retain valid cycles
        respiration = respiration_samp ./ fs; 
        valid_rows = all(isfinite(respiration),2) & ...
            respiration(:,1) < respiration(:,2) & ...
            respiration(:,2) < respiration(:,3) & ...
            respiration(:,3) < respiration(:,4);
        inhaleexhale = respiration(valid_rows,:);   % [inh_on, inh_off, exh_on, exh_off] in seconds

        % Identify cycle-duration outliers
        dur_s       = inhaleexhale(:,4) - inhaleexhale(:,1);
        is_dur_out  = isoutlier(dur_s,'quartiles','ThresholdFactor',5);
        dur_out_idx = find(is_dur_out);
        fprintf('Duration stats (s): min=%.3f  median=%.3f  mean=%.3f  max=%.3f\n', ...
            min(dur_s), median(dur_s), mean(dur_s), max(dur_s));
        fprintf('Duration outliers: %d\n', numel(dur_out_idx));
        if ~isempty(dur_out_idx)
            fprintf('Outlier cycle durations (s):\n');
            fprintf('  cycle %4d: %.3f s\n', [dur_out_idx, dur_s(dur_out_idx)].');
        end

        outlier_indices = [];
        if ~isempty(dur_out_idx)
            segSamp = round([inhaleexhale(dur_out_idx,1) inhaleexhale(dur_out_idx,4)] * fs);
            segSamp(:,1) = max(segSamp(:,1), 1);
            segSamp(:,2) = min(segSamp(:,2), L);
            dur_samp = arrayfun(@(s,e) s:e, segSamp(:,1), segSamp(:,2), 'UniformOutput', false);
            outlier_indices = unique([dur_samp{:}]);
        end

        % Remove duration-outlier cycles
        inhaleexhale(is_dur_out,:) = [];

        % Identify cycle-amplitude outliers
        nCycles    = size(inhaleexhale,1);
        amplitudes = nan(nCycles,1);
        for iC = 1:nCycles
            i_inh  = max(1, round(inhaleexhale(iC,1)*fs));
            i_peak = min(L, round(inhaleexhale(iC,3)*fs));
            if i_peak > i_inh
                seg = resp(i_inh:i_peak); seg = seg(isfinite(seg));
                if ~isempty(seg), amplitudes(iC) = (max(seg)-min(seg))/2; end
            end
        end
        is_amp_out  = isoutlier(amplitudes,'quartiles','ThresholdFactor',15);
        amp_out_idx = find(is_amp_out);
        fprintf('Amplitude outliers: %d\n\n', numel(amp_out_idx));

        if ~isempty(amp_out_idx)
            segSamp = round([inhaleexhale(amp_out_idx,1) inhaleexhale(amp_out_idx,4)] * fs);
            segSamp(:,1) = max(segSamp(:,1), 1);
            segSamp(:,2) = min(segSamp(:,2), L);
            amp_samp = arrayfun(@(s,e) s:e, segSamp(:,1), segSamp(:,2), 'UniformOutput', false);
            outlier_indices = unique([outlier_indices, [amp_samp{:}]]);
        end

        % Remove amplitude-outlier cycles
        inhaleexhale(is_amp_out,:) = [];

        % Final validity mask after cycle and outlier rejection
        mask_keep = false(1,L);
        if ~isempty(inhaleexhale)
            idx = round(inhaleexhale(:,[1 4]) * fs);
            idx(idx<1) = 1; idx(idx>L) = L;
            for i = 1:size(idx,1)
                mask_keep(idx(i,1):idx(i,2)) = true;
            end
        end
        nan_mask_final = ~mask_keep | nan_mask_base;   
        if ~isempty(outlier_indices)
            badS = unique(round(outlier_indices(:)'));
            badS = badS(isfinite(badS) & badS>=1 & badS<=L);
            nan_mask_final(badS) = true;
        end
        nan_mask_final = logical(nan_mask_final(:).'); 

        % Apply final mask to respiration and phase
        resp_row        = resp;                 
        resp_row_masked = nan(1,L);
        resp_row_masked(~nan_mask_final) = resp_row(~nan_mask_final);

        phase_cont   = phase;                       
        phase_masked = nan(size(phase_cont));
        phase_masked(~nan_mask_final) = phase_cont(~nan_mask_final);

        % Store continuous and masked respiration/phase outputs
        respiration_out = struct();
        respiration_out.resp              = resp_row;         % continuous respiration amplitude
        respiration_out.resp_masked       = resp_row_masked;  % masked respiration amplitude
        respiration_out.phase             = phase_cont;       % continuous respiration phase
        respiration_out.phase_masked      = phase_masked;     % masked respiration phase

        data.respiration = respiration_out;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 1.3: EEG band-power extraction
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % -------------------------------------------------------------------------
        % EEG band-power extraction
        %
        % Sliding-window spectral power is computed for selected EEG channels using
        % 2-s windows with 75% overlap. Power spectral density is estimated with
        % spectrogram and integrated within predefined frequency bands.
        %
        % Main output:
        %   - power_per_band: channel-level band-power time series.
        % -------------------------------------------------------------------------

        fs = data.fsample;

        % EEG channels for spectral-power extraction
        sel = {'Fp1','Fp2','F3','F4','C3','C4','P3','P4','O1','O2','F7','F8','T7','T8','P7','P8', ...
            'Fz','Cz','Pz','FC1','FC2','CP1','CP2','FC5','FC6','CP5','CP6','TP9','TP10', ...
            'F1','F2','C1','C2','P1','P2','AF3','AF4','FC3','FC4','CP3','CP4','PO3','PO4', ...
            'F5','F6','C5','C6','P5','P6','AF7','AF8','FT7','FT8','TP7','TP8','PO7','PO8', ...
            'Fpz','CPz','POz','Oz'};
        [tf, idx] = ismember(sel, data.label);
        assert(all(tf), 'Some selected electrodes not found in data.label.');
        EEGsel = double(data.trial{1}(idx, :));                 % [nCh x N]
        labels = sel(:);
        [nCh, ~] = size(EEGsel);

        % Frequency bands (Hz)
        band_defs = { ...
            'delta',  [0.5 3.5]; ...
            'theta',  [3.5 7.5]; ...
            'alpha',  [7.5 12.5]; ...
            'alpha1', [7.5 10]; ...
            'alpha2', [10 12.5]; ...
            'beta',   [12.5 30]; ...
            'beta1',  [12.5 20]; ...
            'beta2',  [20 30]; ...
            'gamma',  [30 45] };
        band_names = band_defs(:,1);
        B = numel(band_names);

        % Windowing: 2-s windows with 75% overlap, yielding one estimate every 0.5 s
        winSamp   = round(2*fs);
        noverlap  = round(0.75*winSamp);
        nfft      = max(1024, 2^nextpow2(winSamp));

        % Build frequency-bin indices for each band
        [~,F,T,~] = spectrogram(EEGsel(1,:), winSamp, noverlap, nfft, fs, 'psd');  % F: [nF x 1], T in s
        band_idx = cell(B,1);
        for b = 1:B
            fr = band_defs{b,2};
            if b < B
                band_idx{b} = find(F >= fr(1) & F < fr(2));
            else
                band_idx{b} = find(F >= fr(1) & F <= fr(2));
            end
        end

        % Allocate: per-band power time series [nCh x nT] for each band
        nT = numel(T);
        PEEG = struct();
        for b = 1:B
            PEEG.(band_names{b}) = zeros(nCh, nT);
        end

        % Estimate PSD and integrate power within each frequency band
        for ch = 1:nCh
            % Demean signal before PSD estimation
            x = EEGsel(ch,:) - mean(EEGsel(ch,:), 'omitnan');
            [~,~,~,Ppsd] = spectrogram(x, winSamp, noverlap, nfft, fs, 'psd'); 
            for b = 1:B
                idxb = band_idx{b};
                if ~isempty(idxb)
                    % Integrate PSD within the frequency band to obtain band power (µV²)
                    PEEG.(band_names{b})(ch,:) = trapz(F(idxb), Ppsd(idxb,:), 1);
                end
            end
        end

        % Store results
        data.np.peeg = struct( ...
            'labels',              {labels}, ...
            'bands',               {band_defs}, ...
            'win_s',               2, ...
            'overlap_pct',         75, ...
            'nfft',                nfft, ...
            'T_s',                 T, ...
            'power_per_band',      PEEG);          % struct of [nCh x nT] band-power matrices

        % Save derived heart, respiration, and EEG time-series outputs
        data_bodybrain = struct('heart',data.heart,'respiration',data.respiration,'np',data.np);
        save(bodyBrainFile, 'data_bodybrain');
    else
        fprintf('Body brain file exists. Skipping subject %s...\n', subjectID);
    end
end


%% Step 2: Build 2-Hz body-brain time series

clear; clc;

addpath('C:\Users\francescob\Desktop\Open Data');
addpath('C:\Users\francescob\Desktop\Open Data\Code\Functions');

% Define directories
mainDir   = 'C:\Users\francescob\Desktop\Open Data\Data';
dataDir = fullfile(mainDir, 'Processed data', 'BodyBrain_DerivedSignals');
outputDir = fullfile(mainDir, 'Processed data', 'BodyBrain_2HzTimeSeries');

validConds = {'lp','hp','pat','hea'};
selectedCond = lower(strtrim(input('Enter condition to process (lp, hp, pat, hea): ', 's')));
if ~ismember(selectedCond, validConds)
    error('Invalid condition selection.');
end
dataFolder = fullfile(dataDir, selectedCond);
outputFolder = fullfile(outputDir, selectedCond);

% Get list of input files
files = dir(fullfile(dataFolder, 'bodyBrain_*.mat'));
if isempty(files)
    fprintf('No bodyBrain files found in %s\n', dataFolder);
    return;
end
fprintf('Found %d bodyBrain file(s) in %s\n', numel(files), dataFolder);

% Loop over files
for iFile = 1:numel(files)
    inPath = fullfile(files(iFile).folder, files(iFile).name);
    [~, baseName] = fileparts(inPath);
    subjectID = erase(baseName, 'bodyBrain_');
    networkPhysiologyFile        = fullfile(outputFolder, sprintf('networkPhysiology_%s.mat', subjectID)); 
    if ~exist(networkPhysiologyFile, 'file')
        fprintf('Preparing Network Physiology data for subject %s...\n', subjectID);
        S = load(inPath);
        data = S.data_bodybrain;

        % Processing options
        RRexclude = 'interpolate';  % 'invalid' or 'interpolate'
        respMode = 'continuous';   % 'masked' or 'continuous'

        % --- Heart: RR intervals on the EEG-derived 2-Hz grid ---
        Tg = data.np.peeg.T_s(:)';        
        tR = data.heart.rpeaks_time_s(:)';    
        RR = data.heart.rr_interval_s(:)';    

        % Identify valid RR intervals
        if isfield(data.heart,'rr_valid') && numel(data.heart.rr_valid)==numel(RR)
            v = data.heart.rr_valid(:)';
            validRR = isfinite(v) & (v~=0);   % 1 = valid, 0 = invalid
        else
            validRR = true(1,numel(RR));
        end

        % Assign RR intervals to their midpoint times
        tRRmid = (tR(1:end-1) + tR(2:end))/2;   

        switch lower(RRexclude)
            case 'interpolate'
                % Interpolate RR intervals from valid RR midpoints.
                % Remaining gaps are filled only when a valid RR is available within maxGap.
                RRI_2Hz = interp1(tRRmid(validRR), RR(validRR), Tg, 'linear', NaN);
                maxGap = 1;  % seconds
                gapIdx = isnan(RRI_2Hz);

                if any(gapIdx)
                    Tgap = Tg(gapIdx);
                    RRI_fill = interp1(tRRmid(validRR), RR(validRR), Tgap, 'nearest', 'extrap');
                    % Keep nearest-neighbor fills only when close to observed data
                    nearestValid = interp1(tRRmid(validRR), tRRmid(validRR), Tgap, 'nearest', 'extrap');
                    dt = abs(Tgap - nearestValid);
                    keep = dt <= maxGap;
                    RRI_2Hz(gapIdx) = RRI_fill;                 
                    idxGap = find(gapIdx);                      
                    RRI_2Hz(idxGap(~keep)) = NaN;               
                    nFilled = nnz(keep);
                    fprintf('   RR: Filled %d values by nearest within %.1f s\n', nFilled, maxGap);
                end

            case 'invalid'
                % Assign RR intervals to the 2-Hz grid and remove samples that fall
                % within invalid RR intervals.
                RRI_2Hz = interp1(tRRmid(validRR), RR(validRR), Tg, 'previous', NaN);
                edges = [-Inf, (tR(2:end) + tR(1:end-1))/2, Inf];
                bins = discretize(Tg, edges); 
                badRR = ~validRR;                               
                idxBad = ~isnan(bins) & badRR(bins);
                RRI_2Hz(idxBad) = NaN;

            otherwise
                error('RRexclude must be ''invalid'' or ''interpolate''.');
        end

        % Store RR time series
        data.np.RRI_s = RRI_2Hz;

        % --- Respiration: amplitude and phase on the 2-Hz grid ---
        % Amplitude is averaged linearly; phase is averaged circularly.
        fs    = 500;
        Tg    = data.np.peeg.T_s(:)';      
        binW  = 0.5;                       
        halfW = binW/2;
        L     = numel(data.respiration.resp);
        tResp = (0:L-1)/fs;

        % Select continuous or masked respiration signals
        if strcmpi(respMode,'masked')
            ampSrc      = double(data.respiration.resp_masked(:))';
            phSrc       = double(data.respiration.phase_masked(:))';

            % Require at least 20% valid samples per bin
            minValid = ceil(0.20 * round(binW*fs));
        else
            ampSrc      = double(data.respiration.resp(:))';
            phSrc       = double(data.respiration.phase(:))';

            % Require at least one valid sample per bin
            minValid = 1;
        end

        amp_2Hz       = nan(size(Tg));     
        ph_2Hz        = nan(size(Tg));    

        for k = 1:numel(Tg)
            t1 = Tg(k) - halfW;
            t2 = Tg(k) + halfW;           
            idx = (tResp >= t1) & (tResp < t2);

            % Amplitude
            a = ampSrc(idx); a = a(isfinite(a));
            if numel(a) >= minValid
                amp_2Hz(k) = mean(a);
            end

            % Phase (circular mean)
            p = phSrc(idx); p = p(isfinite(p));
            if numel(p) >= minValid
                C = mean(cos(p)); S = mean(sin(p));
                ang = atan2(S, C);         
                if ang < 0, ang = ang + 2*pi; end
                ph_2Hz(k) = ang;           
            end
        end

        % Store respiration time series
        data.np.resp_amp          = amp_2Hz;         
        data.np.resp_phase        = ph_2Hz;      

        % Save final 2-Hz body-brain time series
        np = data.np;                                     
        save(networkPhysiologyFile, 'np');
    else
        fprintf('Network Physiology file exists. Skipping subject %s...\n', subjectID);
    end
end

