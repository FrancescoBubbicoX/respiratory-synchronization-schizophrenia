%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                 %%%
%%%                 PREPROCESSING DATA              %%%
%%%                                                 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Written by: Francesco Bubbico
% Last updated: June 2026

% This script processes the cut data files (e.g., "ID_cut.mat") using the following steps:
%
% 1: Filter, downsample, and align EEG/BIOPAC signals
% 2: Channel rejection and conservative ASR 
% 3: ICA decomposition and EOG/ECG component flagging
% 4: ICA components inspection & rejection
% 5: ASR
% 6: Channel interpolation and re-referencing


% Intermediate results are saved with descriptive filenames so that
% steps can be skipped if already completed.
%
% Adapt paths, channel selections, and function calls as necessary.

clear; clc;

%% Setup

% FieldTrip
addpath('C:\Program Files\MATLAB\R2024a\toolbox\fieldtrip-20250106');
ft_defaults;

% EEGLAB 
eeglabRoot = 'C:\Program Files\MATLAB\R2024a\toolbox\eeglab2022.0';
addpath(eeglabRoot); 
cd(eeglabRoot); 
eeglab;    % Adds EEGLAB subfolders to the MATLAB path

% Custom functions
addpath('C:\Users\francescob\Desktop\Open Data\Code\Functions');

% Define directories
mainDir = 'C:\Users\francescob\Desktop\Open Data\Data';
outputDir = fullfile(mainDir, 'Processed data', 'PreprocessedSignals');

% Select condition
validConds = {'lp','hp','pat','hea'};
selectedCond = lower(strtrim(input('Enter condition to process (lp, hp, pat, hea): ', 's')));
if ~ismember(selectedCond, validConds)
    error('Invalid condition selection.');
end

dataFolder = fullfile(mainDir, 'Source data', 'eeg', selectedCond);
outputFolder = fullfile(outputDir, selectedCond);

% Get list of input files
files = dir(fullfile(dataFolder, '*_cut.mat'));
if isempty(files)
    fprintf('No cut files found in %s\n', dataFolder);
    return;
end
fprintf('Found %d cut file(s) in %s\n', numel(files), dataFolder);

% Loop over subjects
for k = 1:numel(files)
    filePath = fullfile(dataFolder, files(k).name);
    [~, name, ~] = fileparts(files(k).name);
    subjectID = regexprep(name, '_cut$', '');

    fprintf('\n-------------------------\nProcessing subject: %s\n', subjectID);

    % Define filenames for each processing step:
    filterFile      = fullfile(outputFolder, sprintf('filter_%s.mat', subjectID));
    chrejFile        = fullfile(outputFolder, sprintf('chrej1_%s.mat', subjectID)); 
    compICAFile      = fullfile(outputFolder, sprintf('compICA_%s.mat', subjectID));
    rejectedCompFile = fullfile(outputFolder, sprintf('rejectedComp_%s.mat', subjectID));
    ASRFile = fullfile(outputFolder, sprintf('ASR2_%s.mat', subjectID)); 
    ReRefFile = fullfile(outputFolder, sprintf('ReRef_%s.mat', subjectID));

    % Skip entire preprocessing if the final output already exists
    if exist(ReRefFile, 'file')
        fprintf('Final preprocessed file already exists for subject %s. Skipping.\n', subjectID);
        continue;
    end

    %% 1: Filter, downsample, and align EEG/BIOPAC signals
    if ~exist(filterFile, 'file')
        fprintf('Filtering, downsampling, and EEG/BIOPAC alignment for subject %s...\n', subjectID);

        S = load(filePath);
        if ~isfield(S, 'data_cut')
            fprintf('File %s does not contain variable "data_cut". Skipping...\n', files(k).name);
            continue;
        end
        data = S.data_cut;

        % BIOPAC ECG is required for EEG/BIOPAC alignment
        if ~isfield(data, 'biopac') || ~isfield(data.biopac, 'ecg')
            msg = sprintf('\nSubject %s: BIOPAC ECG not found. Skipping...', subjectID);
            fprintf('%s\n', msg);
            continue;  % skip this subject
        end

        % --- Filter and downsample all EEG-system channels, including ECG and EOG ---
        if any(isnan(data.trial{1}(:)))
            warning('Subject %s: NaNs found in data_cut.trial. Skipping...', subjectID);
            continue;
        end

        cfg               = [];
        cfg.demean        = 'no';
        cfg.continuous    = 'yes';
        cfg.bpfilter      = 'yes';
        cfg.bpfreq        = [0.5 45];            % band-pass filter (0.5–45 Hz)
        cfg.bpfilttype    = 'firws';
        cfg.reref         = 'no';
        data_all          = ft_preprocessing(cfg, data);

        % Downsample to 500 Hz
        cfg               = [];
        cfg.resamplefs    = 500;
        cfg.detrend       = 'no';
        data_all          = ft_resampledata(cfg, data_all);

        % Keep header consistent
        data_all.hdr.Fs = data_all.fsample;

        % --- Filter and downsample ECG recorded with BIOPAC ---
        fs_target = 500;
        ecg_raw = data.biopac.ecg(:);
        fs_bio  = data.biopac.fs;

        % Temporarily remove edge NaNs before filtering
        nanMask = isnan(ecg_raw);
        if any(nanMask)
            first = find(~nanMask,1,'first');
            last  = find(~nanMask,1,'last');
            if isempty(first) || isempty(last)
                error('BIOPAC ECG is all NaN for subject %s.', subjectID);
            end
            lead_nan  = first-1;
            trail_nan = numel(ecg_raw)-last;
            ecg_seg   = ecg_raw(first:last);
        else
            lead_nan = 0; trail_nan = 0;
            ecg_seg  = ecg_raw;
        end

        % Convert the valid ECG segment to a FieldTrip structure
        ecg_bp = [];
        ecg_bp.label    = {'ECG_biopac'};
        ecg_bp.fsample  = fs_bio;
        ecg_bp.hdr.Fs   = fs_bio;
        ecg_bp.trial    = { ecg_seg' };
        ecg_bp.time     = { (0:numel(ecg_seg)-1)/fs_bio };

        % Apply the same band-pass filter used for the EEG data
        cfg               = [];
        cfg.demean        = 'no';
        cfg.continuous    = 'yes';
        cfg.bpfilter      = 'yes';
        cfg.bpfreq        = [0.5 45];
        cfg.bpfilttype    = 'firws';
        ecg_filt = ft_preprocessing(cfg, ecg_bp);

        % Downsample to 500 Hz
        cfg            = [];
        cfg.resamplefs = fs_target;
        cfg.detrend    = 'no';
        ecg_ds = ft_resampledata(cfg, ecg_filt);

        ecg_clean_seg = ecg_ds.trial{1}(:);

        % Reattach edge NaNs after downsampling
        lead_ds  = round(lead_nan  * fs_target / fs_bio);
        trail_ds = round(trail_nan * fs_target / fs_bio);

        ecg_clean_full = [nan(lead_ds,1); ecg_clean_seg; nan(trail_ds,1)];

        % Enforce the expected duration after resampling
        N_expected = round(numel(ecg_raw) * fs_target / fs_bio);
        if numel(ecg_clean_full) < N_expected
            ecg_clean_full(end+1:N_expected) = NaN;
        elseif numel(ecg_clean_full) > N_expected
            ecg_clean_full = ecg_clean_full(1:N_expected);
        end

        % Store in data structure
        data.biopac.ecg_clean    = ecg_clean_full;
        data.biopac.ecg_clean_fs = fs_target;

       % --- Preprocess respiration recorded with BIOPAC ---

        % BIOPAC respiration is required 
        if ~isfield(data, 'biopac') || ~isfield(data.biopac, 'resp')
            msg = sprintf('\nSubject %s: BIOPAC respiration not found. Skipping...\n', subjectID);
            fprintf('%s\n', msg);
            continue; 
        end

        % Temporarily remove edge NaNs before filtering
        raw_resp = data.biopac.resp(:);
        fs_bio   = data.biopac.fs;

        nanMask = isnan(raw_resp);
        if any(nanMask)
            first = find(~nanMask,1,'first');
            last  = find(~nanMask,1,'last');
            if isempty(first) || isempty(last)
                error('BIOPAC RESP is all NaN for subject %s.', subjectID);
            end
            lead_nan  = first - 1;
            trail_nan = numel(raw_resp) - last;
            resp_seg  = raw_resp(first:last);
        else
            lead_nan  = 0;
            trail_nan = 0;
            resp_seg  = raw_resp;
        end

        % Remove DC offset
        resp_seg = detrend(resp_seg,'constant');

        % Band-pass filter respiration signal
        [b,a] = butter(2,[0.05 2]/(fs_bio/2),'bandpass');
        resp_flt = filtfilt(b,a,resp_seg);

        % Mark extreme spikes as missing values
        thr = 4 * std(resp_flt,'omitnan');
        resp_nan = resp_flt;
        resp_nan(abs(resp_nan) > thr) = NaN;

        resp_clean = resp_nan;

        % Interpolate short gaps using a shape-preserving method
        maxGap = round(0.3 * fs_bio); % Up to 300 ms
        resp_clean = fillmissing(resp_clean,'pchip','MaxGap', maxGap);
        % Fallback for longer gaps
        if any(isnan(resp_clean))
            resp_clean = fillmissing(resp_clean,'makima');
        end
        % Only fill very short remaining gaps 
        if any(isnan(resp_clean))
            resp_clean = fillmissing(resp_clean,'pchip','MaxGap', round(0.05 * fs_bio));
        end
        % Final fallback 
        if any(isnan(resp_clean))
            resp_clean = fillmissing(resp_clean,'nearest');
        end

        % Smooth respiration signal
        resp_smooth = smoothdata(resp_clean,'movmean', round(0.05 * fs_bio));      
   
        % Downsample to 500 Hz
        rsp_bp = [];
        rsp_bp.label    = {'RESP_biopac'};
        rsp_bp.fsample  = fs_bio;
        rsp_bp.hdr.Fs   = fs_bio;
        rsp_bp.trial    = {resp_smooth'};
        rsp_bp.time     = {(0:numel(resp_smooth)-1)/fs_bio};

        cfg = [];
        cfg.resamplefs = fs_target;
        cfg.detrend    = 'no';
        rsp_ds = ft_resampledata(cfg, rsp_bp);

        resp_seg_ds = rsp_ds.trial{1}(:);

        % Reattach edge NaNs after downsampling
        lead_ds  = round(lead_nan  * fs_target / fs_bio);
        trail_ds = round(trail_nan * fs_target / fs_bio);
        resp_ds_full = [nan(lead_ds,1); resp_seg_ds; nan(trail_ds,1)];

        % Enforce the expected duration after resampling
        N_expected = round(numel(raw_resp) * fs_target / fs_bio);
        if numel(resp_ds_full) < N_expected
            resp_ds_full(end+1:N_expected) = NaN;
        elseif numel(resp_ds_full) > N_expected
            resp_ds_full = resp_ds_full(1:N_expected);
        end

        % Z-score valid samples only
        mu = mean(resp_ds_full,'omitnan');
        sd = std(resp_ds_full,[],'omitnan');

        resp_final = (resp_ds_full - mu) ./ sd;

        % Store in data structure
        data.biopac.resp_clean    = resp_final(:);
        data.biopac.resp_clean_fs = fs_target;
       

        % --- Align BIOPAC ECG and respiration to EEG recording ---
        data_biopac = data.biopac;

        % Select the cleanest ECG channel recorded with the EEG system
        [cleanECGchan, ecg_eeg] = pickCleanECG(data_all);
        ecg_eeg = ecg_eeg(:);          

        % Extract the preprocessed BIOPAC ECG signal
        ecg_bp  = data_biopac.ecg_clean(:);  
        fs = data_all.fsample; 

        % Compare EEG-derived ECG and BIOPAC ECG durations after downsampling.
        % Large mismatches are handled with RR-anchor rescue alignment.
        N_eeg   = numel(ecg_eeg);
        N_bp    = numel(ecg_bp);
        timeDiff = abs(N_eeg - N_bp)/fs;
        if timeDiff > 2
            fprintf('Subject %s: large mismatch (%.1f s) → RR-based rescue alignment\n', ...
                subjectID, timeDiff);

            % Detect R-peaks used as anchors for rescue alignment
            [locs_eeg, ~] = detect_R_peaks(ecg_eeg, fs);
            [locs_bp,  ~] = detect_R_peaks(ecg_bp,  fs);

            data_biopac = alignBiopacRRanchor( ...
                ecg_eeg, ...
                data_biopac.ecg_clean(:), ...
                data_biopac.resp_clean(:), ...
                locs_eeg, ...
                locs_bp, ...
                fs, ...
                cleanECGchan);
        else

            % Detect R-peaks in EEG-derived ECG and BIOPAC ECG
            [locs_eeg, peaks_eeg] = detect_R_peaks(ecg_eeg, fs);
            [locs_bp,  peaks_bp]  = detect_R_peaks(ecg_bp,  fs);

            fprintf('\nDetected %d EEG peaks, %d BIOPAC peaks\n', numel(locs_eeg), numel(locs_bp));

            % Estimate lag using cross-correlation between the two ECG signals
            maxLag = 2*fs;

            % Use the common valid duration for lag estimation
            L = min(numel(ecg_bp), numel(ecg_eeg));
            ecg_bp  = ecg_bp(1:L);
            ecg_eeg = ecg_eeg(1:L);

            % Use only overlapping, non-NaN samples
            valid = ~isnan(ecg_bp) & ~isnan(ecg_eeg);
            x = ecg_bp(valid);
            y = ecg_eeg(valid);

            % Robust normalization 
            x = (x - mean(x)) / std(x);
            y = (y - mean(y)) / std(y);

            % Cross-correlation
            [cc, lags] = xcorr(x, y, maxLag, 'coeff');
            [~, idx] = max(cc);
            bestLag = lags(idx);

            % Positive lag: "BIOPAC is late, shift BIOPAC backward."
            % Negative lag: "BIOPAC is early, shift BIOPAC forward."

            % Estimate lag from nearest-neighbour R-peak differences
            matchedLagSamples = nan(size(locs_eeg));
            for i = 1:numel(locs_eeg)
                [~, j] = min(abs(locs_bp - locs_eeg(i)));
                matchedLagSamples(i) = locs_bp(j) - locs_eeg(i);
            end
            
            % Retain plausible beat-to-beat lags within ±0.5 s
            maxSamples = round(0.5 * fs);
            valid = abs(matchedLagSamples) <= maxSamples;
            lagsPP = matchedLagSamples(valid);

            % Store the first peak lags for diagnostic comparison
            nPlot = min(3, numel(locs_eeg));
            firstLags = matchedLagSamples(1:nPlot);

            % Print alignment diagnostics
            fprintf('\n=== Alignment Summary (in samples) ===\n');
            fprintf('1) Cross‑correlation lag: %+d samples\n', bestLag);

            fprintf('\n2) Peak‑to‑peak: N=%d, mean=%+.1f, median=%+.1f\n', ...
                numel(lagsPP), mean(lagsPP), median(lagsPP));

            fprintf('\n3) First %d individual peak lags:\n', nPlot);
            for i = 1:nPlot
                fprintf('   Beat %d: %+d samples\n', i, firstLags(i));
            end
            fprintf('======================================\n\n');

            % Check whether cross-correlation and R-peak-based lag estimates agree
            allSmall = all(abs([bestLag, median(lagsPP), mean(firstLags)]) < 25);
            maxDiff  = max([bestLag, median(lagsPP), mean(firstLags)]) - ...
                min([bestLag, median(lagsPP), mean(firstLags)]);

            if allSmall && (maxDiff <= 10)
                % Automatically use cross-correlation when all lag estimates are consistent
                alignLag = bestLag;
                autoMode = true;
                fprintf('\nAuto-selected cross-correlation lag (%+d samples) since all metrics agreed (±10, <25 samples).\n\n', alignLag);
            else
                % If lag estimates disagree, manually select the most reliable estimate
                fprintf('Which lag should I use to align the BIOPAC signals?\n');
                fprintf('  [1] Cross‑correlation\n');
                fprintf('  [2] Peak‑to‑peak median\n');
                fprintf('  [3] First peak only\n');
                choice = input('Enter 1, 2, or 3: ');

                switch choice
                    case 1
                        alignLag = bestLag;
                        fprintf('→ Using cross‑correlation lag: %+d samples\n', alignLag);

                    case 2
                        alignLag = round(median(lagsPP));
                        fprintf('→ Using peak‑to‑peak median lag: %+d samples\n', alignLag);

                    case 3
                        alignLag = firstLags(1);
                        fprintf('→ Using first‑peak lag: %+d samples\n', alignLag);

                    otherwise
                        warning('Invalid choice—defaulting to cross‑correlation lag.');
                        alignLag = bestLag;
                        fprintf('→ Using cross‑correlation lag: %+d samples\n', alignLag);
                end
                autoMode = false;
            end

            N_eeg      = numel(ecg_eeg);
            clean_ecg  = data_biopac.ecg_clean(:);
            clean_resp = data_biopac.resp_clean(:);

            % Preallocate aligned BIOPAC signals using NaNs for non-overlapping edges
            ecg_bp_aligned  = nan(N_eeg,1);
            resp_bp_aligned = nan(N_eeg,1);

            if alignLag > 0
                % BIOPAC is late: shift it backward (advance by alignLag)
                d = alignLag;
                nCopy = min(N_eeg, numel(clean_ecg) - d);
                if nCopy > 0
                    ecg_bp_aligned(1:nCopy)  = clean_ecg((1:nCopy) + d);
                    resp_bp_aligned(1:nCopy) = clean_resp((1:nCopy) + d);
                end

            elseif alignLag < 0
                % BIOPAC is early: shift it forward (delay by |alignLag|)
                d = -alignLag;
                nCopy = min(N_eeg - d, numel(clean_ecg));
                if nCopy > 0
                    ecg_bp_aligned((1:nCopy) + d)  = clean_ecg(1:nCopy);
                    resp_bp_aligned((1:nCopy) + d) = clean_resp(1:nCopy);
                end

            else
                % No temporal shift required: truncate or pad to EEG length
                nCopy = min(N_eeg, numel(clean_ecg));
                ecg_bp_aligned(1:nCopy)  = clean_ecg(1:nCopy);
                resp_bp_aligned(1:nCopy) = clean_resp(1:nCopy);
            end

            % Diagnostic plot of EEG-derived ECG and aligned BIOPAC ECG
            % Build common time axis in seconds
            t = (0:N_eeg-1) / fs;
            % Z-score both ECG signals for visual comparison
            ecg_eeg_z    = (ecg_eeg    - mean(ecg_eeg))    / std(ecg_eeg);
            ecg_bp_z     = (ecg_bp_aligned - nanmean(ecg_bp_aligned)) / nanstd(ecg_bp_aligned);

            if autoMode
                % 15 s diagnostic plot
                tWin = t <= 15;
                figure;
                plot(t(tWin), ecg_eeg_z(tWin), 'b', 'LineWidth', 1.2); hold on;
                plot(t(tWin), ecg_bp_z(tWin),  'r', 'LineWidth', 1.2);
                xlabel('Time (s)'); ylabel('z-scored amplitude');
                title(sprintf('Subject %s – EEG ECG (%s, blue) vs BIOPAC ECG (red), First 15 s', ...
                    subjectID, cleanECGchan));
                legend({'EEG ECG','BIOPAC ECG'}); grid on;
            else
                % Interactive scroll plot
                eegAlign_scrollPlot(t, ecg_eeg_z, ecg_bp_z, cleanECGchan, 10);
            end

            % Store aligned BIOPAC signals
            data_biopac.ecg_aligned    = ecg_bp_aligned;
            data_biopac.resp_aligned   = resp_bp_aligned;
            data_biopac.aligned_fs     = fs;  
        end

        % --- Prepare filtered, downsampled, and aligned data for saving ---

        % Separate auxiliary ECG and EOG channels from EEG channels
        cfg               = [];
        cfg.channel       = ft_channelselection({'ECG*'}, data_all.label);
        data_ecg          = ft_selectdata(cfg, data_all);
        cfg.channel       = ft_channelselection({'IO'}, data_all.label);
        data_eog          = ft_selectdata(cfg, data_all);

        % Keep only EEG channels in the main FieldTrip structure
        cfg               = [];
        cfg.channel       = ft_channelselection({'all','-ECG*','-IO'}, data_all.label);
        data_eeg          = ft_selectdata(cfg, data_all);

        % Store auxiliary ECG and EOG channels separately
        data_eeg.ecg         = data_ecg;
        data_eeg.ocular      = data_eog;

        % Store aligned BIOPAC ECG and respiration signals
        data_eeg.ecg_biopac     = data_biopac.ecg_aligned(:);
        data_eeg.resp_biopac    = data_biopac.resp_aligned(:);
        data = data_eeg;

        % Ensure FieldTrip trial bookkeeping is present
        if ~isfield(data, 'sampleinfo')
            nsamp = size(data.trial{1}, 2);
            data.sampleinfo = [1 nsamp];
        end

        if ~isfield(data.cfg, 'trl')
            data.cfg.trl = [data.sampleinfo 0];
        end

        % Save filtered, downsampled, and aligned data
        save(filterFile, 'data');
        fprintf('Saved preprocessed file: %s\n', filterFile);

    else
        fprintf('Step 1 output already exists. Skipping filtering, downsampling, and alignment.\n');
    end

    %% 2: Channel rejection and conservative ASR 
        if ~exist(chrejFile, 'file')
            fprintf('Performing channel rejection and conservative ASR for subject %s...\n', subjectID);
            load(filterFile, 'data');

            % Identify and remove noisy channels using variance-based screening
            % followed by visual inspection. A high variance threshold is used
            % to flag only clearly suspicious channels.
            [data, removedChans] = ChannelRejectionManual(data, 5);

            % Store the removed channels for later interpolation
            data.removedChans = removedChans;

            % Update header information after channel removal.
            % This ensures compatibility with the EEGLAB-based ASR implementation.
            data.hdr.nChans = numel(data.label);
            data.hdr.label = data.label;
            data.hdr.Fs = data.fsample;
            if isfield(data.hdr, 'chantype')
                data.hdr.chantype = repmat({'EEG'}, data.hdr.nChans, 1);      
            end
            if isfield(data.hdr, 'chanunit')
                data.hdr.chanunit = repmat({'uV'}, numel(data.label), 1);
            end

            % Perform conservative Artifact Subspace Reconstruction (ASR).
            % The burst threshold was selected after visual comparison of ASR outputs
            % and configured to attenuate only high-amplitude transient artifacts.
            data = use_ASR(data, 120);

            save(chrejFile, 'data');
            fprintf('Saved channel-rejection/ASR file: %s\n', chrejFile);
        else
            fprintf('Step 2 output already exists. Skipping channel rejection and conservative ASR.\n');
        end

        %% 3: ICA decomposition and EOG/ECG component flagging
        if ~exist(compICAFile, 'file')
            fprintf('Running ICA for subject %s...\n', subjectID);
            load(chrejFile, 'data');

            % Estimate the number of ICA components from the data rank
            X = data.trial{1};
            [~, S, ~] = svd(X, 'econ');
            s = diag(S);
            cumVar = cumsum(s.^2) / sum(s.^2);
            nComp = find(cumVar >= 1, 1); % retain full data rank

            % Run ICA and flag components related to EOG/ECG activity
            compFlag = ICA_eog_ecg(data, subjectID, nComp);
            close all;

            save(compICAFile, 'compFlag');
            fprintf('Saved ICA components file: %s\n', compICAFile);
        else
            fprintf('Step 3 output already exists; skipping ICA.\n');
        end

        %% 4: ICA components inspection & rejection
        if ~exist(rejectedCompFile, 'file')
            fprintf('Inspecting ICA components for subject %s...\n', subjectID);
            load(chrejFile, 'data');

            if exist(compICAFile, 'file')
                load(compICAFile, 'compFlag');
            else
                error('ICA components file missing for subject %s.', subjectID);
            end

            % Store data before ICA component rejection for later comparison
            data_before_ica = data;

            % Inspect ICA components and select final components for rejection
            comp2removeFinal = InspectICAcomponents(compFlag);

            % Remove selected ICA components from the EEG data
            if ~isempty(comp2removeFinal)
                cfg = [];
                cfg.component = comp2removeFinal;
                data = ft_rejectcomponent(cfg, compFlag.comp, data);
                fprintf('Removed components: ');
                disp(comp2removeFinal);

                data.removedComp = comp2removeFinal;
            else
                fprintf('No components were removed.\n');
            end

            % Store data after ICA component rejection
            data_after_ica = data;

            % Quantify the percentage of total signal variance removed
            X_before = data_before_ica.trial{1};
            X_after  = data_after_ica.trial{1};
            var_before = sum(X_before(:).^2);
            var_after  = sum(X_after(:).^2);
            percentVarianceRemoved = 100 * (var_before - var_after) / var_before;
            fprintf('→ %.2f%% of the total signal variance was removed by rejecting ICA components.\n', percentVarianceRemoved);
            data.percentVarianceRemoved_ICA = percentVarianceRemoved;

            % Visual comparison of data before and after ICA component rejection
            comparePrePost_ICA(data_before_ica, data_after_ica);

            save(rejectedCompFile, 'data');
            fprintf('Saved ICA-rejected file: %s\n', rejectedCompFile);
        else
            fprintf('Step 4 output already exists. Skipping ICA component rejection.\n');
        end

    %% 5: ASR
        if ~exist(ASRFile, 'file')
            fprintf('Performing ASR for subject %s...\n', subjectID);
            load(rejectedCompFile, 'data');

            % Update header information before EEGLAB-based ASR
            data.hdr.nChans = numel(data.label);
            data.hdr.label = data.label;
            data.hdr.Fs = data.fsample;
            if isfield(data.hdr, 'chantype')
                data.hdr.chantype = repmat({'EEG'}, data.hdr.nChans, 1);     
            end
            if isfield(data.hdr, 'chanunit')
                data.hdr.chanunit = repmat({'uV'}, numel(data.label), 1);
            end

            % Perform a stricter ASR pass after ICA component rejection.
            % This step targets residual artifacts remaining after ICA.
            data = use_ASR(data, 15);

            save(ASRFile, 'data');
                fprintf('Saved ASR file: %s\n', ASRFile);
        else
            fprintf('Step 5 output already exists. Skipping ASR.\n');
        end

        %% 6: Channel interpolation and re‑referencing
        if ~exist(ReRefFile, 'file')
            fprintf('Performing channel interpolation and re-referencing for subject %s...\n', subjectID);
            load(ASRFile, 'data');

            if ~isempty(data.removedChans)
                fprintf('Interpolating %d removed channel(s): %s\n', ...
                    numel(data.removedChans), strjoin(data.removedChans, ', '));

                % Define the complete EEG channel montage before channel removal
                full_labels = {
                    'Fp1', 'Fp2', 'F3', 'F4', 'C3', 'C4', 'P3', 'P4', 'O1', 'O2', ...
                    'F7', 'F8', 'T7', 'T8', 'P7', 'P8', 'Fz', 'Cz', 'Pz', ...
                    'FC1', 'FC2', 'CP1', 'CP2', 'FC5', 'FC6', 'CP5', 'CP6', ...
                    'TP9', 'TP10', 'F1', 'F2', 'C1', 'C2', 'P1', 'P2', ...
                    'AF3', 'AF4', 'FC3', 'FC4', 'CP3', 'CP4', 'PO3', 'PO4', ...
                    'F5', 'F6', 'C5', 'C6', 'P5', 'P6', ...
                    'AF7', 'AF8', 'FT7', 'FT8', 'TP7', 'TP8', 'PO7', 'PO8', ...
                    'Fpz', 'CPz', 'POz', 'Oz'
                    };

                % Build dummy structure containing the full channel montage for neighbour preparation
                data_dummy = [];
                data_dummy.label = full_labels(:); 
                data_dummy.fsample = 1000;
                data_dummy.trial = {nan(length(full_labels), 100)};  
                data_dummy.time = {linspace(0, 0.1, 100)};            
                data_dummy.hdr = [];

                % Prepare channel neighbours from the template layout
                cfg_neighb = [];
                cfg_neighb.method = 'template'; 
                cfg_neighb.layout = 'easycapM11.mat';
                neighbours = ft_prepare_neighbours(cfg_neighb, data_dummy);

                % Load electrode coordinates for spline interpolation
                electrodeFile = 'C:\Users\francescob\Desktop\Open Data\Code\Functions\Standard-10-20-Cap62+ref_EDIT.txt';
                assert(exist(electrodeFile,'file')==2, 'Electrode file not found: %s', electrodeFile);
                Ttxt = readtable(electrodeFile);
                elec.label = Ttxt.labels;
                elec.chanpos = [Ttxt.X, Ttxt.Y, Ttxt.Z];
                elec.elecpos = elec.chanpos;
                elec.unit = 'm';
                elec.coordsys= 'EEG';

                % Interpolate removed EEG channels 
                cfg = [];
                cfg.method     = 'spline';                    
                cfg.badchannel = data.removedChans;
                cfg.neighbours = neighbours;
                cfg.elec       = elec;            
                data = ft_channelrepair(cfg, data);

                % Restore the original channel order after interpolation
                [found, idx] = ismember(full_labels, data.label);
                if ~all(found)
                    missing = full_labels(~found);
                    error('These channels are missing from data.label:\n %s', strjoin(missing, ', '));
                end
                data.label = data.label(idx);
                data.trial{1} = data.trial{1}(idx, :);
                data.hdr.label = data.label;
                data.hdr.nChans = numel(data.hdr.label);
                if isfield(data.hdr, 'chantype')
                    data.hdr.chantype = repmat({'EEG'}, data.hdr.nChans, 1);
                end
                if isfield(data.hdr, 'chanunit')
                    data.hdr.chanunit = repmat({'uV'}, data.hdr.nChans, 1);
                end
                data.sampleinfo = [1, size(data.trial{1},2)];
                data.hdr.nSamples = data.sampleinfo(2);
                data.hdr.nTrials = 1;
            else
                fprintf('No channels to interpolate.\n');
            end

            % Define the list of non-standard fields to preserve
            extraFields = {'ecg', 'ecg_biopac', 'resp_biopac', 'removedChans'};
            extraData = struct();
            for i = 1:numel(extraFields)
                if isfield(data, extraFields{i})
                    extraData.(extraFields{i}) = data.(extraFields{i});
                end
            end

            % Apply common-average re-reference over all EEG channels
            cfg = [];
            cfg.reref       = 'yes';
            cfg.refchannel  = 'all';
            cfg.rerefmethod = 'avg';
            data = ft_preprocessing(cfg, data);
            fprintf('Applied common‑average reference over %d channels.\n', numel(data.label));

            % Reinsert the saved fields into the rereferenced data
            extraFieldNames = fieldnames(extraData);
            for i = 1:numel(extraFieldNames)
                data.(extraFieldNames{i}) = extraData.(extraFieldNames{i});
            end

            if isfield(data, 'cfg')
                data = rmfield(data, 'cfg');
            end

            % Save the re‑referenced data
            save(ReRefFile, 'data');
            fprintf('Saved re‑referenced file: %s\n', ReRefFile);
        else
            fprintf('Step 6 output already exists. Skipping channel interpolation and re-referencing.\n');
        end

end
