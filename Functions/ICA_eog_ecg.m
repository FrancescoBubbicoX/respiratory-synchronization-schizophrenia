function compFlag = ICA_eog_ecg(data, subjectID, nComp)
%   Run ICA on EEG channels and flag components associated with ocular
%   and cardiac activity.
%
%   Components are flagged using:
%   1. cross-correlation with the EOG channel;
%   2. cross-correlation with the cleanest ECG channel;
%   3. frequency-domain coherence with ECG in the cardiac frequency range.
%
% INPUTS:
%   data      - FieldTrip structure containing EEG data, with optional
%               auxiliary fields data.ocular, data.ecg, and data.ecg_biopac.
%   subjectID - Subject identifier used for printed messages/output.
%   nComp     - Number of ICA components to compute. Default = 15.
%
% OUTPUT:
%   compFlag  - Structure containing the ICA decomposition, EOG/ECG
%               correlations, ECG coherence, and flagged component indices.

%% --- Settings ---
if nargin < 3
    nComp = 15;   % default number of ICA components
end
thresholdOcular = 0.2;  % correlation threshold for ocular
thresholdECG    = 0.2;  % correlation threshold for ECG (time-domain)
cohThresholdECG = 0.2;   % coherence threshold for ECG (frequency)

%% Prepare BIOPAC ECG as FieldTrip structure, if available
if isfield(data,'ecg_biopac') && ~isempty(data.ecg_biopac)
    ecg_bp = [];
    ecg_bp.label    = {'ECG_biopac'};
    ecg_bp.fsample  = data.fsample;                  % same sampling as EEG
    ecg_bp.trial    = {data.ecg_biopac(:)'};         % row vector
    ecg_bp.time     = {data.time{1}};                % reuse time axis
else
    ecg_bp = [];
end

%% Select auxiliary channels used for component flagging
% Ocular
if isfield(data, 'ocular') && ~isempty(data.ocular)
    dOcular = data.ocular;
else
    fprintf('Subject %s: No ocular channel available.\n', subjectID);
    dOcular = [];
end

% ECG
if isfield(data,'ecg') && ~isempty(data.ecg)
    dECG = data.ecg;
else
    dECG = [];
end
% If we also have biopac ECG, append it
if ~isempty(ecg_bp)
    if isempty(dECG)
        dECG = ecg_bp;
    else
        dECG = ft_appenddata([], dECG, ecg_bp);
    end
end


%% Run ICA on EEG channels
cfg = [];
cfg.channel = data.label;  % use all EEG channels
cfg.continuous = 'yes';
cfg.method = 'fastica';
cfg.numcomponent = nComp;
comp = ft_componentanalysis(cfg, data);

%% Extract ICA component time courses
cfg = [];
cfg.keeptrials = 'yes';
tlComp = ft_timelockanalysis(cfg, comp);
if iscell(tlComp.time)
    t_comp = tlComp.time{1};
else
    t_comp = tlComp.time;
end
if iscell(tlComp.trial)
    compMatrix = tlComp.trial{1};
else
    compMatrix = squeeze(tlComp.trial(1,:,:)); % [nComponents x nTimePoints]
end
numComponents = size(compMatrix, 1);


%% EOG x ICA Components cross-correlation
if ~isempty(dOcular)
    % Timelock to get continuous trial & time vector
    cfg = [];
    cfg.keeptrials = 'yes';
    tlOcular = ft_timelockanalysis(cfg, dOcular);
    if iscell(tlOcular.trial)
        ocularMatrix = tlOcular.trial{1};
    else
        ocularMatrix = tlOcular.trial;
    end
    xOcular = ocularMatrix(1, :)';            % take first EOG channel
    if iscell(tlOcular.time)
        t_oco = tlOcular.time{1};
    else
        t_oco = tlOcular.time;
    end
    % Interpolate EOG onto ICA time axis if needed
    if length(t_oco) ~= length(t_comp)
        fprintf('Interpolating ocular signal...\n');
        xOcular = interp1(t_oco, xOcular, t_comp, 'linear', 'extrap');
    end

    % Zero-mean the EOG once
    xOcular = xOcular(:) - mean(xOcular);

    % Cross-correlate each ICA component with EOG
    maxLagSec  = 0.1;                       % ±100 ms window
    fs         = dOcular.fsample;
    maxLagSamp = round(maxLagSec * fs);

    % Prepare outputs
    rOcular   = zeros(numComponents,1);
    lagOcular = zeros(numComponents,1);

    for c = 1:numComponents
        y = compMatrix(c, :)';
        y = y - mean(y);

        % normalized cross-correlation
        [xc, lags] = xcorr(y, xOcular, maxLagSamp, 'coeff');

        % find peak absolute correlation
        [rOcular(c), idx] = max(abs(xc));
        lagOcular(c)     = lags(idx) / fs;  % in seconds
    end
end

%% ECG x ICA Components cross-correlation
if ~isempty(dECG)
    % 1) Pick the cleanest ECG channel automatically
    [cleanECGchan, ecg_eeg] = pickCleanECG2(dECG);

    % Prepare a new FT struct with only the cleanest ECG channel
    dECG_best = dECG;
    dECG_best.trial{1} = ecg_eeg(:)';  % row vector
    dECG_best.label    = {cleanECGchan};
    ecgChannel         = cleanECGchan;

    %--- 2) Extract ECG data & time
    cfg_tl = [];
    cfg_tl.keeptrials = 'yes';
    tlECG = ft_timelockanalysis(cfg_tl, dECG_best);

    if iscell(tlECG.trial)
        ecgMatrix = tlECG.trial{1};
    else
        ecgMatrix = tlECG.trial;
    end
    xECG = ecgMatrix(1,:);

    if iscell(tlECG.time)
        t_ecg = tlECG.time{1};
    else
        t_ecg = tlECG.time;
    end

    % 3) Remove leading/trailing NaNs (no interpolation)
    finiteMask = isfinite(xECG);
    if any(~finiteMask)
        firstValid = find(finiteMask,1,'first');
        lastValid  = find(finiteMask,1,'last');
        xECG       = xECG(firstValid:lastValid);
        t_ecg      = t_ecg(firstValid:lastValid);
    end
    xECG = xECG(:);

    % 4) Restrict ICA components to the overlapping time window
    validIdx = (t_comp >= t_ecg(1)) & (t_comp <= t_ecg(end));
    compMatrix_valid = compMatrix(:, validIdx);
    t_comp_valid     = t_comp(validIdx);

    % Match lengths if needed (truncate to shortest)
    minLen = min(length(t_comp_valid), length(xECG));
    compMatrix_valid = compMatrix_valid(:,1:minLen);
    xECG             = xECG(1:minLen);

    %--- 5) Cross-correlation
    maxLagSec  = 0.1;                     % ±100 ms window
    fs         = dECG.fsample;
    maxLagSamp = round(maxLagSec * fs);

    % Zero-mean ECG once
    xECG = xECG - mean(xECG);

    % Preallocate
    rECG   = zeros(numComponents,1);
    lagECG = zeros(numComponents,1);

    for c = 1:numComponents
        y = compMatrix_valid(c,:)';
        y = y - mean(y);

        [xc,lags] = xcorr(y, xECG, maxLagSamp, 'coeff');
        [rECG(c), idx] = max(abs(xc));
        lagECG(c)      = lags(idx) / fs;
    end
end


%% Flag components based on ocular correlation
comp2remove_eog = find(abs(rOcular) >= thresholdOcular);

%% Flag components based on ECG correlation
comp2remove_ecg_corr = find(abs(rECG) >= thresholdECG);

%% Additional: ECG artifact detection using frequency coherence
compWithECG = ft_appenddata([], comp, dECG_best);

cfg = [];
cfg.method = 'mtmfft';
cfg.output = 'fourier';
cfg.foilim = [0 100];
cfg.taper = 'hanning';
cfg.pad = 'maxperlen';
freq = ft_freqanalysis(cfg, compWithECG);

cfg_trl = [];
cfg_trl.length = 5;
data_epoched = ft_redefinetrial(cfg_trl, compWithECG);

cfg = [];
cfg.method = 'mtmfft';
cfg.output = 'fourier';
cfg.foilim = [0 100];
cfg.taper = 'dpss';
cfg.tapsmofrq = 2;
cfg.pad = 'maxperlen';
freq = ft_freqanalysis(cfg, data_epoched);

icaLabels = freq.label(1:numComponents);
ecgLabel  = freq.label{end};

cfg = [];
cfg.channelcmb = {icaLabels, {ecgLabel}};
cfg.jackknife = 'no';
cfg.method    = 'coh';
fdcomp = ft_connectivityanalysis(cfg, freq);

ecgBand = [0.5 3];
bandIdx = find(fdcomp.freq >= ecgBand(1) & fdcomp.freq <= ecgBand(2));
meanCoh = mean(abs(fdcomp.cohspctrm(:, bandIdx)), 2);

% Plot the coherence spectrum.
% figure;
% subplot(2,1,1);
% plot(fdcomp.freq, abs(fdcomp.cohspctrm));
% xlim([0.5 3]);
% xlabel('Frequency (Hz)');
% ylabel('Coherence');
% subplot(2,1,2);
% imagesc(fdcomp.freq, 1:size(fdcomp.cohspctrm,1), abs(fdcomp.cohspctrm));
% axis xy;
% xlim([0.5 3]);
% xlabel('Frequency (Hz)');
% ylabel('ICA Component');
% colorbar;
% title('Coherence Spectrum');
% drawnow;

comp2remove_ecg_coh = find(meanCoh >= cohThresholdECG);

%% Prepare output structure.
compFlag.subjectID = subjectID;
compFlag.dOcular = dOcular;
compFlag.dECG = dECG_best; 
compFlag.comp = comp;               
compFlag.rOcular = rOcular;
compFlag.rECG = rECG;
compFlag.cohECG = meanCoh;
compFlag.comp2remove_eog = comp2remove_eog;
compFlag.comp2remove_ecg = comp2remove_ecg_corr;
compFlag.comp2remove_ecgCoh = comp2remove_ecg_coh;
end

