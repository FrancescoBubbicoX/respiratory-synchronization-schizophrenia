function comp2removeFinal = InspectICAcomponents(compFlag)
%   Inspect ICA components and manually select components for rejection.
%
%   The function summarizes components associated with EOG and ECG activity,
%   displays spectral and time-domain screening metrics, opens the FieldTrip
%   component browser, and asks the user to select and confirm the final
%   components to remove.
%
%   Component rejection is not automatic. The final decision is made after
%   visual inspection of component time courses, topographies, and their
%   relationship with EOG/ECG signals.
%
% INPUT:
%   compFlag - Structure produced by ICA_eog_ecg containing:
%                .comp     : ICA component structure
%                .rOcular  : EOG-component correlation values
%                .rECG     : ECG-component correlation values
%                .cohECG   : ECG-component coherence values
%
% OUTPUT:
%   comp2removeFinal - Vector of ICA component indices selected for removal.

%% 1. Print Top 3 Components by Correlation/Coherence
% For ocular (EOG)
[~, sortedOcularIdx] = sort(abs(compFlag.rOcular), 'descend');
nOcular = min(3, numel(sortedOcularIdx));
top3Ocular = sortedOcularIdx(1:nOcular);

% For ECG
[~, sortedECGIdx] = sort(abs(compFlag.rECG), 'descend');
nECG = min(3, numel(sortedECGIdx));
top3ECG = sortedECGIdx(1:nECG);

% For ECG coherence
[~, sortedECGcohIdx] = sort(abs(compFlag.cohECG), 'descend');
nECG = min(3, numel(sortedECGcohIdx));
top3ECGcoh = sortedECGcohIdx(1:nECG);

%% 2. Compute spectral/time-domain screening metrics
% Compute spectral and time-domain metrics for ICA components
cfg = [];
cfg.method = 'mtmfft'; cfg.output = 'pow';
cfg.taper  = 'hanning'; cfg.foilim = [0 45]; cfg.pad = 'maxperlen';
freqComp = ft_freqanalysis(cfg, compFlag.comp);

nComp      = numel(compFlag.comp.label);
hfLfRatio  = zeros(nComp,1);
kurtVals   = zeros(nComp,1);
for ic = 1:nComp
    psd = squeeze(freqComp.powspctrm(ic,:)); f = freqComp.freq;
    lf = trapz(f(f>=1 & f<=20),    psd(f>=1 & f<=20));
    hf = trapz(f(f>20 & f<=45),     psd(f>20 & f<=45));
    hfLfRatio(ic) = hf / lf;
    tc = compFlag.comp.trial{1}(ic,:);
    kurtVals(ic) = kurtosis(tc);
end

% -----------------------
% Screening thresholds:
% identify components with relatively high HF/LF ratio and/or kurtosis.
% These flags are used only to guide visual inspection, not for automatic rejection.
pHigh = 50; % percentile
hfLfThresh = prctile(hfLfRatio, pHigh);
kurtThresh = prctile(kurtVals,   pHigh);

flagHF = hfLfRatio >= hfLfThresh;
flagK  = kurtVals   >= kurtThresh;
% Require both HF/LF and kurtosis for strict artifact flagging
flags  = find(flagHF & flagK);
% If none meet both criteria, highlight components high on either metric
if isempty(flags)
    warning('No components met both HF/LF and kurtosis thresholds: relaxing to either criterion.');
    flags = find(flagHF | flagK);
end

% Plot PSD of Flagged Components Only
nF  = numel(flags);
nCol = 4;
nRow = ceil(nF/nCol);
figure('Name','Flagged ICA PSDs','Color','w','WindowState','maximized');
for i = 1:nF
    ic = flags(i);
    subplot(nRow,nCol,i);
    plot(freqComp.freq, squeeze(freqComp.powspctrm(ic,:)));
    xlim([0 45]);
    title(sprintf('Comp %d: HF/LF=%.2f, Kurt=%.1f', ic, hfLfRatio(ic), kurtVals(ic)));
    grid on;
    if i > (nRow-1)*nCol, xlabel('Hz'); end
    if mod(i-1,nCol)==0, ylabel('Power'); end
end
sgtitle('Components flagged by HF/LF or Kurtosis');
fprintf('Inspect flagged PSDs. Close figure when ready.\n');

%% 3. Visual Inspection via Data Browser
% Print the flagged components along with the top 3 for EOG and ECG.
fprintf('\n\n---------- Flagged Components ----------\n');
fprintf('EOG: ');
for i = 1:length(top3Ocular)
    fprintf('%d (%.2f) ', top3Ocular(i), compFlag.rOcular(top3Ocular(i)));
end
fprintf('\n');
fprintf('ECG: ');
for i = 1:length(top3ECG)
    fprintf('%d (%.2f) ', top3ECG(i), compFlag.rECG(top3ECG(i)));
end
fprintf('\n');
fprintf('ECG coherence: ');
for i = 1:length(top3ECGcoh)
    fprintf('%d (%.2f) ', top3ECGcoh(i), compFlag.cohECG(top3ECGcoh(i)));
end
fprintf('\n');
fprintf('----------------------------------------\n');

% Launch the data browser for further inspection.
cfg = [];
cfg.layout = 'easycapM11.mat'; 
cfg.viewmode = 'component';
cfg.blocksize = 30;
cfg.ylim = [-250 250];
% Suppress ft_databrowser's internal instructions.
dummy = evalc('ft_databrowser(cfg, compFlag.comp);');
set(gcf, 'Color','white', 'WindowState','maximized');

fprintf('\nInspect the components in the data browser.\n');
fprintf('Note the components you consider problematic.\n');
str = input('Enter components to remove (e.g. 1 3 5):\n','s');
tmp = str2num(str);  %#ok<ST2NM> converts '1 3 5' → [1 3 5]
tmp = sort(unique(tmp));

close(gcf);
close(gcf);

%% 4. Topographical and Time Series Inspection and Final Confirmation
if ~isempty(tmp)
    confirmed = false;
    while ~confirmed
        % Show topoplots
        figureHandle = figure('Position', [300, 300, 600, 400]);
        cfg = [];
        cfg.component = tmp;
        cfg.layout = 'easycapM11.mat';
        evalc('ft_topoplotIC(cfg, compFlag.comp);');
        set(figureHandle, 'color', 'white');
        drawnow;

        % Show stacked ICA/Ocular/ECG in a separate figure
        figureHandle2 = plot_ICAcomp2remove_EOG_ECG(compFlag, tmp);
        drawnow;

        % Prompt user for confirmation
        fprintf('\n---------- Selected Components for Removal ----------\n');
        disp(tmp);
        confirmStr = input('Definitely remove these components?\n Type "yes" to confirm or "no" to re-select: ', 's');

        % Close figures
        if ishandle(figureHandle), close(figureHandle); end
        if ishandle(figureHandle2),   close(figureHandle2);   end

        if strcmpi(confirmStr, 'yes')
            confirmed = true;
            comp2removeFinal = tmp;
        else
            % Restart inspection if the selected components are not confirmed
            fprintf('Restarting the entire ICA component inspection...\n\n');
            comp2removeFinal = InspectICAcomponents(compFlag);  
            confirmed = true;
        end
    end
else
    comp2removeFinal = [];
end

close all;

end


