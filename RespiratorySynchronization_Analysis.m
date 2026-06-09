%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                 %%%
%%%              RESPIRATORY SYNC ANALYSIS          %%%
%%%                                                 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Step 1: Compute relative phase, circular mean phase, and PLV
% Step 2: Compute windowed cross-correlation indices
% Step 3: Export final respiratory synchronization table

% Synchronization performance outputs:
%   PLV
%   mean_relative_phase_deg
%   WCC_mean_maxCorr
%   WCC_mean_maxLag_sec

% Written by: Francesco Bubbico
% Last updated: June 2026

%% Step 1: Relative Phase and Phase Locking Values

clear; clc;

addpath('C:\Users\francescob\Desktop\Open Data');
addpath('C:\Users\francescob\Desktop\Open Data\Code\Functions');

mainDir = 'C:\Users\francescob\Desktop\Open Data\Data\Processed data\RespiratorySync';

participantFile = 'C:\Users\francescob\Desktop\Open Data\Data\Participants.xlsx';

datasets = { ...
    'RespSync_LP_data', ...
    'RespSync_HP_data', ...
    'RespSync_PAT_data', ...
    'RespSync_HEA_data'};

% Load participant information
C = readcell(participantFile);
headers = matlab.lang.makeValidName(string(C(1,:)));
participantTable = cell2table(C(2:end,:), 'VariableNames', cellstr(headers));
participantTable.ID = str2double(string(participantTable.ID));
participantTable.MatlabID = arrayfun(@(x) sprintf('ID_%03d', x), ...
    participantTable.ID, 'UniformOutput', false);

% Compute metrics
for d = 1:numel(datasets)
    datasetName = datasets{d};
    filePath = fullfile(mainDir, [datasetName '.mat']);
    if ~isfile(filePath)
        warning('File not found: %s', filePath);
        continue
    end
    S = load(filePath, datasetName);
    if ~isfield(S, datasetName)
        warning('Variable %s not found in %s', datasetName, filePath);
        continue
    end
    data = S.(datasetName);
    participantIDs = fieldnames(data);
    fprintf('\n===== Processing %s =====\n', datasetName);

    for p = 1:numel(participantIDs)
        participantID = participantIDs{p};

        % Assign group
        match = strcmp(participantTable.MatlabID, participantID);
        if any(match)
            data.(participantID).group = string(participantTable.gruppo{find(match, 1)});
        else
            data.(participantID).group = "";
            warning('Group not found for %s', participantID);
        end

        % Compute relative phase and PLV
        requiredFields = {'participant_phase_masked', 'stimulus_phase'};
        if ~all(isfield(data.(participantID), requiredFields))
            warning('Missing phase fields for %s in %s', participantID, datasetName);
            continue
        end
        participant_phase = data.(participantID).participant_phase_masked(:)';
        stimulus_phase    = data.(participantID).stimulus_phase(:)';

        % Match length as a safety check
        n = min(numel(participant_phase), numel(stimulus_phase));
        participant_phase = participant_phase(1:n);
        stimulus_phase    = stimulus_phase(1:n);

        valid = isfinite(participant_phase) & isfinite(stimulus_phase);

        relative_phase_rad = NaN(1, n);
        relative_phase_deg = NaN(1, n);
        mean_relative_phase_rad = NaN;
        mean_relative_phase_deg = NaN;
        PLV = NaN;

        if any(valid)

            % Relative phase:
            % positive values = participant leads stimulus
            % negative values = participant lags stimulus
            dphi = angle(exp(1i * ...
                (participant_phase(valid) - stimulus_phase(valid))));
            relative_phase_rad(valid) = dphi;
            relative_phase_deg(valid) = rad2deg(dphi);

            % Circular mean relative phase
            mean_relative_phase_rad = angle(mean(exp(1i * dphi)));
            mean_relative_phase_deg = rad2deg(mean_relative_phase_rad);

            % Phase-locking value
            PLV = abs(mean(exp(1i * dphi)));
        end

        data.(participantID).relative_phase_rad = relative_phase_rad;
        data.(participantID).relative_phase_deg = relative_phase_deg;
        data.(participantID).mean_relative_phase_rad = mean_relative_phase_rad;
        data.(participantID).mean_relative_phase_deg = mean_relative_phase_deg;
        data.(participantID).PLV = PLV;
    end

    outStruct = struct();
    outStruct.(datasetName) = data;
    save(filePath, '-struct', 'outStruct', '-v7.3');
    fprintf('Saved updated file: %s\n', filePath);
end


% Plot group-level relative-phase distributions
plotDir = 'C:\Users\francescob\Desktop\Open Data\Plots\RelativePhaseDistribution';
if ~exist(plotDir, 'dir'); mkdir(plotDir); end

edges = -180:2:180;
binCenters = edges(1:end-1) + diff(edges)/2;

comparisonPairs = { ...
    {'RespSync_HP_data',  'RespSync_LP_data',  'High predictability', 'Low predictability'} ...
    {'RespSync_HEA_data', 'RespSync_PAT_data', 'Healthy bias',        'Patient bias'}};

groups = {'CONTROLLI','PAZIENTI'};
groupDisplayMap = containers.Map({'CONTROLLI','PAZIENTI'}, {'Controls','Patients'});

for g = 1:numel(groups)
    groupName = groups{g};
    groupLabel = groupDisplayMap(groupName);

    for c = 1:numel(comparisonPairs)
        dataset1 = comparisonPairs{c}{1};
        dataset2 = comparisonPairs{c}{2};
        label1   = comparisonPairs{c}{3};
        label2   = comparisonPairs{c}{4};

        S1 = load(fullfile(mainDir, [dataset1 '.mat']), dataset1);
        S2 = load(fullfile(mainDir, [dataset2 '.mat']), dataset2);
        data1 = S1.(dataset1);
        data2 = S2.(dataset2);

        participantIDs = intersect(fieldnames(data1), fieldnames(data2));
        hists1 = [];
        hists2 = [];

        for p = 1:numel(participantIDs)
            participantID = participantIDs{p};

            hasGroup = isfield(data1.(participantID), 'group') && isfield(data2.(participantID), 'group');
            hasPhase = isfield(data1.(participantID), 'relative_phase_deg') && isfield(data2.(participantID), 'relative_phase_deg');
            if ~hasGroup || ~hasPhase; continue; end

            if string(data1.(participantID).group) ~= string(groupName) || ...
               string(data2.(participantID).group) ~= string(groupName)
                continue
            end

            rp1 = data1.(participantID).relative_phase_deg;
            rp2 = data2.(participantID).relative_phase_deg;
            rp1 = rp1(isfinite(rp1));
            rp2 = rp2(isfinite(rp2));
            if isempty(rp1) || isempty(rp2); continue; end

            counts1 = histcounts(rp1, edges);
            counts2 = histcounts(rp2, edges);
            hists1 = [hists1; 100 * counts1 / sum(counts1)]; %#ok<AGROW>
            hists2 = [hists2; 100 * counts2 / sum(counts2)]; %#ok<AGROW>
        end

        if isempty(hists1) || isempty(hists2)
            warning('No valid data for %s vs %s, group %s', dataset1, dataset2, groupName);
            continue
        end

        mean1 = mean(hists1, 1, 'omitnan');
        mean2 = mean(hists2, 1, 'omitnan');
        sem1  = std(hists1, [], 1, 'omitnan') ./ sqrt(size(hists1,1));
        sem2  = std(hists2, [], 1, 'omitnan') ./ sqrt(size(hists2,1));

        fig = figure('Color','w','Position',[100 100 1200 700]); hold on;

        fill([binCenters fliplr(binCenters)], [mean1-sem1 fliplr(mean1+sem1)], ...
            [0.10 0.45 0.85], 'FaceAlpha',0.18, 'EdgeColor','none');
        p1 = plot(binCenters, mean1, 'Color',[0.10 0.45 0.85], 'LineWidth',2.5);

        fill([binCenters fliplr(binCenters)], [mean2-sem2 fliplr(mean2+sem2)], ...
            [0.90 0.45 0.10], 'FaceAlpha',0.18, 'EdgeColor','none');
        p2 = plot(binCenters, mean2, 'Color',[0.90 0.45 0.10], 'LineWidth',2.5);

        xline(0, '--k', 'LineWidth',1.1);
        xline(-180, ':k', 'LineWidth',1.0);
        xline(180, ':k', 'LineWidth',1.0);

        xlabel('Relative phase (degrees)');
        ylabel('Occurrence (%)');
        title(sprintf('Relative-phase distribution: %s', groupLabel), 'Interpreter','none');
        subtitle(sprintf('%s vs %s | n = %d', label1, label2, size(hists1,1)));
        legend([p1 p2], {label1, label2}, 'Location','best', 'Box','off');

        xlim([-180 180]);
        grid on; box off;
        set(gca, 'FontSize',12, 'LineWidth',1.1);

        saveName = sprintf('RelativePhaseDistribution_%s_vs_%s_%s.png', dataset1, dataset2, groupLabel);
        saveName = strrep(saveName, ' ', '');
        exportgraphics(fig, fullfile(plotDir, saveName), 'Resolution',600);
        close(fig);

        fprintf('Saved plot: %s\n', fullfile(plotDir, saveName));
    end
end

%% Step 2: Windowed cross-correlation
% Computes WCC between participant and stimulus respiration.
%
% Outputs:
%   mean_maxCorr     : average local synchronization strength at optimal lag
%   mean_maxLag_sec  : average signed lag at strongest coupling
%
% Lag interpretation:
%   Positive lag = participant lags behind stimulus
%   Negative lag = participant leads / anticipates stimulus

clear; clc;

mainDir = 'C:\Users\francescob\Desktop\Open Data\Data\Processed data\RespiratorySync';

datasets = {'RespSync_LP_data','RespSync_HP_data','RespSync_PAT_data','RespSync_HEA_data'};

fs = 500;              % Hz
window_size_sec = 8;   % sliding-window length
max_lag_sec     = 2;   % maximum lag tested
window_inc_sec  = 2;   % window step
lag_inc_sec     = 0.05;

window_size = round(window_size_sec * fs);
max_lag     = round(max_lag_sec * fs);
window_inc  = round(window_inc_sec * fs);
lag_inc     = max(1, round(lag_inc_sec * fs));

for d = 1:numel(datasets)
    datasetName = datasets{d};
    filePath = fullfile(mainDir, [datasetName '.mat']);

    if ~isfile(filePath)
        warning('File not found: %s', filePath);
        continue
    end

    S = load(filePath, datasetName);
    if ~isfield(S, datasetName)
        warning('Variable %s not found in %s', datasetName, filePath);
        continue
    end

    data = S.(datasetName);
    participantIDs = fieldnames(data);
    fprintf('\n===== WCC: %s =====\n', datasetName);
    
    for p = 1:numel(participantIDs)
        participantID = participantIDs{p};
        requiredFields = {'participant_resp_masked','stimulus_resp'};

        if ~all(isfield(data.(participantID), requiredFields))
            warning('Missing respiration fields for %s in %s', participantID, datasetName);
            continue
        end

        participant_resp = data.(participantID).participant_resp_masked(:);
        stimulus_resp    = data.(participantID).stimulus_resp(:);

        n = min(numel(participant_resp), numel(stimulus_resp));
        participant_resp = participant_resp(1:n);
        stimulus_resp    = stimulus_resp(1:n);

        if n < window_size
            warning('Signal too short for WCC: %s in %s', participantID, datasetName);
            continue
        end

        nWindows = floor((n - window_size) / window_inc) + 1;
        window_maxCorr = NaN(nWindows,1);
        window_maxLag_sec = NaN(nWindows,1);

        for w = 1:nWindows

            idx1 = (w - 1) * window_inc + 1;
            idx2 = idx1 + window_size - 1;
            x = participant_resp(idx1:idx2);
            y = stimulus_resp(idx1:idx2);

            if any(~isfinite(x)) || any(~isfinite(y))
                continue
            end

            [xc, lags] = xcorr(x, y, max_lag, 'coeff');
            xc = xc(1:lag_inc:end);
            lags = lags(1:lag_inc:end);
            [window_maxCorr(w), idxMax] = max(xc);
            window_maxLag_sec(w) = lags(idxMax) / fs;
        end

        data.(participantID).WCC.mean_maxCorr = mean(window_maxCorr, 'omitnan');
        data.(participantID).WCC.mean_maxLag_sec = mean(window_maxLag_sec, 'omitnan');
    end

    outStruct = struct();
    outStruct.(datasetName) = data;
    save(filePath, '-struct', 'outStruct', '-v7.3');
    fprintf('Saved updated file: %s\n', filePath);
end

%% Step 3: Export respiratory synchronization summary table

clear; clc;

inputDir  = 'C:\Users\francescob\Desktop\Open Data\Data\Processed data\RespiratorySync';
outputDir = 'C:\Users\francescob\Desktop\Open Data\Data\Final datasets';

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

datasets = {'RespSync_LP_data','RespSync_HP_data','RespSync_PAT_data','RespSync_HEA_data'};

results = struct( ...
    'ID', {}, ...
    'Condition', {}, ...
    'group', {}, ...
    'PLV', {}, ...
    'mean_relative_phase_deg', {}, ...
    'WCC_mean_maxCorr', {}, ...
    'WCC_mean_maxLag_sec', {} );

idx = 1;

for d = 1:numel(datasets)

    datasetName = datasets{d};
    filePath = fullfile(inputDir, [datasetName '.mat']);

    if ~isfile(filePath)
        warning('File not found: %s', filePath);
        continue
    end

    S = load(filePath, datasetName);

    if ~isfield(S, datasetName)
        warning('Variable %s not found in %s', datasetName, filePath);
        continue
    end

    data = S.(datasetName);
    participantIDs = fieldnames(data);

    for p = 1:numel(participantIDs)
        participantID = participantIDs{p};
        numID = sscanf(participantID, 'ID_%d');
        results(idx).ID = sprintf('%03d', numID);
        results(idx).Condition = datasetName;

        if isfield(data.(participantID), 'group')
            results(idx).group = data.(participantID).group;
        else
            results(idx).group = "";
        end

        if isfield(data.(participantID), 'PLV')
            results(idx).PLV = data.(participantID).PLV;
        else
            results(idx).PLV = NaN;
        end

        if isfield(data.(participantID), 'mean_relative_phase_deg')
            results(idx).mean_relative_phase_deg = data.(participantID).mean_relative_phase_deg;
        else
            results(idx).mean_relative_phase_deg = NaN;
        end

        if isfield(data.(participantID), 'WCC')
            W = data.(participantID).WCC;

            if isfield(W, 'mean_maxCorr')
                results(idx).WCC_mean_maxCorr = W.mean_maxCorr;
            else
                results(idx).WCC_mean_maxCorr = NaN;
            end

            if isfield(W, 'mean_maxLag_sec')
                results(idx).WCC_mean_maxLag_sec = W.mean_maxLag_sec;
            else
                results(idx).WCC_mean_maxLag_sec = NaN;
            end
        else
            results(idx).WCC_mean_maxCorr = NaN;
            results(idx).WCC_mean_maxLag_sec = NaN;
        end

        idx = idx + 1;
    end
end

summaryTable = struct2table(results);
outFile = fullfile(outputDir, 'RespiratorySync_dataset.csv');
writetable(summaryTable, outFile);
fprintf('\nRespiratory synchronization table exported to:\n%s\n', outFile);
