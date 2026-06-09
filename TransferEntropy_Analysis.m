%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                 %%%
%%%              TRANSFER ENTROPY ANALYSIS          %%%
%%%                                                 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Computes transfer entropy between heart, respiration,
% and EEG-derived signals.

% Step 1 computes subject-level transfer entropy estimates.
% Step 2 performs global averaging, cluster-based permutation testing,
% and result export.

% Written by: Francesco Bubbico
% Last updated: June 2026

%% Step 1: Network Physiology analysis - Transfer Entropy

clear; clc;

% Parallel + toolbox setup
toolboxDir = 'C:\Program Files\MATLAB\R2024a\toolbox\Network Physiology\TOOLS\ITScode_v2.1\ITScode_v2.1\functions';
addpath(toolboxDir);

assert(~isempty(which('its_BTElin'))   && ~isempty(which('its_BTEknn'))   && ...
    ~isempty(which('its_BTElin_V')) && ~isempty(which('its_BTEknn_V')) && ...
    ~isempty(which('its_SetLag')), ...
    'ITS TE functions not on path');

delete(gcp('nocreate'));
parpool('local');
pctRunOnAll addpath('C:\Program Files\MATLAB\R2024a\toolbox\Network Physiology\TOOLS\ITScode_v2.1\ITScode_v2.1\functions');
pctRunOnAll rehash path

% Select manually datasets to process (paired conditions)
datasets = [ ...
    struct('dataFolder','C:\Users\francescob\Desktop\Open Data\Data\Processed data\BodyBrain_2HzTimeSeries\hp', 'selectedCond','hp'), ...
    struct('dataFolder','C:\Users\francescob\Desktop\Open Data\Data\Processed data\BodyBrain_2HzTimeSeries\lp', 'selectedCond','lp') ...
    %struct('dataFolder','C:\Users\francescob\Desktop\Open Data\Data\Processed data\BodyBrain_2HzTimeSeries\hea', 'selectedCond','hea'), ...
    %struct('dataFolder','C:\Users\francescob\Desktop\Open Data\Data\Processed data\BodyBrain_2HzTimeSeries\pat', 'selectedCond','pat'), ...
    ];

% Output folder
outputDir = 'C:\Users\francescob\Desktop\Open Data\Data\Final datasets';

% ============================================================
% TRANSFER ENTROPY ANALYSIS
% EEG <-> RESP and EEG <-> RR
% ============================================================

% Parameters
p        = 4;    % fixed embedding lags 1:p
kNN      = 10;   % k-nearest-neighbour estimator
minValid = 50;   % minimum valid samples

bandNames = {'delta','theta','alpha','alpha1','alpha2', ...
             'beta','beta1','beta2','gamma'};

for iDS = 1:numel(datasets)

    dataFolder   = datasets(iDS).dataFolder;
    selectedCond = datasets(iDS).selectedCond;

    dsName  = lower(string(selectedCond));
    condStr = string(selectedCond);

    fprintf('\n=== DATASET %d/%d: %s ===\n', iDS, numel(datasets), dsName);

    files = dir(fullfile(dataFolder, 'networkPhysiology_*.mat'));

    if isempty(files)
        fprintf('No Network Physiology files found — skipping.\n');
        continue;
    end

    % Output table
    data_TE = table( ...
        strings(0,1), strings(0,1), ...
        zeros(0,1), strings(0,1), strings(0,1), ...
        nan(0,1), nan(0,1), ...
        'VariableNames', {'subjectID','selectedCond', ...
                          'channel','chanLabel','interaction', ...
                          'te_lin','te_knn'});

    tStart = tic;

    % ========================================================
    % SUBJECT LOOP
    % ========================================================
    for iFile = 1:numel(files)

        inPath = fullfile(files(iFile).folder, files(iFile).name);
        [~, baseName] = fileparts(inPath);
        subjectID = erase(baseName, 'networkPhysiology_');

        S = load(inPath);
        data = S;

        % -------------------------
        % Extract physiology
        % -------------------------
        RR_ms      = 1000 * data.np.RRI_s(:);
        RESP_amp   = data.np.resp_amp(:);
        RESP_phase = data.np.resp_phase(:);

        % Respiration model:
        % amplitude + sin(phase) + cos(phase)
        RESP_data = [ ...
            RESP_amp, ...
            sin(RESP_phase), ...
            cos(RESP_phase) ...
        ];

        % -------------------------
        % Extract EEG band power
        % -------------------------
        labels = string(data.np.peeg.labels(:));
        P      = data.np.peeg.power_per_band;

        nCh    = numel(labels);
        nBands = numel(bandNames);

        fprintf('Subject %d/%d: %s | channels=%d\n', ...
            iFile, numel(files), subjectID, nCh);

        outCell = cell(nCh,1);

        % ====================================================
        % CHANNEL LOOP
        % ====================================================
        parfor ch = 1:nCh

            EEG_cols = zeros(numel(RESP_amp), nBands);

            for b = 1:nBands
                EEG_cols(:,b) = P.(bandNames{b})(ch,:).';
            end

            subj    = string(subjectID);
            chanLab = labels(ch);

            subjTab_ch = table();

            for b = 1:nBands

                EEGk = EEG_cols(:,b);

                % =====================================================
                % EEG <-> RESP
                % SERIES = [RESP_amp RESP_sin RESP_cos EEG]
                % RESP indices = 1:3
                % EEG index    = 4
                % =====================================================
                SERIES_er = [RESP_data EEGk];
                valid_er  = all(isfinite(SERIES_er), 2);

                if nnz(valid_er) >= minValid

                    Yo_er = SERIES_er(valid_er,:);

                    % Linear TE: demean only
                    Y_er = zeros(size(Yo_er));
                    for m = 1:size(Yo_er,2)
                        Y_er(:,m) = Yo_er(:,m) - mean(Yo_er(:,m));
                    end

                    % kNN TE: z-score
                    Yz_er = zscore(Yo_er, 0, 1);

                    nResp  = size(RESP_data,2);
                    respIdx = 1:nResp;
                    eegIdx  = nResp + 1;
                    M_er    = size(Yo_er,2);

                    % -------------------------
                    % EEG -> RESP
                    % interaction label: band_resp
                    % -------------------------
                    ii = eegIdx;
                    jj = respIdx;

                    % Linear TE
                    if numel(jj) == 1
                        zerolag = zeros(1,M_er);
                        zerolag(ii) = 1;

                        pV = p * ones(1,M_er);
                        VL = its_SetLag(pV, ones(1,M_er), ones(1,M_er), zerolag);

                        out_lin = its_BTElin(Y_er, ii, jj, VL);
                        te_lin = out_lin.Txy;
                    else
                        out_lin = its_BTElin_V(Y_er, ii, jj, p);
                        te_lin = out_lin.TE;
                    end

                    % kNN TE
                    zerolag_knn = zeros(1,M_er);
                    zerolag_knn(ii) = 1;

                    pV_knn = p * ones(1,M_er);
                    VL_knn = its_SetLag(pV_knn, ones(1,M_er), ones(1,M_er), zerolag_knn);

                    if numel(jj) == 1
                        out_knn = its_BTEknn(Yz_er, VL_knn, ii, jj, kNN);
                    else
                        out_knn = its_BTEknn_V(Yz_er, VL_knn, ii, jj, kNN);
                    end

                    te_knn = out_knn.Txy;

                    row_eeg_resp = table( ...
                        subj, condStr, ...
                        ch, chanLab, strcat(string(bandNames{b}), "_resp"), ...
                        te_lin, te_knn, ...
                        'VariableNames', data_TE.Properties.VariableNames);

                    % -------------------------
                    % RESP -> EEG
                    % interaction label: resp_band
                    % -------------------------
                    ii = respIdx;
                    jj = eegIdx;

                    % Linear TE
                    zerolag = zeros(1,M_er);
                    zerolag(ii) = 1;

                    pV = p * ones(1,M_er);
                    VL = its_SetLag(pV, ones(1,M_er), ones(1,M_er), zerolag);

                    out_lin = its_BTElin(Y_er, ii, jj, VL);
                    te_lin = out_lin.Txy;

                    % kNN TE
                    zerolag_knn = zeros(1,M_er);
                    zerolag_knn(ii) = 1;

                    pV_knn = p * ones(1,M_er);
                    VL_knn = its_SetLag(pV_knn, ones(1,M_er), ones(1,M_er), zerolag_knn);

                    out_knn = its_BTEknn(Yz_er, VL_knn, ii, jj, kNN);
                    te_knn = out_knn.Txy;

                    row_resp_eeg = table( ...
                        subj, condStr, ...
                        ch, chanLab, strcat("resp_", string(bandNames{b})), ...
                        te_lin, te_knn, ...
                        'VariableNames', data_TE.Properties.VariableNames);

                    subjTab_ch = [subjTab_ch; row_eeg_resp; row_resp_eeg];
                end

                % =====================================================
                % EEG <-> RR
                % SERIES = [RR EEG]
                % RR index  = 1
                % EEG index = 2
                % =====================================================
                SERIES_rr = [RR_ms EEGk];
                valid_rr  = all(isfinite(SERIES_rr), 2);

                if nnz(valid_rr) >= minValid

                    Yo_rr = SERIES_rr(valid_rr,:);

                    % Linear TE: demean only
                    Y_rr = zeros(size(Yo_rr));
                    for m = 1:size(Yo_rr,2)
                        Y_rr(:,m) = Yo_rr(:,m) - mean(Yo_rr(:,m));
                    end

                    % kNN TE: z-score
                    Yz_rr = zscore(Yo_rr, 0, 1);

                    M_rr = 2;

                    % -------------------------
                    % EEG -> RR
                    % interaction label: band_rr
                    % -------------------------
                    ii = 2;   % EEG driver
                    jj = 1;   % RR target

                    % Linear TE
                    zerolag = zeros(1,M_rr);
                    zerolag(ii) = 1;

                    pV = p * ones(1,M_rr);
                    VL = its_SetLag(pV, ones(1,M_rr), ones(1,M_rr), zerolag);

                    out_lin = its_BTElin(Y_rr, ii, jj, VL);
                    te_lin = out_lin.Txy;

                    % kNN TE
                    zerolag_knn = zeros(1,M_rr);
                    zerolag_knn(ii) = 1;

                    pV_knn = p * ones(1,M_rr);
                    VL_knn = its_SetLag(pV_knn, ones(1,M_rr), ones(1,M_rr), zerolag_knn);

                    out_knn = its_BTEknn(Yz_rr, VL_knn, ii, jj, kNN);
                    te_knn = out_knn.Txy;

                    row_eeg_rr = table( ...
                        subj, condStr, ...
                        ch, chanLab, strcat(string(bandNames{b}), "_rr"), ...
                        te_lin, te_knn, ...
                        'VariableNames', data_TE.Properties.VariableNames);

                    % -------------------------
                    % RR -> EEG
                    % interaction label: rr_band
                    % -------------------------
                    ii = 1;   % RR driver
                    jj = 2;   % EEG target

                    % Linear TE
                    zerolag = zeros(1,M_rr);
                    zerolag(ii) = 1;

                    pV = p * ones(1,M_rr);
                    VL = its_SetLag(pV, ones(1,M_rr), ones(1,M_rr), zerolag);

                    out_lin = its_BTElin(Y_rr, ii, jj, VL);
                    te_lin = out_lin.Txy;

                    % kNN TE
                    zerolag_knn = zeros(1,M_rr);
                    zerolag_knn(ii) = 1;

                    pV_knn = p * ones(1,M_rr);
                    VL_knn = its_SetLag(pV_knn, ones(1,M_rr), ones(1,M_rr), zerolag_knn);

                    out_knn = its_BTEknn(Yz_rr, VL_knn, ii, jj, kNN);
                    te_knn = out_knn.Txy;

                    row_rr_eeg = table( ...
                        subj, condStr, ...
                        ch, chanLab, strcat("rr_", string(bandNames{b})), ...
                        te_lin, te_knn, ...
                        'VariableNames', data_TE.Properties.VariableNames);

                    subjTab_ch = [subjTab_ch; row_eeg_rr; row_rr_eeg];
                end
            end

            outCell{ch} = subjTab_ch;
        end

        subjTab = vertcat(outCell{:});
        data_TE = [data_TE; subjTab];

        fprintf('  done subject %d/%d | elapsed: %.1fs\n', ...
            iFile, numel(files), toc(tStart));
    end

    % Save output
    writetable(data_TE, fullfile(outputDir, sprintf('TE_%s.csv', dsName)));
    save(fullfile(outputDir, sprintf('TE_%s.mat', dsName)), 'data_TE');

    fprintf('Saved TE_%s | elapsed %.1fs\n', dsName, toc(tStart));

    clear data_TE files
end


%% Step 2: Global averages, cluster-based permutation and export tables

clear; clc;

addpath('C:\Users\francescob\Desktop\Open Data\Code\Functions');
addpath('C:\Program Files\MATLAB\R2024a\toolbox\Network Physiology\TOOLS\ITScode_v2.1\ITScode_v2.1\functions');

% Directories 
mainDir   = 'C:\Users\francescob\Desktop\Open Data';
outputDir = fullfile(mainDir, 'Data', 'Final datasets');
plotDir   = fullfile(mainDir, 'Plots', 'TE_Maps');
if ~exist(outputDir,'dir'), mkdir(outputDir); end
if ~exist(plotDir,'dir'),   mkdir(plotDir);   end


% 2.1) -------- SETUP --------

% Inputs: manually select two files at time (paired conditions)
%file1 = fullfile(outputDir, 'TE_hp.mat');  cond1 = 'hp';
%file2 = fullfile(outputDir, 'TE_lp.mat');  cond2 = 'lp';
file1 = fullfile(outputDir, 'TE_hea.mat');  cond1 = 'hea';
file2 = fullfile(outputDir, 'TE_pat.mat');  cond2 = 'pat';

electrodeFile    = 'C:\Users\francescob\Desktop\Open Data\Code\Functions\Standard-10-20-Cap62+ref_EDIT.txt';
participantsFile = 'C:\Users\francescob\Desktop\Open Data\Data\Participants.xlsx';

% Interactions:
% all EEG<->RR and EEG<->RESP interactions
interactions = { ...
    'delta_rr','theta_rr','alpha_rr','alpha1_rr','alpha2_rr','beta_rr','beta1_rr','beta2_rr','gamma_rr', ...
    'rr_delta','rr_theta','rr_alpha','rr_alpha1','rr_alpha2','rr_beta','rr_beta1','rr_beta2','rr_gamma', ...
    'delta_resp','theta_resp','alpha_resp','alpha1_resp','alpha2_resp','beta_resp','beta1_resp','beta2_resp','gamma_resp', ...
    'resp_delta','resp_theta','resp_alpha','resp_alpha1','resp_alpha2','resp_beta','resp_beta1','resp_beta2','resp_gamma', ...
    };

interactions = string(interactions(:));

inter_rr_out   = interactions(endsWith(interactions,"_rr"));
inter_rr_in    = interactions(startsWith(interactions,"rr_"));
inter_resp_out = interactions(endsWith(interactions,"_resp"));
inter_resp_in  = interactions(startsWith(interactions,"resp_"));

interactions_to_use = interactions;   % use to subset to specific interactions; here we use all

% Load and prepare data
% Channel location
L = readtable(electrodeFile,'FileType','text','Delimiter',{' ','\t'},'MultipleDelimsAsOne',true);
locations = table(string(L.labels), L.theta, L.radius, 'VariableNames',{'labels','theta','radius'});
locations.labels = strtrim(string(locations.labels));

% Participants and groups
P = readtable(participantsFile);
P.subjectID = compose('%03d', str2double(regexprep(string(P.ID), '\D', '')));

% Load TE data
S1 = load(file1);
S2 = load(file2);

% Tag
conditionTag = sprintf('%s_%s', cond1, cond2);

% Normalize subject IDs
S1.data_TE.subjectID = compose('%03d', str2double(regexprep(string(S1.data_TE.subjectID), '\D', '')));
S2.data_TE.subjectID = compose('%03d', str2double(regexprep(string(S2.data_TE.subjectID), '\D', '')));

% Keep only selected interactions
S1.data_TE = S1.data_TE(ismember(string(S1.data_TE.interaction), interactions_to_use), :);
S2.data_TE = S2.data_TE(ismember(string(S2.data_TE.interaction), interactions_to_use), :);

% Attach gruppo (group)
S1.data_TE = outerjoin(S1.data_TE, P(:,{'subjectID','gruppo'}), 'Keys','subjectID','MergeKeys',true,'Type','left');
S2.data_TE = outerjoin(S2.data_TE, P(:,{'subjectID','gruppo'}), 'Keys','subjectID','MergeKeys',true,'Type','left');

% Relabel gruppo (Italian to English)
if ~iscategorical(S1.data_TE.gruppo), S1.data_TE.gruppo = categorical(string(S1.data_TE.gruppo)); end
if ~iscategorical(S2.data_TE.gruppo), S2.data_TE.gruppo = categorical(string(S2.data_TE.gruppo)); end
old = {'PAZIENTI','CONTROLLI'};
new = {'Patients','Controls'};
c = categories(S1.data_TE.gruppo); tf = ismember(old,c); if any(tf), S1.data_TE.gruppo = renamecats(S1.data_TE.gruppo, old(tf), new(tf)); end
c = categories(S2.data_TE.gruppo); tf = ismember(old,c); if any(tf), S2.data_TE.gruppo = renamecats(S2.data_TE.gruppo, old(tf), new(tf)); end
S1.data_TE.Properties.VariableNames{'gruppo'} = 'group';
S2.data_TE.Properties.VariableNames{'gruppo'} = 'group';

% Groups order
gNames = unique([string(categories(S1.data_TE.group)); string(categories(S2.data_TE.group))]);
gNames = gNames(~strcmp(gNames,"<undefined>"));
if all(ismember(["Controls","Patients"], gNames))
    gNames = ["Controls","Patients"];
end

% FieldTrip neighbours
% Montage = locations ∩ data channels (ONLY channels in data) 
dataCh = unique(strtrim(string([S1.data_TE.chanLabel; S2.data_TE.chanLabel])));
dataCh(strcmpi(dataCh,"global")) = [];
locAll = locations;
locAll.labels = strtrim(string(locAll.labels));
useMask = ismember(locAll.labels, dataCh);
locTab  = locAll(useMask, :);
labels_all = string(locTab.labels);
nCh = numel(labels_all);

addpath('C:\Program Files\MATLAB\R2024a\toolbox\fieldtrip-20250106');
ft_defaults;

% 2D projection 
[chx, chy] = pol2cart(deg2rad((-1).*locTab.theta + 90), locTab.radius);

% FieldTrip layout struct 
lay = [];
lay.label  = cellstr(labels_all);     % cell array of char vectors
lay.pos    = [chx(:) chy(:)];         % Nx2

% Width/height for FieldTrip layout visualizations
lay.width  = 0.03 * ones(nCh,1);
lay.height = 0.03 * ones(nCh,1);

% Compute neighbours via triangulation
cfg = [];
cfg.method  = 'triangulation';
cfg.layout  = lay;
neigh = ft_prepare_neighbours(cfg);

% Convert FieldTrip neighbours to adjacency matrix
A = false(nCh);
lab2idx = containers.Map(cellstr(labels_all), 1:nCh);

for i = 1:numel(neigh)
    if ~isKey(lab2idx, neigh(i).label), continue; end
    ii = lab2idx(neigh(i).label);

    nb = neigh(i).neighblabel;
    for j = 1:numel(nb)
        if ~isKey(lab2idx, nb{j}), continue; end
        jj = lab2idx(nb{j});
        A(ii,jj) = true;
        A(jj,ii) = true;  % enforce symmetry
    end
end

% Visualize neighbours 
cfg = [];
cfg.layout     = lay;
cfg.neighbours = neigh;
ft_neighbourplot(cfg);

% Check number of neighbours
deg_ft = sum(A,2);
disp(table(labels_all, deg_ft));


% 2.2) -------- ADD GLOBAL ELECTRODES AVERAGE AND EXPORT IN LONG FORMAT --------

T_electrodes = export_electrode_table_NP( ...
    S1.data_TE, S2.data_TE, interactions_to_use, labels_all, ...
    'metricVars', {'te_lin','te_knn'}, ...
    'addGlobal', true);

writetable(T_electrodes, fullfile(outputDir, sprintf('TE_electrode_%s.csv', conditionTag)));


% 2.3) -------- CLUSTER BASED PERMUTATION --------
% (Group difference in condition-related TE changes)
%
% For each interaction and TE metric (te_lin, te_knn), subject-level
% change scores are computed as:
%
%       delta = TE_condition2 − TE_condition1
%
% At each channel, delta values are compared between groups using
% an independent-samples t-test. Channels exceeding an uncorrected
% threshold (alphaCluster = 0.05) are spatially clustered using the
% electrode adjacency matrix.
%
% For each cluster, cluster mass is computed as the sum of t-values.
% A permutation null distribution (nPerm = 5000) is generated by
% randomly shuffling group labels and storing the maximum cluster
% mass observed across the scalp at each iteration.
%
% Cluster-level p-values are computed as the proportion of permutations
% in which the maximum null cluster mass exceeds the observed cluster mass.
%
% This procedure controls the family-wise error rate across channels
% and provides spatially corrected inference for group differences in
% condition-related TE changes.

ClusterPerm = cluster_perm_interaction_NP( ...
    S1.data_TE, S2.data_TE, interactions_to_use, labels_all, A, ...
    'metricVars', {'te_lin','te_knn'}, ...
    'groupNames', gNames, ...
    'nPerm', 5000, ...
    'alphaCluster', 0.05, ...
    'alphaCorrected', 0.05, ...
    'tail', 'two', ...
    'minCluster', 2, ...
    'seed', 1, ...
    'verbose', true);

save(fullfile(outputDir, ['clusterPerm_TE_' conditionTag '.mat']), 'ClusterPerm');

% Print significant clusters for linear and nonlinear TE
metrics = {'te_lin','te_knn'};
alpha = 0.05;

for m = 1:numel(metrics)
    metricName = metrics{m};
    R = ClusterPerm.(metricName);
    fprintf('\n\nSIGNIFICANT CLUSTERS: %s\n', metricName);
    fprintf('====================================\n');
    foundAny = false;
    for i = 1:numel(R)
        if isempty(R(i).clusters)
            continue
        end
        sigClusters = find([R(i).clusters.pCluster] < alpha);

        for c = sigClusters
            foundAny = true;
            cl = R(i).clusters(c);
            fprintf('\nInteraction: %s\n', R(i).interaction);
            fprintf('Cluster %d | sign = %+d | size = %d | mass = %.3f | p = %.4f\n', ...
                c, cl.sign, cl.size, cl.mass, cl.pCluster);
            fprintf('Electrodes: %s\n', strjoin(cellstr(cl.labels), ', '));
        end
    end
    if ~foundAny
        fprintf('No significant clusters found.\n');
    end
end


% Export cluster TE table
T_cluster = cluster_table_NP( ...
    ClusterPerm, S1.data_TE, S2.data_TE, labels_all, {'te_lin','te_knn'}, 0.05, ...
    'conditionTag', conditionTag);

writetable(T_cluster, fullfile(outputDir, sprintf('TE_cluster_%s.csv', conditionTag)));


% 2.4) -------- PLOT TE MAPS + CLUSTER --------

TEPlotDir = fullfile(plotDir, 'TE_per_interaction_withCluster');
if ~exist(TEPlotDir, 'dir')
    mkdir(TEPlotDir);
end

cfgTE = struct();
cfgTE.outDir = TEPlotDir;
cfgTE.conditionTag = conditionTag;
cfgTE.groupNames = gNames;
cfgTE.condNames  = [string(cond1) string(cond2)];
cfgTE.sigClusterAlpha = 0.05;
cfgTE.metricFields = {'te_lin','te_knn'};
cfgTE.metricLabelY = '\Delta TE';
cfgTE.diffLabelY   = '\Delta difference';

plot_maps_NP(ClusterPerm, labels_all, chx, chy, locTab, cfgTE);


