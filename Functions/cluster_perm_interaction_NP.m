function OUT = cluster_perm_interaction_NP(T1, T2, interactions_to_use, labels_all, A, varargin)
% Cluster-based permutation test on subject-level change scores.
%
% For each subject, condition-specific values are reduced to a change score:
%
%   delta = condition2 - condition1
%
% At each channel, change scores are compared between groups using
% Welch's independent-samples t-test. Spatially adjacent supra-threshold
% channels are combined into clusters, and cluster significance is
% evaluated using a permutation-based maximum cluster-mass statistic.
%
% Written by: Francesco Bubbico
% Last updated: June 2026
%
%
% Inputs
% ------
% T1, T2 : tables
%   Tables for condition 1 and condition 2.
%   Required columns:
%       - interaction   : interaction label
%       - chanLabel     : channel label
%       - subjectID     : subject identifier
%       - group         : group label
%       - metric vars   : e.g. mi_lin, mi_knn, te_lin, te_knn
%
% interactions_to_use : cellstr or string
%   Interaction labels to analyse, e.g. {'alpha_resp','resp_alpha'}
%
% labels_all : cellstr or string
%   Ordered list of channels used in the adjacency matrix
%
% A : logical/double matrix [nCh x nCh]
%   Symmetric channel adjacency matrix
%
%
% Options:
% ------------------
% 'metricVars'       : metric variables to analyse
%                      (default {'mi_lin','mi_knn'})
% 'groupNames'       : 2 group names in desired contrast order [g1 g2]
%                      effect is always group2 - group1
% 'nPerm'            : number of permutations (default 5000)
% 'alphaCluster'     : cluster-forming alpha (default 0.05)
% 'alphaCorrected'   : corrected cluster alpha (default 0.05)
% 'tail'             : 'two' | 'pos' | 'neg' (default 'two')
% 'minCluster'       : minimum channels in a cluster (default 2)
% 'seed'             : RNG seed (default 1)
% 'verbose'          : true/false (default true)
%
%
% Output
% ------
% OUT.(metricVar)(k)
%   Struct array, one entry per interaction.
%
%
% Notes
% -----
% - Condition-related effects are quantified as subject-level change scores:
%       delta = condition2 - condition1
% - Group differences are tested on these change scores.
% - Welch's t-test is performed independently at each channel.
% - Group labels are randomly permuted to generate the null distribution.
% - Cluster mass is computed as the sum of t-values within each cluster.
% - Statistical significance is assessed using the maximum cluster mass
%   observed across all permutations.
% --------------------------- options ---------------------------

p = inputParser;
p.addParameter('metricVars', {'mi_lin','mi_knn'}, @(x) iscell(x) || isstring(x));
p.addParameter('groupNames', [], @(x) isempty(x) || iscell(x) || isstring(x));
p.addParameter('nPerm', 5000, @(x) isnumeric(x) && isscalar(x) && x > 0);
p.addParameter('alphaCluster', 0.05, @(x) isnumeric(x) && isscalar(x) && x > 0 && x < 1);
p.addParameter('alphaCorrected', 0.05, @(x) isnumeric(x) && isscalar(x) && x > 0 && x < 1);
p.addParameter('tail', 'two', @(x) ischar(x) || isstring(x));
p.addParameter('minCluster', 2, @(x) isnumeric(x) && isscalar(x) && x >= 1);
p.addParameter('seed', 1, @(x) isnumeric(x) && isscalar(x));
p.addParameter('verbose', true, @(x) islogical(x) && isscalar(x));
p.parse(varargin{:});
cfg = p.Results;

% --------------------------- normalize ---------------------------
interactions_to_use = string(interactions_to_use(:));
labels_all          = upper(strtrim(string(labels_all(:))));
metricVars          = string(cfg.metricVars(:));
nCh                 = numel(labels_all);

requiredCols = ["interaction","chanLabel","subjectID","group"];
for c = requiredCols
    if ~ismember(c, string(T1.Properties.VariableNames))
        error('T1 is missing required column "%s".', c);
    end
    if ~ismember(c, string(T2.Properties.VariableNames))
        error('T2 is missing required column "%s".', c);
    end
end

T1 = normalize_input_table(T1);
T2 = normalize_input_table(T2);

% --------------------------- groups ---------------------------
gNames = string(cfg.groupNames(:));
if isempty(gNames)
    gNames = unique([T1.group; T2.group], 'stable');
end
gNames = lower(strtrim(gNames));

if numel(gNames) ~= 2
    error('Exactly 2 group names are required. Found: %s', strjoin(gNames, ', '));
end

g1 = gNames(1);
g2 = gNames(2);

% --------------------------- subjects in both conditions ---------------------------
subjBoth = intersect(unique(T1.subjectID), unique(T2.subjectID));
T1 = T1(ismember(T1.subjectID, subjBoth), :);
T2 = T2(ismember(T2.subjectID, subjBoth), :);

if isempty(T1) || isempty(T2)
    error('No overlapping subjects between the two conditions.');
end

% --------------------------- adjacency checks ---------------------------
A = logical(A);
if ~isequal(size(A), [nCh nCh])
    error('Adjacency matrix A must be %d x %d.', nCh, nCh);
end
A = A | A';
A(1:nCh+1:end) = false;

label2idx = containers.Map(cellstr(labels_all), num2cell(1:nCh));

rng(cfg.seed);

OUT = struct();

% =========================== main loop ===========================
for mv = 1:numel(metricVars)
    metricVar = metricVars(mv);

    if ~ismember(metricVar, string(T1.Properties.VariableNames)) || ...
       ~ismember(metricVar, string(T2.Properties.VariableNames))
        error('Metric variable "%s" not found in both T1 and T2.', metricVar);
    end

    RES = repmat(struct( ...
        'interaction', "", ...
        'metricVar', "", ...
        'groupNames', [g1 g2], ...
        'nSubjGroup1', NaN, ...
        'nSubjGroup2', NaN, ...
        'deltaMatrix', [], ...
        'groupVector', strings(0,1), ...
        'subjectIDs', strings(0,1), ...
        'tObs', nan(nCh,1), ...
        'dfObs', nan(nCh,1), ...
        'pObs_unc', nan(nCh,1), ...
        'tCrit', NaN, ...
        'supraMaskPos', false(nCh,1), ...
        'supraMaskNeg', false(nCh,1), ...
        'clusters', struct('idx',{},'labels',{},'sign',{},'mass',{},'size',{},'pCluster',{}), ...
        'sigClusterMask', false(nCh,1), ...
        'permMaxMass', nan(cfg.nPerm,1), ...
        'cfg', cfg), numel(interactions_to_use), 1);

    for kInt = 1:numel(interactions_to_use)
        interactionLabel = string(strtrim(interactions_to_use(kInt)));
        RES(kInt).interaction = interactionLabel;
        RES(kInt).metricVar   = metricVar;

        [Y, G, subjIDs, ok, msg] = build_delta_matrix( ...
            T1, T2, interactionLabel, metricVar, labels_all, label2idx, g1, g2);

        if ~ok
            if cfg.verbose
                warning('[%s | %s] %s', char(metricVar), char(interactionLabel), msg);
            end
            continue;
        end

        idx1 = (G == g1);
        idx2 = (G == g2);

        RES(kInt).deltaMatrix = Y;
        RES(kInt).groupVector = G;
        RES(kInt).subjectIDs  = subjIDs;
        RES(kInt).nSubjGroup1 = sum(idx1);
        RES(kInt).nSubjGroup2 = sum(idx2);

        % ---------- observed Welch t-map ----------
        [tObs, dfObs, pObs] = welch_t_map(Y, idx2, idx1);   % group2 - group1
        tObs  = tObs(:);
        dfObs = dfObs(:);
        pObs  = pObs(:);

        validDf = dfObs(isfinite(dfObs) & dfObs > 0);
        if isempty(validDf)
            if cfg.verbose
                warning('[%s | %s] no valid df for observed map.', char(metricVar), char(interactionLabel));
            end
            continue;
        end

        % Conservative single threshold across channels
        dfMin = min(validDf);
        tCrit = get_tcrit(dfMin, cfg.alphaCluster, cfg.tail);

        [maskPos, maskNeg] = threshold_t_map(tObs, tCrit, cfg.tail);

        % ---------- observed clusters ----------
        clPos = find_clusters(maskPos, A, cfg.minCluster);
        clNeg = find_clusters(maskNeg, A, cfg.minCluster);

        clusters = struct('idx',{},'labels',{},'sign',{},'mass',{},'size',{},'pCluster',{});
        kCount = 0;

        for k = 1:numel(clPos)
            idx = clPos{k};
            kCount = kCount + 1;
            clusters(kCount).idx      = idx(:)';
            clusters(kCount).labels   = labels_all(idx(:))';
            clusters(kCount).sign     = +1;
            clusters(kCount).mass     = sum(tObs(idx), 'omitnan');
            clusters(kCount).size     = numel(idx);
            clusters(kCount).pCluster = NaN;
        end

        for k = 1:numel(clNeg)
            idx = clNeg{k};
            kCount = kCount + 1;
            clusters(kCount).idx      = idx(:)';
            clusters(kCount).labels   = labels_all(idx(:))';
            clusters(kCount).sign     = -1;
            clusters(kCount).mass     = abs(sum(tObs(idx), 'omitnan'));
            clusters(kCount).size     = numel(idx);
            clusters(kCount).pCluster = NaN;
        end

        % ---------- permutation null ----------
        nS = size(Y,1);
        permMaxMass = zeros(cfg.nPerm,1);

        for perm = 1:cfg.nPerm
            Gp = G(randperm(nS));
            idx1p = (Gp == g1);
            idx2p = (Gp == g2);

            [tPerm, ~, ~] = welch_t_map(Y, idx2p, idx1p);
            tPerm = tPerm(:);

            [mPos, mNeg] = threshold_t_map(tPerm, tCrit, cfg.tail);

            clp = find_clusters(mPos, A, cfg.minCluster);
            cln = find_clusters(mNeg, A, cfg.minCluster);

            maxMass = 0;

            for k = 1:numel(clp)
                maxMass = max(maxMass, sum(tPerm(clp{k}), 'omitnan'));
            end
            for k = 1:numel(cln)
                maxMass = max(maxMass, abs(sum(tPerm(cln{k}), 'omitnan')));
            end

            permMaxMass(perm) = maxMass;
        end

        % ---------- cluster p-values ----------
        for k = 1:numel(clusters)
            clusters(k).pCluster = (sum(permMaxMass >= clusters(k).mass) + 1) / (cfg.nPerm + 1);
        end

        % ---------- significant cluster mask ----------
        sigMask = false(nCh,1);
        for k = 1:numel(clusters)
            if clusters(k).pCluster < cfg.alphaCorrected
                sigMask(clusters(k).idx) = true;
            end
        end

        % ---------- store ----------
        RES(kInt).tObs            = tObs;
        RES(kInt).dfObs           = dfObs;
        RES(kInt).pObs_unc        = pObs;
        RES(kInt).tCrit           = tCrit;
        RES(kInt).supraMaskPos    = maskPos;
        RES(kInt).supraMaskNeg    = maskNeg;
        RES(kInt).clusters        = clusters;
        RES(kInt).sigClusterMask  = sigMask;
        RES(kInt).permMaxMass     = permMaxMass;
        RES(kInt).meanDeltaGroup1 = mean(Y(idx1,:), 1, 'omitnan')';
        RES(kInt).meanDeltaGroup2 = mean(Y(idx2,:), 1, 'omitnan')';
        RES(kInt).sdDeltaGroup1   = std(Y(idx1,:), 0, 1, 'omitnan')';
        RES(kInt).sdDeltaGroup2   = std(Y(idx2,:), 0, 1, 'omitnan')';
        RES(kInt).diffDelta       = RES(kInt).meanDeltaGroup2 - RES(kInt).meanDeltaGroup1;
        RES(kInt).nValidGroup1    = sum(~isnan(Y(idx1,:)), 1)';
        RES(kInt).nValidGroup2    = sum(~isnan(Y(idx2,:)), 1)';
        RES(kInt).labels_all      = labels_all;

        if cfg.verbose
            nClusters = numel(clusters);
            if isempty(clusters)
                nSig = 0;
            else
                nSig = sum([clusters.pCluster] < cfg.alphaCorrected);
            end
            fprintf('\n[%s | %s]\n', char(metricVar), char(interactionLabel));
            fprintf('n(%s)=%d | n(%s)=%d | clusters=%d | sig=%d\n', ...
                g1, sum(idx1), g2, sum(idx2), nClusters, nSig);
        end
    end

    OUT.(matlab.lang.makeValidName(metricVar)) = RES;
end
end

% ========================================================================
% HELPERS
% ========================================================================

function T = normalize_input_table(T)
T.interaction = strtrim(string(T.interaction));
T.chanLabel   = upper(strtrim(string(T.chanLabel)));
T.subjectID   = strtrim(string(T.subjectID));
T.group       = lower(strtrim(string(T.group)));
end

function [Y, G, subjIDs, ok, msg] = build_delta_matrix(T1, T2, interactionLabel, metricVar, labels_all, label2idx, g1, g2)

Y = [];
G = strings(0,1);
subjIDs = strings(0,1);
ok = false;
msg = "";

X1 = T1(T1.interaction == interactionLabel & ismember(T1.chanLabel, labels_all), :);
X2 = T2(T2.interaction == interactionLabel & ismember(T2.chanLabel, labels_all), :);

if isempty(X1) || isempty(X2)
    msg = "empty data in one condition for this interaction";
    return;
end

K1 = X1(:, {'subjectID','chanLabel','group', char(metricVar)});
K2 = X2(:, {'subjectID','chanLabel','group', char(metricVar)});

K1 = renamevars(K1, char(metricVar), 'val1');
K2 = renamevars(K2, char(metricVar), 'val2');

% ----- normalize again just to be safe -----
K1.subjectID = strtrim(string(K1.subjectID));
K2.subjectID = strtrim(string(K2.subjectID));
K1.chanLabel = upper(strtrim(string(K1.chanLabel)));
K2.chanLabel = upper(strtrim(string(K2.chanLabel)));
K1.group     = lower(strtrim(string(K1.group)));
K2.group     = lower(strtrim(string(K2.group)));

% ----- check duplicates -----
assert_unique_subject_channel(K1, 'T1');
assert_unique_subject_channel(K2, 'T2');

% ----- rename group columns BEFORE join -----
K1 = renamevars(K1, 'group', 'group1');
K2 = renamevars(K2, 'group', 'group2');

% ----- join -----
J = innerjoin(K1, K2, 'Keys', {'subjectID','chanLabel'});

if isempty(J)
    msg = "no overlapping subject/channel rows after join";
    return;
end

% ----- check expected variables exist -----
neededVars = {'subjectID','chanLabel','group1','group2','val1','val2'};
missingVars = neededVars(~ismember(neededVars, J.Properties.VariableNames));
if ~isempty(missingVars)
    error('After innerjoin, missing variables: %s', strjoin(missingVars, ', '));
end

% ----- check group consistency across conditions -----
sameGroup = (J.group1 == J.group2);
if ~all(sameGroup)
    bad = J(~sameGroup, {'subjectID','chanLabel','group1','group2'});
    disp(bad(1:min(height(bad),10), :));
    error('Group mismatch between conditions for some subject/channel rows.');
end

J.group = J.group1;
J.delta = J.val2 - J.val1;

% ----- keep only requested groups -----
J = J((J.group == g1) | (J.group == g2), :);

if isempty(J)
    msg = "no rows remain after group filtering";
    return;
end

subjIDs = unique(J.subjectID, 'stable');
nS = numel(subjIDs);
nCh = numel(labels_all);

Y = nan(nS, nCh);
G = strings(nS,1);

for s = 1:nS
    sid = subjIDs(s);
    Js = J(J.subjectID == sid, :);

    subjGroup = unique(Js.group);
    if numel(subjGroup) ~= 1
        error('Subject %s has inconsistent group assignment.', sid);
    end
    G(s) = subjGroup;

    for r = 1:height(Js)
        lab = char(Js.chanLabel(r));
        if isKey(label2idx, lab)
            Y(s, label2idx(lab)) = Js.delta(r);
        end
    end
end

if size(Y,1) < 4
    msg = sprintf('too few subjects after join (n=%d)', size(Y,1));
    return;
end

if ~any(G == g1) || ~any(G == g2)
    msg = sprintf('one group missing after join (%s=%d, %s=%d)', ...
        g1, sum(G == g1), g2, sum(G == g2));
    return;
end

ok = true;
end

function assert_unique_subject_channel(K, tag)
keys = K(:, {'subjectID','chanLabel'});
[G, keyTable] = findgroups(keys);
counts = splitapply(@numel, K.subjectID, G);

if any(counts > 1)
    badRows = keyTable(counts > 1, :);
    disp(badRows(1:min(height(badRows),10), :));
    error('%s contains duplicate subjectID x chanLabel rows.', tag);
end
end

function [t, df, p] = welch_t_map(Y, idxA, idxB)
% Channelwise Welch t-test with NaN handling.
% Returns t = mean(A) - mean(B)

YA = Y(idxA,:);
YB = Y(idxB,:);

nA = sum(~isnan(YA), 1);
nB = sum(~isnan(YB), 1);

mA = mean(YA, 1, 'omitnan');
mB = mean(YB, 1, 'omitnan');

vA = var(YA, 0, 1, 'omitnan');
vB = var(YB, 0, 1, 'omitnan');

se2 = vA ./ nA + vB ./ nB;
t = (mA - mB) ./ sqrt(se2);

df_num = se2 .^ 2;
df_den = ((vA ./ nA).^2) ./ (nA - 1) + ((vB ./ nB).^2) ./ (nB - 1);
df = df_num ./ df_den;

bad = (nA < 2) | (nB < 2) | ~isfinite(t) | ~isfinite(df) | (df <= 0) | (se2 <= 0);
t(bad)  = NaN;
df(bad) = NaN;

p = nan(size(t));
good = ~bad;
p(good) = 2 * tcdf(-abs(t(good)), df(good));
end

function tCrit = get_tcrit(df, alpha, tail)
switch lower(string(tail))
    case "two"
        tCrit = tinv(1 - alpha/2, df);
    case {"pos","neg"}
        tCrit = tinv(1 - alpha, df);
    otherwise
        error('tail must be "two", "pos", or "neg".');
end
end

function [maskPos, maskNeg] = threshold_t_map(t, tCrit, tail)
t = t(:);

switch lower(string(tail))
    case "two"
        maskPos = isfinite(t) & (t >  tCrit);
        maskNeg = isfinite(t) & (t < -tCrit);
    case "pos"
        maskPos = isfinite(t) & (t > tCrit);
        maskNeg = false(size(t));
    case "neg"
        maskPos = false(size(t));
        maskNeg = isfinite(t) & (t < -tCrit);
    otherwise
        error('tail must be "two", "pos", or "neg".');
end
end

function clusters = find_clusters(mask, A, minCluster)
mask = logical(mask(:));
n = numel(mask);
visited = false(n,1);
clusters = {};

for i = 1:n
    if ~mask(i) || visited(i)
        continue;
    end

    stack = i;
    visited(i) = true;
    comp = i;

    while ~isempty(stack)
        v = stack(end);
        stack(end) = [];

        nb = find(A(v,:));
        nb = nb(mask(nb) & ~visited(nb));

        if ~isempty(nb)
            visited(nb) = true;
            stack = [stack; nb(:)]; %#ok<AGROW>
            comp  = [comp; nb(:)];  %#ok<AGROW>
        end
    end

    comp = unique(comp);
    if numel(comp) >= minCluster
        clusters{end+1} = comp; %#ok<AGROW>
    end
end
end