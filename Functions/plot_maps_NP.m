function plot_maps_NP(ClusterPerm, labels_all, chx, chy, locTab, cfg)

if nargin < 6 || isempty(cfg), cfg = struct(); end

if ~isfield(cfg,'outDir'),            cfg.outDir = pwd; end
if ~exist(cfg.outDir,'dir'), mkdir(cfg.outDir); end

if ~isfield(cfg,'conditionTag'), cfg.conditionTag = ""; end
if ~isfield(cfg,'groupNames'),        cfg.groupNames = ["Controls","Patients"]; end
if ~isfield(cfg,'condNames') || isempty(cfg.condNames), cfg.condNames = ["cond1","cond2"]; end
if ~isfield(cfg,'sigClusterAlpha'),   cfg.sigClusterAlpha = 0.05; end

if ~isfield(cfg,'make_contour'),      cfg.make_contour = false; end
if ~isfield(cfg,'plot_channels'),     cfg.plot_channels = false; end
if ~isfield(cfg,'plot_clabels'),      cfg.plot_clabels = false; end
if ~isfield(cfg,'Np'),                cfg.Np = 1000; end

if ~isfield(cfg,'elecMS'),            cfg.elecMS = 24; end
if ~isfield(cfg,'elecLW'),            cfg.elecLW = 1.0; end
if ~isfield(cfg,'roiMS'),             cfg.roiMS = 90; end
if ~isfield(cfg,'roiLW'),             cfg.roiLW = 2.0; end
if ~isfield(cfg,'roiOutlineLW'),      cfg.roiOutlineLW = 2.4; end

if ~isfield(cfg,'colClusterSig'),     cfg.colClusterSig = [0.00 0.60 0.10]; end

if ~isfield(cfg,'metricLabelY'),      cfg.metricLabelY = '\Delta TE'; end
if ~isfield(cfg,'diffLabelY'),        cfg.diffLabelY   = '\Delta difference'; end

labels_all = string(labels_all(:));
nCh = numel(labels_all);

if isfield(cfg,'metricFields') && ~isempty(cfg.metricFields)
    metricFields = string(cfg.metricFields(:));
else
    fn = fieldnames(ClusterPerm);
    keep = false(size(fn));
    for i = 1:numel(fn)
        keep(i) = isstruct(ClusterPerm.(fn{i}));
    end
    metricFields = string(fn(keep));
end

if isempty(metricFields)
    error('ClusterPerm must contain at least one struct metric field.');
end

for f = 1:numel(metricFields)

    metricName = char(metricFields(f));
    R = ClusterPerm.(metricName);
    metricLabel = get_metric_label(metricName, cfg);

    for k = 1:numel(R)

        S = R(k);

needed = {'interaction','meanDeltaGroup1','meanDeltaGroup2', ...
    'diffDelta','clusters'};

        ok = all(cellfun(@(x) isfield(S,x), needed));

        if ~ok
            warning('Skipping %s | entry %d because required fields are missing.', metricName, k);
            continue;
        end

        interaction = string(S.interaction);

        map1 = S.meanDeltaGroup1(:);
        map2 = S.meanDeltaGroup2(:);
        map3 = S.diffDelta(:);

        if numel(map1) ~= nCh || numel(map2) ~= nCh || numel(map3) ~= nCh
            warning('Skipping %s | %s because map length does not match labels_all.', metricName, interaction);
            continue;
        end

        clusters  = S.clusters;

        allVals = [map1; map2; map3];
        allVals = allVals(isfinite(allVals));

        if isempty(allVals)
            clim = [-1 1];
        else
            hi = quantile(abs(allVals), 0.98);
            if hi < 1e-6, hi = 1; end
            clim = [-hi hi];
        end

        fig = figure('Color','w','Position',[80 80 1550 520]);

        tl = tiledlayout(1,3, ...
            'TileSpacing','compact', ...
            'Padding','compact');

        title(tl, build_super_title(metricLabel, interaction, cfg), ...
            'Interpreter','tex', ...
            'FontWeight','bold', ...
            'FontSize',22);

        % ---------------- Group 1 ----------------
        nexttile;

        plot_topography('all', map1, cfg.make_contour, locTab, ...
            cfg.plot_channels, cfg.plot_clabels, cfg.Np, ...
            '', cfg.metricLabelY, clim, 1);

        title(string(cfg.groupNames(1)), ...
            'Interpreter','none', ...
            'FontWeight','bold', ...
            'FontSize',18);

        hold on;
        draw_overlays(chx, chy, clusters, cfg);
        hold off;

        % ---------------- Group 2 ----------------
        nexttile;

        plot_topography('all', map2, cfg.make_contour, locTab, ...
            cfg.plot_channels, cfg.plot_clabels, cfg.Np, ...
            '', cfg.metricLabelY, clim, 1);

        title(string(cfg.groupNames(2)), ...
            'Interpreter','none', ...
            'FontWeight','bold', ...
            'FontSize',18);

        hold on;
        draw_overlays(chx, chy, clusters, cfg);
        hold off;

        % ---------------- Difference ----------------
        nexttile;

        plot_topography('all', map3, cfg.make_contour, locTab, ...
            cfg.plot_channels, cfg.plot_clabels, cfg.Np, ...
            '', cfg.diffLabelY, clim, 1);

        title(sprintf('%s - %s', ...
            string(cfg.groupNames(2)), string(cfg.groupNames(1))), ...
            'Interpreter','none', ...
            'FontWeight','bold', ...
            'FontSize',18);

        hold on;
        draw_overlays(chx, chy, clusters, cfg);
        hold off;

        safeInteraction = regexprep(char(interaction), '[^\w\-]', '_');
        safeTag = regexprep(char(cfg.conditionTag), '[^\w\-]', '_');

        if isempty(safeTag)
            safeTag = 'condition';
        end

        fname = fullfile(cfg.outDir, ...
            sprintf('Groups_%s_%s_%s.png', metricName, safeInteraction, safeTag));

        exportgraphics(fig, fname, 'Resolution', 600);
        close(fig);

    end
end
end

% =======================================================================
function draw_overlays(chx, chy, clusters, cfg)

scatter(chx, chy, cfg.elecMS, 'o', ...
    'MarkerEdgeColor', [0 0 0], ...
    'MarkerFaceColor', 'none', ...
    'LineWidth', cfg.elecLW);

for k = 1:numel(clusters)

    idx = clusters(k).idx(:);
    idx = idx(idx >= 1 & idx <= numel(chx));

    if isempty(idx)
        continue;
    end

    isSig = isfield(clusters(k),'pCluster') && ...
        ~isnan(clusters(k).pCluster) && ...
        clusters(k).pCluster < cfg.sigClusterAlpha;

    if isSig
        C = cfg.colClusterSig;
    else
        continue;
    end

    xr = chx(idx);
    yr = chy(idx);

    if numel(idx) >= 3
        try
            K = convhull(xr, yr);
            plot(xr(K), yr(K), '-', ...
                'Color', C, ...
                'LineWidth', cfg.roiOutlineLW);
        catch
            plot(xr, yr, '-', ...
                'Color', C, ...
                'LineWidth', cfg.roiOutlineLW);
        end
    elseif numel(idx) == 2
        plot(xr, yr, '-', ...
            'Color', C, ...
            'LineWidth', cfg.roiOutlineLW);
    end

    scatter(xr, yr, cfg.roiMS, 'o', ...
        'MarkerEdgeColor', C, ...
        'MarkerFaceColor', 'none', ...
        'LineWidth', cfg.roiLW);
end
end

% =======================================================================
function txt = pretty_condition_label(x)

x = lower(string(x));

switch x
    case {"hea","healthy"}
        txt = "Healthy bias";

    case {"pat"}
        txt = "Patient bias";

    case {"hp"}
        txt = "High pred.";

    case {"lp"}
        txt = "Low pred.";

    otherwise
        txt = x;
end

end

% =======================================================================
function txt = build_super_title(metricLabel, interaction, cfg)

interactionLabel = format_interaction_label(interaction);

c1 = pretty_condition_label(cfg.condNames(1));
c2 = pretty_condition_label(cfg.condNames(2));

txt = sprintf('%s (%s TE) | \\Delta Condition (%s - %s)', ...
    interactionLabel, ...
    string(metricLabel), ...
    c2, ...
    c1);

end

% =======================================================================
function interactionLabel = format_interaction_label(interaction)

interaction = string(interaction);
parts = split(interaction, "_");

if numel(parts) >= 2

    source = pretty_signal_label(parts(1));
    target = pretty_signal_label(parts(2));

    interactionLabel = sprintf('%s\\rightarrow %s', source, target);

else
    interactionLabel = char(interaction);
end

end

% =======================================================================
function label = pretty_signal_label(x)

x = lower(string(x));

switch x
    case "rr"
        label = "RR";
    case "resp"
        label = "Resp";
    case "delta"
        label = "Delta";
    case "theta"
        label = "Theta";
    case "alpha"
        label = "Alpha";
    case "alpha1"
        label = "Alpha1";
    case "alpha2"
        label = "Alpha2";
    case "beta"
        label = "Beta";
    case "beta1"
        label = "Beta1";
    case "beta2"
        label = "Beta2";
    case "gamma"
        label = "Gamma";
    otherwise
        label = char(x);
end

end

% =======================================================================
function metricLabel = get_metric_label(metricName, cfg)

metricName = string(metricName);

if isfield(cfg,'metricLabels') && ~isempty(cfg.metricLabels)

    if isa(cfg.metricLabels, 'containers.Map')

        if isKey(cfg.metricLabels, char(metricName))
            metricLabel = string(cfg.metricLabels(char(metricName)));
            return;
        end

    elseif isstruct(cfg.metricLabels)

        f = matlab.lang.makeValidName(char(metricName));

        if isfield(cfg.metricLabels, f)
            metricLabel = string(cfg.metricLabels.(f));
            return;
        end
    end
end

switch lower(metricName)
    case {"mi_lin","te_lin"}
        metricLabel = "Linear";
    case {"mi_knn","te_knn"}
        metricLabel = "Nonlinear";
    otherwise
        metricLabel = upper(metricName);
end
end