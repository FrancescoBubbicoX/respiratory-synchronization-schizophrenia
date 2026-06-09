function T_out = export_electrode_table_NP(S1, S2, interactions_to_use, labels_all, varargin)

p = inputParser;
p.addParameter('metricVars', {'te_lin','te_knn'}, @(x)iscell(x)||isstring(x));
p.addParameter('addGlobal', true, @(x)islogical(x)&&isscalar(x));
p.parse(varargin{:});
cfg = p.Results;

metricVars          = string(cfg.metricVars(:));
interactions_to_use = string(interactions_to_use(:));
labels_all          = upper(strtrim(string(labels_all(:))));

% Stack conditions
T = [S1; S2];

% Normalize key columns
T.subjectID    = strtrim(string(T.subjectID));
T.selectedCond = lower(strtrim(string(T.selectedCond)));
T.interaction  = strtrim(string(T.interaction));
T.chanLabel    = upper(strtrim(string(T.chanLabel)));

% Optional group column
if ismember("group", string(T.Properties.VariableNames))
    T.group = lower(strtrim(string(T.group)));
else
    T.group = strings(height(T),1);
end

% Keep only requested interactions and scalp electrodes
T = T(ismember(T.interaction, interactions_to_use) & ...
      ismember(T.chanLabel, labels_all), :);

% Build long table
rows = table();

for mv = 1:numel(metricVars)

    v = metricVars(mv);

    if ~ismember(v, string(T.Properties.VariableNames))
        error('Metric variable "%s" not found in input tables.', v);
    end

    Ti = table();

    Ti.subjectID   = T.subjectID;
    Ti.group       = T.group;
    Ti.condition   = T.selectedCond;
    Ti.interaction = T.interaction;
    Ti.electrode   = T.chanLabel;
    Ti.metricVar   = repmat(v, height(T), 1);
    Ti.metricValue = T.(v);

    rows = [rows; Ti]; %#ok<AGROW>
end

% Add GLOBAL rows
if cfg.addGlobal

    [gid, sid, grp, cond, interaction, metricVar] = findgroups( ...
        rows.subjectID, ...
        rows.group, ...
        rows.condition, ...
        rows.interaction, ...
        rows.metricVar);

    metricGlobal = splitapply(@(x) mean(x,'omitnan'), rows.metricValue, gid);

    Tg = table();

    Tg.subjectID   = sid;
    Tg.group       = grp;
    Tg.condition   = cond;
    Tg.interaction = interaction;
    Tg.electrode   = repmat("GLOBAL", size(metricGlobal));
    Tg.metricVar   = metricVar;
    Tg.metricValue = metricGlobal;

    rows = [rows; Tg];
end

% Final output
T_out = rows(:, { ...
    'subjectID', ...
    'group', ...
    'condition', ...
    'interaction', ...
    'electrode', ...
    'metricVar', ...
    'metricValue'});

end