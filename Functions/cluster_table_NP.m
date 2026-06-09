function T_clu = cluster_table_NP(ClusterPerm, T1, T2, labels_all, metricVars, sigAlpha, varargin)
% Build a long cluster-level table:
% one row per subject x condition x cluster x metricVar
%
% Each cluster value is the mean metric across the electrodes
% belonging to that cluster.

labels_all = upper(strtrim(string(labels_all(:))));
metricVars = string(metricVars(:));
sigAlpha   = double(sigAlpha);

% Stack conditions
T = [T1; T2];

% Normalize
T.subjectID   = strtrim(string(T.subjectID));
T.group       = lower(strtrim(string(T.group)));
T.condition   = lower(strtrim(string(T.selectedCond)));
T.interaction = strtrim(string(T.interaction));
T.chanLabel   = upper(strtrim(string(T.chanLabel)));

rows = table();

for mv = 1:numel(metricVars)

    metricVar = metricVars(mv);
    permField = matlab.lang.makeValidName(metricVar);

    if ~isfield(ClusterPerm, permField)
        warning('ClusterPerm has no field %s. Skipping.', permField);
        continue;
    end

    R = ClusterPerm.(permField);

    for k = 1:numel(R)

        if ~isfield(R(k), 'interaction') || isempty(R(k).interaction)
            continue;
        end

        interaction = string(R(k).interaction);

        if interaction == "" || ~isfield(R(k), 'clusters') || isempty(R(k).clusters)
            continue;
        end

        cl = R(k).clusters;

        for cIdx = 1:numel(cl)

            chIdx = cl(cIdx).idx(:);
            chIdx = chIdx(chIdx >= 1 & chIdx <= numel(labels_all));

            if isempty(chIdx)
                continue;
            end

            elec  = labels_all(chIdx);
            nElec = numel(elec);

            X = T(T.interaction == interaction & ismember(T.chanLabel, elec), :);

            if isempty(X)
                continue;
            end

            if ~ismember(metricVar, string(X.Properties.VariableNames))
                continue;
            end

            % Average within cluster per subject x condition x group
            [gid, sid, cond, grp] = findgroups( ...
                X.subjectID, X.selectedCond, X.group);

            metricMean   = splitapply(@(x) mean(x,'omitnan'),   X.(metricVar), gid);

            Tc = table();

            Tc.subjectID    = sid;
            Tc.group        = grp;
            Tc.condition    = cond;
            Tc.interaction  = repmat(interaction, size(metricMean));
            Tc.metricVar    = repmat(metricVar, size(metricMean));
            Tc.clusterID    = repmat(metricVar + "_" + interaction + "_clu_" + string(cIdx), size(metricMean));
            Tc.clusterP     = repmat(cl(cIdx).pCluster, size(metricMean));
            Tc.isSigCluster = repmat(cl(cIdx).pCluster < sigAlpha, size(metricMean));

            if isfield(cl(cIdx), 'sign')
                Tc.clusterSign = repmat(cl(cIdx).sign, size(metricMean));
            else
                Tc.clusterSign = nan(size(metricMean));
            end

            if isfield(cl(cIdx), 'mass')
                Tc.clusterMass = repmat(cl(cIdx).mass, size(metricMean));
            else
                Tc.clusterMass = nan(size(metricMean));
            end

            Tc.nElectrodes   = repmat(nElec, size(metricMean));
            Tc.metric_mean   = metricMean;
            Tc.electrodes    = repmat(strjoin(elec, '/'), size(metricMean));

            rows = [rows; Tc]; %#ok<AGROW>
        end
    end
end

T_clu = rows;
end