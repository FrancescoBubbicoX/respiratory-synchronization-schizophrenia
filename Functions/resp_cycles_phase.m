function [phase, T, P, half_cycles, half_cycles_sec] = resp_cycles_phase(resp, fs, params, plotDir, participant_id)
% Detect respiratory peaks and troughs, build respiratory half-cycles
% (inhale/exhale), and compute respiratory phase in [0, 2π).
% Optional QC tools allow visual inspection and manual correction of
% unusually long cycles.
%
% Written by: Francesco Bubbico
% Last updated: June 2026
%
% Inputs
%   resp           : respiration signal (row or column vector)
%   fs             : sampling rate (Hz)
%   params         : optional struct with fields:
%       .frac_iqr          : prominence threshold expressed as a fraction
%                            of the signal interquartile range (IQR)
%       .min_dist_pk_s     : minimum distance between detected peaks (s)
%       .min_dist_tr_s     : minimum distance between detected troughs (s)
%       .min_width_s       : minimum width required for detected landmarks (s)
%       .snap_win_s        : local search window used to refine landmark
%                            positions to nearby extrema (s)
%       .guard_s           : exclusion zone around troughs when searching
%                            for peaks within a respiratory cycle (s)
%       .interactive_split : enables manual correction of unusually long
%                            respiratory cycles
%       .make_plots        : generates static QC plots of landmarks,
%                            half-cycles, and phase
%       .browser_plot      : generates an interactive browser for visual
%                            inspection of respiration and phase
%   plotDir        : optional folder for saving static QC plots
%   participant_id : optional string for plot title/filename
%
% Outputs
%   phase            : 1×L respiratory phase in radians [0, 2π)
%   T                : trough indices (1-based samples)
%   P                : peak indices (1-based samples)
%   half_cycles      : [start_idx, end_idx, type], type: 1=inhale, 2=exhale
%   half_cycles_sec  : same as half_cycles, with boundaries in seconds
%
% Notes
%   * Respiratory phase is estimated by mapping trough→peak to [0, π]
%     and peak→trough to [π, 2π].
%   * Default parameters were used for most recordings. For noisier traces,
%     detection parameters, particularly frac_iqr, may require visual QC.

%% 0) Inputs & defaults ----------------------------------------------------

if nargin < 2 || isempty(fs), error('Provide resp and fs.'); end
resp = resp(:)';                 
L    = numel(resp);
t    = (0:L-1)/fs;               % zero-based time (first sample at 0 s)

% Default landmark-detection parameters
def.frac_iqr          = 0.3;
def.min_dist_pk_s     = 1;
def.min_dist_tr_s     = 1;
def.min_width_s       = 0.5;
def.snap_win_s        = 0.15;
def.guard_s           = 0.02;
def.interactive_split = true;
def.make_plots        = true;
def.browser_plot      = true;

if nargin < 3 || isempty(params), params = struct; end
fns = fieldnames(def);
for k = 1:numel(fns)
    if ~isfield(params, fns{k}) || isempty(params.(fns{k}))
        params.(fns{k}) = def.(fns{k});
    end
end

if nargin < 4, plotDir = []; end
if nargin < 5, participant_id = 'participant'; end

frac_iqr      = params.frac_iqr;
min_dist_pk_s = params.min_dist_pk_s;
min_dist_tr_s = params.min_dist_tr_s;
min_width_s   = params.min_width_s;
snap_win_s    = params.snap_win_s;
guard_s       = params.guard_s;
interactive_split = params.interactive_split;

%% 1) Convert time-based params to samples --------------------------------

min_dist_pk = round(min_dist_pk_s*fs);  
min_dist_tr = round(min_dist_tr_s*fs);   
min_width   = round(min_width_s*fs);     
snap        = max(1, round(snap_win_s*fs)); 
guard       = max(1, round(guard_s*fs));   

% Adaptive prominence threshold based on signal variability (IQR).
finite  = isfinite(resp);
iqr_val = iqr(resp(finite));
if ~isfinite(iqr_val) || iqr_val == 0
    iqr_val = max(std(resp(finite)), 0.1);      
end
min_peak_prom = frac_iqr * iqr_val;             

%% 2) Detect candidate peaks and troughs ----------------------------------
% Peaks are detected on the respiration signal; troughs are detected by
% applying the same procedure to the inverted signal.

[~, pklocs] = findpeaks(resp, ...
    'MinPeakDistance',   min_dist_pk, ...
    'MinPeakProminence', min_peak_prom, ...
    'MinPeakWidth',      min_width);

[~, trlocs] = findpeaks(-resp, ...
    'MinPeakDistance',   min_dist_tr, ...
    'MinPeakProminence', min_peak_prom, ...
    'MinPeakWidth',      min_width);

T = sort(trlocs(:)');                             % trough candidates
T = T([true, diff(T) > min_dist_tr]);             
P = sort(pklocs(:)');                             % peak candidates

% Boundaries (safety)
T = T(T>=1 & T<=L);
P = P(P>=1 & P<=L);

%% 3) Stabilize cycles: snap troughs; one peak per trough-to-trough cycle  -------------
% Refine landmark positions and enforce one peak per trough-to-trough cycle.
% Troughs are snapped to the nearest local minimum, and peaks are selected
% as the strongest local maximum within each consecutive trough interval.

if ~isempty(T)
    % Snap troughs to nearby local minima
    for ii = 1:numel(T)
        a = max(1, T(ii)-snap);
        b = min(L, T(ii)+snap);
        [~, j] = min(resp(a:b));
        T(ii) = a + j - 1;                        
    end
    T = unique(T(:));                            

    % Select one peak per trough-to-trough interval
    TT = numel(T)-1;
    Pk = nan(TT,1);

    for ii = 1:TT
        i1 = T(ii); i2 = T(ii+1);
        a  = min(L, max(1, i1 + guard));          
        b  = min(L, max(1, i2 - guard));          

        if b > a
            cand = P(P > a & P < b);              
            if ~isempty(cand)
                [~, j] = max(resp(cand));       
                pk = cand(j);
            else
                [~, j] = max(resp(a:b));          
                pk = a + j - 1;
            end
        else
            pk = i1 + floor((i2 - i1)/2);        
        end

        % Snap selected peak to nearby local maximum
        aa = max(1, pk - snap);
        bb = min(L, pk + snap);
        [~, j2] = max(resp(aa:bb));
        Pk(ii) = aa + j2 - 1;
    end
    P = unique(Pk(~isnan(Pk)));                   
else
     % If no troughs are detected, preserve peak candidates for later handling.
    T = T(:);
    P = P(:);
end

%% 4) Optional interactive correction of long cycles -----------------------
% Very long trough-to-trough intervals may indicate a missed trough.
% When enabled, this step displays each suspect interval and allows manual
% insertion of a missing trough and, if needed, one additional peak.
%
% This QC step is optional and intended for visual correction of clear
% cycle-detection failures in noisy recordings.

if interactive_split && numel(T) >= 2 && numel(P) == numel(T)-1

    % Ensure landmarks are sorted and unique
    T = sort(unique(T(:)));
    P = sort(unique(P(:)));

    % Safety parameters for manual splitting
    min_sub_s = max(0.8, 0.5*min_dist_tr_s);     
    min_sub   = round(min_sub_s*fs);
    pad       = round(0.20*fs);                  

    % Keep track of skipped intervals
    skipped_pairs = []; 

    % Recompute suspect intervals after each edit
    max_passes = 50;  % hard stop to avoid infinite loops
    pass = 0;
    while pass < max_passes
        pass = pass + 1;

        % Recompute suspects 
        T = sort(unique(T(:)));
        if numel(T) < 2, break; end
        tt_s = diff(T)/fs;

        if numel(tt_s) >= 5 && all(isfinite(tt_s))
            thr_robust = median(tt_s) + 2.5*mad(tt_s,1);   % s
        else
            thr_robust = inf;
        end
        hard_cap_s = 9.0;                       % s
        thr_final  = max(hard_cap_s, thr_robust);
        suspects = find(tt_s > thr_final);

        % Drop suspects you already skipped and that still have identical bounds
        if ~isempty(skipped_pairs)
            keep = true(size(suspects));
            for k = 1:numel(suspects)
                i1 = T(suspects(k)); i2 = T(suspects(k)+1);
                if any(all(skipped_pairs == [i1 i2],2))
                    keep(k) = false;
                end
            end
            suspects = suspects(keep);
        end

        if isempty(suspects)
            break; % nothing to fix
        end

        % Work on the first current suspect
        idx = suspects(1);
        i1  = T(idx);
        i2  = T(idx+1);

        % Find the existing peak inside this interval (strongest if many)
        pk_in = P(P > i1 & P < i2);
        if isempty(pk_in)
            % No peak inside → nothing to split meaningfully; mark as skipped
            skipped_pairs(end+1,:) = [i1 i2]; %#ok<AGROW>
            continue;
        end
        [~, jmax] = max(resp(pk_in));
        pk = pk_in(jmax);

        % Prepare segment view
        aVis = max(1, i1 - pad);
        bVis = min(L, i2 + pad);
        tSeg = t(aVis:bVis);
        rSeg = resp(aVis:bVis);

        fig = figure('Name','Split long cycle','Color','w','Position',[120,120,1100,420]); hold on;
        plot(tSeg, rSeg, 'k', 'LineWidth', 1.0);
        yL = min(rSeg); yU = max(rSeg); if yU==yL, yU = yL + 1; end
        patch([t(i1) t(i2) t(i2) t(i1)], [yL yL yU yU], [0.92 0.96 1.00], ...
            'FaceAlpha',0.25, 'EdgeColor','none');
        plot(t(i1), resp(i1), 'bv', 'MarkerFaceColor','b', 'MarkerSize',6);
        plot(t(i2), resp(i2), 'bv', 'MarkerFaceColor','b', 'MarkerSize',6);
        plot(t(pk), resp(pk), 'r^', 'MarkerFaceColor','r', 'MarkerSize',6);
        grid on; xlabel('Time (s)'); ylabel('Respiration (a.u.)');
        title(sprintf('Subject: %s | Suspect cycle %d (dur = %.2f s). Click a trough (split) or press S to skip.', ...
            strrep(participant_id,'_','\_'), idx, (i2-i1)/fs));
        fprintf(['\nCycle [%d]  [%.2f–%.2f s, dur=%.2f s]\n',...
            'Click a breakpoint (trough) inside the interval,\n',...
            'or press ''S'' / right-click / close to skip.\n'], ...
            idx, t(i1), t(i2), (i2-i1)/fs);

        % Get user click/keypress
        tc = NaN; skipped = false;
        try
            [xClick, ~, button] = ginput(1);
            if isempty(xClick) || any(button == [3 's' 'S'])
                skipped = true;
            else
                tc = min(max(round(xClick*fs) + 1, i1+1), i2-1); % clamp strictly inside
            end
        catch
            skipped = true;
        end

        if skipped
            if isvalid(fig), close(fig); end
            skipped_pairs(end+1,:) = [i1 i2]; %#ok<AGROW>
            continue; % recompute suspects on next while-iteration
        end

        % Snap trough to local MIN
        lo = max(i1+1, tc - snap);
        hi = min(i2-1, tc + snap);

        if hi >= lo
            [~, rel] = min(resp(lo:hi));
            tc = lo + rel - 1;

            % Validate split isn't creating micro-intervals
            if (tc - i1) < min_sub || (i2 - tc) < min_sub
                xline(t(tc), '-', 'Invalid split (too close to boundary)', 'Color',[0.8 0 0], 'LineWidth',1.2);
                drawnow; pause(0.8);
                if isvalid(fig), close(fig); end
                skipped_pairs(end+1,:) = [i1 i2]; %#ok<AGROW>
                continue;
            end

            % Apply split: insert trough
            xline(t(tc), '-', 'Chosen trough', 'Color',[0.3 0.0 0.7], 'LineWidth',1.4);
            drawnow; pause(0.2);
            T = sort(unique([T(:); tc]));

            % Decide which half needs a new peak 
            need_left  = (pk >= tc);   % left half [i1, tc] missing a peak if old peak is on right
            need_right = (pk <= tc);   % right half [tc, i2] missing a peak if old peak is on left

            % Simple guard so user clicks inside the correct half (avoid edges)
            gL = max(guard, 1);
            a1 = max(1, i1 + gL); b1 = min(L, tc - gL);   % left half
            a2 = max(1, tc + gL); b2 = min(L, i2 - gL);   % right half

            % Ask for missing peak (user decides the exact spot)
            if need_left && (b1 > a1)
                xline(t(a1),':','', 'Color',[.7 .7 .7]);
                xline(t(b1),':','', 'Color',[.7 .7 .7]);
                fprintf('  -> Click a PEAK for the LEFT half (%.2f–%.2f s), or press S to skip.\n', t(a1), t(b1));
                try
                    [xPk, ~, btn] = ginput(1);
                catch
                    btn = 'S'; xPk = [];
                end
                if ~isempty(xPk) && ~any(btn==[3 's' 'S'])
                    pk_new = min(max(round(xPk*fs)+1, a1), b1);

                    % Snap manually selected peak to nearby local maximum
                    lo = max(a1, pk_new - snap);
                    hi = min(b1, pk_new + snap);
                    if hi >= lo
                        [~, rel] = max(resp(lo:hi));
                        pk_new = lo + rel - 1;
                    end
                    P = unique([P(:); pk_new]);
                    plot(t(pk_new), resp(pk_new), 'r^', 'MarkerFaceColor','r', 'MarkerSize',6);
                else
                    fprintf('    (left peak skipped)\n');
                end
            end

            if need_right && (b2 > a2)
                xline(t(a2),':','', 'Color',[.7 .7 .7]);
                xline(t(b2),':','', 'Color',[.7 .7 .7]);
                fprintf('  -> Click a PEAK for the RIGHT half (%.2f–%.2f s), or press S to skip.\n', t(a2), t(b2));
                try
                    [xPk, ~, btn] = ginput(1);
                catch
                    btn = 'S'; xPk = [];
                end
                if ~isempty(xPk) && ~any(btn == [3 's' 'S'])
                    pk_new = min(max(round(xPk*fs)+1, a2), b2);
                    
                    % Snap manually selected peak to nearby local maximum
                    lo = max(a2, pk_new - snap);
                    hi = min(b2, pk_new + snap);
                    if hi >= lo
                        [~, rel] = max(resp(lo:hi));
                        pk_new = lo + rel - 1;
                    end
                    P = unique([P(:); pk_new]);
                    plot(t(pk_new), resp(pk_new), 'r^', 'MarkerFaceColor','r', 'MarkerSize',6);
                else
                    fprintf('    (right peak skipped)\n');
                end
            end

            % Wait for user inspection before moving to the next suspect interval
            if isvalid(fig), waitfor(fig); end
        end
    end

    % Final sorting and duplicate removal
    T = sort(unique(T(:)));
    P = sort(unique(P(:)));
end


%% 5) Edge handling --------------------------------------------------------
% Ensure that the beginning and end of the signal are covered by alternating
% respiratory landmarks. If the first/last detected landmark is far from the
% signal boundary, the function searches for one missing opposite landmark
% within the edge segment. If no reliable landmark is found, a boundary
% landmark is inserted.
%
% This avoids undefined phase values at the start or end of the recording.

short_s = 0.50;                          % edge proximity threshold
short   = max(1, round(short_s*fs));     

if isempty(T), T = unique([1; L]); end   % fallback when no troughs are detected

% ---------- START edge: iterate until covered ----------
while true
    % Find earliest landmark and gap from signal start
    firstP = inf;    if ~isempty(P), firstP = P(1); end
    firstT = T(1);
    first_mark = min(firstT, firstP);           
    gap_start  = first_mark - 1;                 

    if gap_start <= 0
        break;                                   
    end

    if gap_start <= short
        % If the first landmark is close to the boundary, insert the opposite
        % landmark at sample 1.
        if first_mark == firstT         
            P = unique([1; P(:)]);
            T(T == 1) = [];             
        else                            
            T = unique([1; T(:)]);
            P(P == 1) = [];
        end
        break;
    else
        % Otherwise, search for the opposite landmark within the start segment.
        a = 1 + guard; 
        b = first_mark - guard;
        if b <= a
            if first_mark == firstT, P = unique([1; P(:)]); T(T == 1) = [];
            else,                   T = unique([1; T(:)]); P(P == 1) = [];
            end
            break;
        end

        if first_mark == firstT
            [~, j] = max(resp(a:b)); opp = a + j - 1;   
            if (opp - 1) <= short
                % If the candidate is still close to the boundary, use the boundary instead.
                P = unique([1; P(:)]); T(T == 1) = [];
                break;
            else
                % If the candidate is far enough from the boundary, insert it and repeat.
                P = unique([P(:); opp]);
            end
        else
            [~, j] = min(resp(a:b)); opp = a + j - 1;   
            if (opp - 1) <= short
                T = unique([1; T(:)]); P(P == 1) = [];
                break;
            else
                T = unique([T(:); opp]);
            end
        end
    end
end

% ---------- END edge: iterate until covered ----------
while true
    % Find latest landmark and gap from signal end
    lastP = -inf;   if ~isempty(P), lastP = P(end); end
    lastT = T(end);
    last_mark = max(lastT, lastP);               
    gap_end   = L - last_mark;                   

    if gap_end <= 0
        break;                                  
    end

    if gap_end <= short
        % If the last landmark is close to the boundary, insert the opposite
        % landmark at the final sample.
        if last_mark == lastT           
            P = unique([P(:); L]);
            T(T == L) = [];             
        else                           
            T = unique([T(:); L]);
            P(P == L) = [];
        end
        break;
    else
        % Otherwise, search for the opposite landmark within the end segment.
        a = last_mark + guard;
        b = L - guard;
        if b <= a
            if last_mark == lastT, P = unique([P(:); L]); T(T == L) = [];
            else,                 T = unique([T(:); L]); P(P == L) = [];
            end
            break;
        end

        if last_mark == lastT
            [~, j] = max(resp(a:b)); opp = a + j - 1;
            if (L - opp) <= short
                % Candidate still too close to the boundary: use a boundary landmark.
                P = unique([P(:); L]); T(T == L) = [];
                break;
            else
                % Candidate is sufficiently far from the boundary: insert and repeat.
                P = unique([P(:); opp]);
            end
        else
            [~, j] = min(resp(a:b)); opp = a + j - 1;
            if (L - opp) <= short
                T = unique([T(:); L]); P(P == L) = [];
                break;
            else
                T = unique([T(:); opp]);
            end
        end
    end
end

% Final sorting and duplicate removal
T = unique(sort(T(:)));
P = unique(sort(P(:)));


%% 6) Build respiratory half-cycles ---------------------------------------
% Construct inhalation (trough→peak) and exhalation (peak→trough) segments
% from the ordered sequence of respiratory landmarks.

half_cycles = [];  % [start_idx, end_idx, type], type: 1=inhale, 2=exhale

if ~isempty(T) || ~isempty(P)
    marks  = [T(:);           P(:)];
    types  = [ones(numel(T),1); 2*ones(numel(P),1)];   % 1=trough, 2=peak
    [marks, ord] = sort(marks);
    types = types(ord);

    % Build half-cycles from consecutive landmark pairs
    for i = 1:numel(marks)-1
        a = marks(i);
        b = marks(i+1);
        if b <= a, continue; end                    % skip invalid/non-increasing

        if     (types(i)==1 && types(i+1)==2)       % trough -> peak
            half_cycles(end+1,:) = [a, b, 1];       % inhale
        elseif (types(i)==2 && types(i+1)==1)       % peak -> trough
            half_cycles(end+1,:) = [a, b, 2];       % exhale
        else
            % same-type adjacency (peak->peak or trough->trough): ignore
        end
    end
end

% Convert half-cycle boundaries from samples to seconds
half_cycles_sec = half_cycles;
if ~isempty(half_cycles_sec)
    half_cycles_sec(:,1:2) = (half_cycles_sec(:,1:2) - 1) / fs;
end

%% 7) Respiratory phase estimation ----------------------------------------
% Assign respiratory phase to each half-cycle using linear interpolation.
% Inhalation (trough→peak) is mapped to [0, π], whereas exhalation
% (peak→trough) is mapped to [π, 2π].

phase = nan(1,L);
for ii = 1:size(half_cycles,1)
    i1 = half_cycles(ii,1);
    i2 = half_cycles(ii,2);
    if i2 <= i1+1, continue; end

    if half_cycles(ii,3) == 1                 % inhale
        xi = [i1 i2]; yi = [0 pi];
    else                                       % exhale
        xi = [i1 i2]; yi = [pi 2*pi];
    end

    xx = i1:i2;
    phase(i1:i2) = interp1(xi, yi, xx, 'linear'); 
end

% Interpolate any remaining undefined samples and wrap phase to [0, 2π)
if any(isnan(phase))
    good = find(~isnan(phase)); 
    bad  = find(isnan(phase));
    if ~isempty(good)
        phase(bad) = interp1(good, phase(good), bad, 'linear', 'extrap');
    end
end
phase = mod(phase, 2*pi);

%% 8) Quality-control plots -----------------------------------------------
% Generate optional figures to inspect respiratory landmarks, half-cycles,
% and phase estimation.

if params.make_plots

    % Static overview plot
    W = 2400; H = 1000;                       
    fig = figure('Color','w','Units','pixels','Position',[50 50 W H]);
    try, set(fig,'WindowState','normal'); end 
    hold on;

    yL = min(resp); yU = max(resp);
    if ~isfinite(yL) || ~isfinite(yU) || yL==yU, yL=-1; yU=1; end
    pad = 0.05*(yU - yL + eps); yL = yL - pad; yU = yU + pad;

    for ii = 1:size(half_cycles,1)
        a = (half_cycles(ii,1)-1)/fs; 
        b = (half_cycles(ii,2)-1)/fs;
        c = (half_cycles(ii,3)==1) * [0.7 0.9 0.7] + (half_cycles(ii,3)==2) * [0.9 0.7 0.7];
        patch([a b b a],[yL yL yU yU], c, 'FaceAlpha',0.35,'EdgeColor','none');
    end

    plot(t, resp, 'k');
    if ~isempty(T), plot((T-1)/fs, resp(T), 'bv','MarkerFaceColor','b','MarkerSize',5); end
    if ~isempty(P), plot((P-1)/fs, resp(P), 'r^','MarkerFaceColor','r','MarkerSize',5); end
    xlabel('Time (s)'); ylabel('Z-scored Respiration');
    title(['Respiration with Peaks, Troughs, and Half-Cycles - ' strrep(participant_id,'_','\_')]);
    ylim([yL yU]); grid on;

    ax = gca;
    set(ax,'Units','normalized','Position',[0.06 0.10 0.92 0.82]); 
    set(ax,'LooseInset', max(get(ax,'TightInset'), 0.02*[1 1 1 1]));

    % Save static QC plot
    if ~isempty(plotDir)
        filename = fullfile(plotDir, [participant_id '_respiratory_cycles.png']);
        drawnow;                                
        dpi = 300;                               
        try
            exportgraphics(fig, filename, 'Resolution', dpi, 'BackgroundColor','white');
        catch
            set(fig,'PaperPositionMode','auto'); 
            print(fig, filename, '-dpng', ['-r' num2str(dpi)]);
        end
    end
end

if isfield(params,'browser_plot') && params.browser_plot
    % Interactive browser for inspecting respiration and phase over time.
    windowDur = 30;  plotDecim = 2;
    time_d    = t(1:plotDecim:end);
    resp_d    = resp(1:plotDecim:end);
    phase_d   = phase(1:plotDecim:end);

    fig2 = figure('Name','Respiration Browser','Position',[100,100,1200,600]);
    col_inhale = [0.8 1 0.8]; 
    col_exhale = [1 0.8 0.8];

    % Respiration panel
    ax1 = subplot(2,1,1,'Parent',fig2); hold on;
    yL = min(resp_d); yU = max(resp_d);
    if ~isfinite(yL) || ~isfinite(yU) || yL==yU, yL=-1; yU=1; end
    pad = 0.05*(yU - yL + eps); yL = yL - pad; yU = yU + pad;

    for ii = 1:size(half_cycles,1)
        x1 = (half_cycles(ii,1)-1)/fs; x2 = (half_cycles(ii,2)-1)/fs;
        c = col_inhale; if half_cycles(ii,3)==2, c = col_exhale; end
        patch([x1 x2 x2 x1],[yL yL yU yU],c,'Parent',ax1,'FaceAlpha',0.2,'EdgeColor','none');
    end
    plot(ax1, time_d, resp_d, 'k');
    if ~isempty(T), plot(ax1,(T-1)/fs, resp(T), 'bv','MarkerFaceColor','b','MarkerSize',5); end
    if ~isempty(P), plot(ax1,(P-1)/fs, resp(P), 'r^','MarkerFaceColor','r','MarkerSize',5); end
    ylabel(ax1,'Respiration'); ylim(ax1,[yL yU]);
    title(ax1,'Respiration and Phase (sliding window)');

    % Phase panel
    ax2 = subplot(2,1,2,'Parent',fig2); hold on;
    for ii = 1:size(half_cycles,1)
        x1 = (half_cycles(ii,1)-1)/fs; x2 = (half_cycles(ii,2)-1)/fs;
        c = col_inhale; if half_cycles(ii,3)==2, c = col_exhale; end
        patch([x1 x2 x2 x1],[0 0 2*pi 2*pi],c,'Parent',ax2,'FaceAlpha',0.2,'EdgeColor','none');
    end
    plot(ax2, time_d, phase_d, 'k','LineWidth',1.2);
    ylabel(ax2,'Phase (rad)'); ylim(ax2,[0 2*pi]);
    yticks(ax2,[0 pi 2*pi]); yticklabels(ax2,{'0','\pi','2\pi'});
    xlabel(ax2,'Time (s)');

    % Link axes and add time-window slider
    linkaxes([ax1,ax2],'x');
    if (time_d(end)-time_d(1)) > windowDur
        set([ax1,ax2],'XLim',[time_d(1), time_d(1)+windowDur]);
        uicontrol('Style','slider','Parent',fig2,...
            'Min',time_d(1),'Max',time_d(end)-windowDur,'Value',time_d(1),...
            'Units','normalized','Position',[0.1,0.02,0.7,0.04],...
            'Callback',@(src,~) set([ax1,ax2],'XLim',...
            [get(src,'Value'), get(src,'Value')+windowDur]));
    end
end

if params.make_plots
% Optional: uncomment to inspect the figure before closing.
% waitfor(fig);
close all

end
