function [respTrim_corr, qc] = resp_start_envelope_correction(respTrim, fs, subjectID, opts)
%
% Inputs
%   respTrim   : respiration segment (trimmed finite), row/column vector
%   fs         : sampling rate (Hz) of respTrim
%   subjectID  : string/char used for titles/filenames
%   opts       : struct (optional) with fields:
%       .start_s        window to evaluate "start" MAD
%       .thr_hi         correct if MAD_start/MAD_ref > thr_hi
%       .thr_lo        
%       .envWin_s      envelope smoothing window (s)
%       .maxFix_s      max duration that can be corrected (s)
%       .ramp_s        fade-out duration (s)
%       .ref_t1        reference window start time (s)
%       .ref_t2        reference window end time (s)
%       .cap_lo        min scale cap
%       .cap_hi        max scale cap
%       .tol_lo        stability band lower bound on env ratio
%       .tol_hi        stability band upper bound on env ratio
%       .req_s         required consecutive stable seconds
%       .make_plots    show QC plot
%       .save_plot     save QC plot to disk
%       .plotDir       folder to save plot (required if save_plot=true)
%       .qc_zoom_s     zoom duration for panel 4
%
% Outputs
%   respTrim_corr : corrected respiration (same size as input, row vector)
%   qc           : struct with diagnostics (gating, refEnv, scale, etc.)

if nargin < 4 || isempty(opts), opts = struct(); end
if nargin < 3 || isempty(subjectID), subjectID = "subject"; end

% ---- defaults ----
d.start_s    = 4;
d.thr_hi     = 3;
d.envWin_s   = 5;
d.maxFix_s   = 20;
d.ramp_s     = 2;
d.ref_t1     = 20;
d.ref_t2     = 160;
d.cap_lo     = 0.33;
d.cap_hi     = 1.66;
d.tol_lo     = 0.80;
d.tol_hi     = 1.20;
d.req_s      = 5;
d.make_plots = true;
d.save_plot  = false;
d.plotDir    = '';
d.qc_zoom_s  = 60;

fn = fieldnames(d);
for k = 1:numel(fn)
    if ~isfield(opts, fn{k}) || isempty(opts.(fn{k}))
        opts.(fn{k}) = d.(fn{k});
    end
end
if ~isfield(opts,'thr_lo') || isempty(opts.thr_lo)
    opts.thr_lo = 1/opts.thr_hi;
end

% ---- enforce vector shape (row) ----
x = respTrim(:).';
x0 = x;

N = numel(x);
t = (0:N-1)/fs;

% ---- reference mask: 60–150 s (fallbacks) ----
refMask = (t >= opts.ref_t1) & (t <= opts.ref_t2);
if ~any(refMask)
    warning('No samples in reference window %.0f–%.0f s. Using t>=%.0f s.', ...
        opts.ref_t1, opts.ref_t2, opts.ref_t1);
    refMask = (t >= opts.ref_t1);
end
if nnz(refMask) < round(5*fs)
    warning('Reference window too short. Using entire signal as reference.');
    refMask = true(size(t));
end

% ---- gating by robust amplitude (MAD) ----
Nstart = min(N, max(1, round(opts.start_s*fs)));
amp_start = mad(x(1:Nstart), 1);
amp_ref   = mad(x(refMask), 1);
ratio_mad = amp_start / max(amp_ref, eps);

doCorr = (ratio_mad > opts.thr_hi) || (ratio_mad < opts.thr_lo);

% ---- envelope (always compute if plotting OR correction) ----
env   = abs(hilbert(x));
env_s = movmedian(env, max(3, round(opts.envWin_s*fs)), 'omitnan');

refEnv = median(env_s(refMask), 'omitnan');
refEnv = max(refEnv, eps);

% ---- default outputs ----
respTrim_corr = x;
scale_blend   = ones(size(x));
scale_raw     = ones(size(x));
Nfix_s_eff    = NaN;
iStable       = NaN;

% ---- apply correction: start-only, adaptive stop, ramp ----
if doCorr
    % scale to normalize envelope to ref
    scale_raw = refEnv ./ env_s;
    scale_raw(~isfinite(scale_raw)) = 1;

    % cap
    scale_raw = min(max(scale_raw, opts.cap_lo), opts.cap_hi);

    % detect stabilization time (within band for req_s seconds)
    ratio_env = env_s / refEnv;
    reqN = max(1, round(opts.req_s*fs));
    stableMask = (ratio_env >= opts.tol_lo) & (ratio_env <= opts.tol_hi);
    stableRun  = conv(double(stableMask), ones(1,reqN), 'same') >= reqN;

    maxFixN = min(N, max(1, round(opts.maxFix_s*fs)));
    iStable = find(stableRun(1:maxFixN), 1, 'first');

    if isempty(iStable)
        iStable = min(N, round(min(opts.start_s, opts.maxFix_s)*fs));
    end

    Nfix_s_eff = (iStable-1)/fs;
    Nfix  = iStable;
    Nramp = min(N - Nfix, max(0, round(opts.ramp_s*fs)));

    alpha = zeros(size(x));
    alpha(1:Nfix) = 1;
    if Nramp > 0
        alpha(Nfix+(1:Nramp)) = linspace(1,0,Nramp);
    end

    scale_blend = 1 + alpha .* (scale_raw - 1);
    respTrim_corr = x .* scale_blend;
end

% ---- QC struct ----
qc = struct();
qc.subjectID    = string(subjectID);
qc.doCorr       = doCorr;
qc.ratio_mad    = ratio_mad;
qc.amp_start    = amp_start;
qc.amp_ref      = amp_ref;
qc.ref_t1       = opts.ref_t1;
qc.ref_t2       = opts.ref_t2;
qc.refEnv       = refEnv;
qc.Nfix_s       = Nfix_s_eff;
qc.ramp_s       = opts.ramp_s;
qc.cap_lo       = opts.cap_lo;
qc.cap_hi       = opts.cap_hi;
qc.t            = t;
qc.env_s        = env_s;
qc.scale_raw    = scale_raw;
qc.scale_blend  = scale_blend;
qc.iStable      = iStable;

% ---- console log ----
if doCorr
    fprintf('Resp corr | %s | MAD ratio(start/ref)=%.3f | correction applied | Nfix≈%.1fs\n', ...
        qc.subjectID, ratio_mad, Nfix_s_eff);
else
    fprintf('Resp corr | %s | MAD ratio(start/ref)=%.3f | no correction needed\n', ...
        qc.subjectID, ratio_mad);
end

% ---- QC plot ----
if opts.make_plots && doCorr

    fig = figure('Color','w','Position',[100 100 1200 800]);

    % Panel 1: original
    subplot(3,2,1);
    plot(t, x0, 'k');
    title('Original respTrim'); xlabel('Time (s)'); ylabel('Amp (a.u.)'); grid on

    % Panel 2: after
    subplot(3,2,2);
    plot(t, respTrim_corr, 'r');
    if doCorr
        title(sprintf('Corrected (start-only) | Nfix≈%.1fs + ramp %ds', Nfix_s_eff, opts.ramp_s));
    else
        title('No correction applied (passed gating)');
    end
    xlabel('Time (s)'); ylabel('Amp (a.u.)'); grid on

    % Panel 3: overlay fixed scale
    subplot(3,2,3);
    yl = prctile(x0,[1 99]);
    plot(t, x0,'k'); hold on; plot(t, respTrim_corr,'r'); ylim(yl);
    title('Overlay (fixed y-limits)'); legend('Original','After');
    xlabel('Time (s)'); ylabel('Amp (a.u.)'); grid on

    % Panel 4: zoom
    subplot(3,2,4);
    idxZ = t <= min(opts.qc_zoom_s, t(end));
    plot(t(idxZ), x0(idxZ),'k'); hold on; plot(t(idxZ), respTrim_corr(idxZ),'r');
    if doCorr && ~isnan(Nfix_s_eff)
        xline(Nfix_s_eff,'k--','End correction');
        xline(Nfix_s_eff+opts.ramp_s,'k:','End ramp');
    end
    title(sprintf('First %.0f s (zoom)', min(opts.qc_zoom_s, t(end))));
    legend('Original','After'); xlabel('Time (s)'); ylabel('Amp (a.u.)'); grid on

    % Panel 5: envelope + reference window
    subplot(3,2,5);
    plot(t, env_s, 'b'); hold on; yline(refEnv,'r--','refEnv');
    xline(opts.ref_t1,'k:'); xline(opts.ref_t2,'k:');
    if doCorr && ~isnan(Nfix_s_eff), xline(Nfix_s_eff,'k--'); end
    title(sprintf('Smoothed envelope (%.0fs) | ref %.0f–%.0fs', opts.envWin_s, opts.ref_t1, opts.ref_t2));
    xlabel('Time (s)'); ylabel('Envelope'); grid on

    % Panel 6: applied scale
    subplot(3,2,6);
    plot(t, scale_blend, 'm'); hold on; yline(1,'k--');
    if doCorr && ~isnan(Nfix_s_eff)
        xline(Nfix_s_eff,'k--'); xline(Nfix_s_eff+opts.ramp_s,'k:');
    end
    title(sprintf('Applied scale (capped %.2f–%.2f)', opts.cap_lo, opts.cap_hi));
    xlabel('Time (s)'); ylabel('Multiplier'); grid on

    sgtitle(sprintf('Resp amplitude correction QC | %s | MAD ratio=%.2f', qc.subjectID, ratio_mad));

    if opts.save_plot
        if isempty(opts.plotDir)
            error('opts.plotDir must be set when opts.save_plot = true');
        end
        if ~exist(opts.plotDir,'dir'), mkdir(opts.plotDir); end
        fname = fullfile(opts.plotDir, sprintf('%s_resp_ampcorr_qc.png', regexprep(char(qc.subjectID),'[^\w-]','_')));
        exportgraphics(fig, fname, 'Resolution', 300);
        close(fig);
        qc.qc_plot_file = fname;
    else
        qc.qc_plot_file = "";
    end
else
    qc.qc_plot_file = "";
end
end