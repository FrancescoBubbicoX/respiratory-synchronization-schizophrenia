function [locs, peaks] = detect_R_peaks(ecg, fs)
% DETECT_R_PEAKS  Robust R-peak detection (Pan–Tompkins style), NaN-safe at edges
%   [locs, peaks] = detect_R_peaks(ecg, fs) returns sample indices locs and
%   amplitudes peaks of R-peaks in the vector ecg, sampled at fs Hz.
%
%   NaN handling:
%     - Supports NaNs ONLY at the beginning/end (common after alignment/padding)
%     - If NaNs are found inside the valid segment, the function errors
%
%   Steps:
%     1) Butterworth band-pass filter (5–15 Hz)
%     2) Differentiate
%     3) Square
%     4) Moving‐window integration (~150 ms)
%     5) Adaptive thresholding + findpeaks
%     6) Refine each R-peak on band-passed ECG within ±50 ms

    % --- force column ---
    ecg = ecg(:);

    % --- NaN trimming (edges only) ---
    validIdx = ~isnan(ecg);
    if ~any(validIdx)
        locs = [];
        peaks = [];
        return;
    end

    firstValid = find(validIdx, 1, 'first');
    lastValid  = find(validIdx, 1, 'last');

    % If NaNs exist inside the valid region, refuse (to avoid silent bad output)
    if any(isnan(ecg(firstValid:lastValid)))
        error('detect_R_peaks: NaNs found inside signal. Only edge-NaNs supported.');
    end

    ecg_trim = ecg(firstValid:lastValid);

    % 1) band-pass 5–15 Hz
    [b,a] = butter(2, [5 15]/(fs/2), 'bandpass');
    ecg_bp = filtfilt(b, a, ecg_trim);

    % 2) differentiate
    diff_ecg = [0; diff(ecg_bp)];

    % 3) square
    sq_ecg = diff_ecg .^ 2;

    % 4) moving‐window integration (~150 ms)
    win = round(0.150 * fs);
    mwi = conv(sq_ecg, ones(win,1)/win, 'same');

    % 5) adaptive threshold + findpeaks
    th = mean(mwi) + 0.5*std(mwi);
    [~, locs] = findpeaks(mwi, ...
                         'MinPeakHeight', th, ...
                         'MinPeakDistance', round(0.3*fs));

    % 6) refine peaks on the filtered ECG
    halfwin = round(0.050 * fs);
    peaks   = zeros(size(locs));
    for i = 1:numel(locs)
        lo = max(1,   locs(i)-halfwin);
        hi = min(numel(ecg_bp), locs(i)+halfwin);
        [pk, idx] = max(ecg_bp(lo:hi));
        peaks(i)  = pk;
        locs(i)   = lo + idx - 1;
    end

    % --- map back to original indexing ---
    locs = locs + (firstValid - 1);
end
