function comparePrePost_ICA(data_before, data_after)
%   Visually compare FieldTrip data before and after ASR/ICA correction.
%
% USAGE:
%   comparePrePost_ICA(data_before, data_after);
%
% INPUTS:
%   data_before - FieldTrip continuous data before correction.
%   data_after  - FieldTrip continuous data after correction.
%
% The function converts both datasets to EEGLAB format and opens
% vis_artifacts to inspect the difference between the original and
% corrected signals.

% Convert FieldTrip data to EEGLAB format
    EEGraw = fieldtrip2eeglab(data_before);
    EEGclean = fieldtrip2eeglab(data_after);

    fprintf('\nControls for Artifact Viewer:\n');
    fprintf('  [n] : display just the new time series\n');
    fprintf('  [o] : display just the old time series\n');
    fprintf('  [b] : display both time series super-imposed\n');
    fprintf('  [d] : display the difference between both time series\n');
    fprintf('  [+] : increase signal scale\n');
    fprintf('  [-] : decrease signal scale\n');
    fprintf('  [*] : expand time range\n');
    fprintf('  [/] : reduce time range\n\n');

    % Launch EEGLAB viewer
    vis_artifacts(EEGclean, EEGraw);
    hFig = gcf;   
    waitfor(hFig);

end
