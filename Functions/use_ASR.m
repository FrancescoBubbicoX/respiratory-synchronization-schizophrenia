function data_clean = use_ASR(data, burstCriterion)
% use_ASR  Perform Artifact Subspace Reconstruction (ASR) on continuous FieldTrip data
%
% USAGE
%   data_clean = use_ASR(data);            
%   data_clean = use_ASR(data, 30);        % set BurstCriterion = 30 SD
%
% INPUT
%   data           - FieldTrip continuous struct with:
%                      .trial{1} = [nCh × nSamp]
%                      .time{1}  = [1 × nSamp]
%   burstCriterion - (optional) SD threshold for ASR bursts (default = 20)
%
% OUTPUT
%   data_clean     - same FieldTrip struct, with .trial cleaned by ASR
%                    and .time restored to original

    % Default threshold
    if nargin<2 || isempty(burstCriterion)
        burstCriterion = 20;
    end

    % Save original time vector
    origTime = data.time;

    % Convert FieldTrip → EEGLAB
    EEG = fieldtrip2eeglab(data);

    % Fix the EEGLAB header so it matches your FieldTrip data
    EEG.pnts   = size(EEG.data,2);          % number of samples per trial

    % these should be equal, if not fix.
    % EEG.nbchan
    % size(EEG.data,1)

    % these should be equal, if not fix.
    % EEG.srate
    % data.fsample

    % these should be equal, if not fix.
    % EEG.times
    % origTime{1};               % full time vector

    % Keep a copy of the raw data for visualization
    rawEEG = EEG;

    % Run ASR: only remove bursts above threshold
    EEG = pop_clean_rawdata( EEG, ...
        'FlatlineCriterion',   'off', ...
        'ChannelCriterion',    'off', ...
        'LineNoiseCriterion',  'off', ...
        'Highpass',            'off', ...
        'BurstCriterion',      burstCriterion, ...
        'WindowCriterion',     'off', ...    
        'BurstRejection',      'off', ...   
        'Distance',            'Euclidian' ...
    );

    fprintf('\nControls for Artifact Viewer:\n');
    fprintf('  [n] : display just the new time series\n');
    fprintf('  [o] : display just the old time series\n');
    fprintf('  [b] : display both time series super-imposed\n');
    fprintf('  [d] : display the difference between both time series\n');
    fprintf('  [+] : increase signal scale\n');
    fprintf('  [-] : decrease signal scale\n');
    fprintf('  [*] : expand time range\n');
    fprintf('  [/] : reduce time range\n\n');
        
    % Check if anything changed
    diffMat = EEG.data - rawEEG.data;
    if max(abs(diffMat(:))) < 1e-12
        fprintf('ASR made no changes (no bursts exceeded %g SD.)\n Skipping visualization.\n', burstCriterion);
    else
        % 8) Visualize what ASR removed, and wait for the user to close it
        vis_artifacts(EEG, rawEEG);      % plot before/after
        hFig = gcf;    % get handle to the artifact‐comparison window
        waitfor(hFig);
    end

    % Convert cleaned data back to FieldTrip with 'raw' to preserve structure
    FT = eeglab2fieldtrip(EEG, 'raw', 'none');

    % 1Rebuild final FieldTrip struct, restoring original time
    data_clean       = data;        % copy all other fields
    data_clean.trial = FT.trial;    % cleaned data
    data_clean.time  = origTime;    % original time vector

end


