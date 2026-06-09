function [dataClean, removedChans] = ChannelRejectionManual(data, threshMultiplier)
% ChannelRejectionManual
%   Manual rejection of noisy channels based on variance screening
%   and visual inspection.
%
% PROCEDURE:
%   1. Compute channel-wise variance.
%   2. Flag channels exceeding threshMultiplier × median variance.
%   3. Display a variance plot highlighting flagged channels.
%   4. Open FieldTrip databrowser for visual inspection of flagged channels.
%   5. Ask the user to enter the final list of channels to remove.
%
% USAGE:
%   [dataClean, removedChans] = ChannelRejectionManual(data);
%   [dataClean, removedChans] = ChannelRejectionManual(data, 2.5);
%
% INPUTS:
%   data             - FieldTrip continuous data structure.
%   threshMultiplier - Optional variance threshold multiplier.
%                      Default = 2.5.
%
% OUTPUTS:
%   dataClean        - FieldTrip data after channel removal.

    if nargin < 2 || isempty(threshMultiplier)
        threshMultiplier = 2.5;  % default value
    end

    %% Compute variance & threshold
    sig       = data.trial{1};
    chanVar   = var(sig, [], 2, 'omitnan');
    medianVar = median(chanVar);
    thresh    = threshMultiplier * medianVar;

    %% Sort channels descending
    [sv, idx]      = sort(chanVar, 'descend');
    sortedLabels   = data.label(idx);
    isFlag         = sv > thresh;
    flaggedLabels  = sortedLabels(isFlag);
    nCh            = numel(sv);

    %% No flagged channels? Skip
    if isempty(flaggedLabels)
        fprintf('No bad channels. Skipping channel rejection.\n');
        dataClean    = data;
        removedChans = {};
        return;
    end

    %% 1) Variance scatter
    scr   = get(0,'ScreenSize');
    figW  = scr(3)*0.6; figH = scr(4)*0.8;
    posX  = (scr(3)-figW)/2; posY = (scr(4)-figH)/2;
    hVarFig = figure('Name','Channel Variance','NumberTitle','off', ...
                     'WindowStyle','normal', ...
                     'Position',[posX,posY,figW,figH]);
    scatter(sv, 1:nCh, 36, 'b','filled'); hold on;
      scatter(sv(isFlag), find(isFlag), 64, 'r','filled');
      xline(thresh, 'r--','LineWidth',1.5);
    hold off;
    ylim([0.5, nCh+0.5]);
    set(gca, 'YTick', 1:nCh, 'YTickLabel', sortedLabels, ...
             'YDir','reverse','FontSize',9,'TickLength',[0 0]);
    xlabel('Variance'); ylabel('Channel');
    title(sprintf('Channel variance (red > %.2f× median = %.2g)', threshMultiplier, thresh));
    legend({'All','Flagged','Thresh'}, 'Location','eastoutside');

    drawnow; shg; pause(0.5);

    fprintf('\nFlagged channels (>%g):\n', thresh);
    disp(flaggedLabels);

    %% 2) Time-course browser of flagged channels
    fprintf('\nOpening FT databrowser for flagged channels.\n');
    cfg = [];
    cfg.viewmode      = 'vertical';
    cfg.continuous    = 'yes';
    cfg.artifactalpha = 0.8;
    cfg.blocksize     = 190;
    cfg.channel       = flaggedLabels;
    cfg.ylim          = 'maxmin';

    hTimeFig = figure('Name','Flagged Channels Timecourses','NumberTitle','off', ...
                      'WindowStyle','normal', ...
                      'Position',[posX,posY,figW,figH]);
    cfg.figure = hTimeFig;
    ft_databrowser(cfg, data);

    drawnow; shg; pause(0.5);

    %% 3) Prompt for final removal (with re-prompt loop)
    validInput = false;
    removedChans = {};

    while ~validInput
        promptStr = '\nChannels to remove (comma separated) or "none":\n';
        userInput = input(promptStr, 's');

        if strcmpi(strtrim(userInput), 'none')
            removedChans = {};
            fprintf('No channels removed.\n');
            validInput = true;

        else
            tokens = regexp(userInput, '\s*,\s*', 'split');
            invalid = setdiff(tokens, data.label);
            if isempty(invalid)
                removedChans = tokens;
                fprintf('Removed channels: %s\n', strjoin(removedChans, ', '));
                validInput = true;
            else
                fprintf('\nInvalid channel(s): %s\n', strjoin(invalid, ', '));
                fprintf('Please re-enter valid channel names from the data.\n');
            end
        end
    end


    %% 4) Apply drop
    if isempty(removedChans)
        dataClean = data;
    else
        keepChans  = setdiff(data.label, removedChans, 'stable');
        cfg2       = [];
        cfg2.channel = keepChans;
        dataClean  = ft_selectdata(cfg2, data);
    end

    %% 5) Cleanup
    if ishandle(hVarFig),  close(hVarFig);  end
    close all
end
