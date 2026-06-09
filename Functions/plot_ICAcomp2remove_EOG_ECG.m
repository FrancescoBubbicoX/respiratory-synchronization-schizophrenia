function figHandle = plot_ICAcomp2remove_EOG_ECG(compFlag, tmp)
% PLOTDATABROWSERSTYLE - Plots ICA components + Ocular + ECG in subplots,
% with horizontal (time) scrolling and vertical (amplitude) scaling sliders.
%
% Features:
%   - Ocular always on top (dark green)
%   - ECG always on bottom (red)
%   - All signals z-scored for comparability
%   - Fixed y-limits [-5 5]
%   - Optional downsampling for faster scrolling
%
% Usage:
%   figHandle = plotDatabrowserStyle(compFlag, tmp)
%
% Inputs:
%   compFlag - Struct with:
%       .comp    : FieldTrip ICA comp data (comp.trial{1}, comp.time{1})
%       .dOcular : Ocular data (dOcular.trial{1})
%       .dECG    : ECG data (dECG.trial{1})
%   tmp      - Vector of selected ICA component indices (e.g. [2,6,29])

    %% --- Settings ---
    plotDecim = 5;   % Plot every 10th sample (set =1 for full resolution)
    fixedYLim = [-3 3];  % same scale for all channels

    %% 1) Extract Data
    time     = compFlag.comp.time{1};
    icaData  = compFlag.comp.trial{1}(tmp, :);
    nICA     = size(icaData, 1);
    ocular   = compFlag.dOcular.trial{1};
    ecg      = compFlag.dECG.trial{1}(1, :);

    % Downsample for plotting only
    time     = time(1:plotDecim:end);
    icaData  = icaData(:,1:plotDecim:end);
    ocular   = ocular(1:plotDecim:end);
    ecg      = ecg(1:plotDecim:end);

    %% 2) Reorder: Ocular (top), ICA (middle), ECG (bottom)
    nTotal = nICA + 2;
    channelNames  = cell(nTotal,1);
    channelsData  = cell(nTotal,1);
    channelColors = cell(nTotal,1);

    % Ocular first
    channelNames{1}  = 'Ocular';
    channelsData{1}  = ocular;
    channelColors{1} = [0 0.5 0]; % dark green

    % ICA in the middle
    for i = 1:nICA
        channelNames{i+1}  = sprintf('ICA %d', tmp(i));
        channelsData{i+1}  = icaData(i,:);
        channelColors{i+1} = [0 0 1]; % blue
    end

    % ECG last
    channelNames{end}  = 'ECG';
    channelsData{end}  = ecg;
    channelColors{end} = [1 0 0]; % red

    %% 3) Create Figure & Subplots
    figHandle = figure('Name','Databrowser-Style','Position',[100,100,1200,900]);
    axAll = gobjects(nTotal,1);

    for c = 1:nTotal
        axAll(c) = subplot(nTotal,1,c,'Parent',figHandle);

        % --- Normalize: z-score each channel ---
        x = channelsData{c};
        x = (x - mean(x,'omitnan')) ./ std(x,'omitnan');
        channelsData{c} = x; % overwrite normalized

        % Plot
        plot(axAll(c), time, x, 'LineWidth',1, 'Color', channelColors{c});
        ylabel(axAll(c), channelNames{c}, 'FontWeight','bold');

        % Fixed y-limits for comparability
        ylim(axAll(c), fixedYLim);

        % Clean look (no grid)
        if c < nTotal
            set(axAll(c),'XTickLabel',[]);
        else
            xlabel(axAll(c),'Time (s)');
        end
    end

    % Link x-axes so scrolling is synchronized
    linkaxes(axAll,'x');

    %% 4) Horizontal scrolling
    windowDuration = 10; % seconds per view
    if (time(end)-time(1)) > windowDuration
        set(axAll,'XLim',[time(1), time(1)+windowDuration]);
        uicontrol('Style','slider','Parent',figHandle,...
            'Min',time(1),'Max',time(end)-windowDuration,'Value',time(1),...
            'Units','normalized','Position',[0.1,0.02,0.7,0.03],...
            'Callback',@(src,~) set(axAll,'XLim',...
                [get(src,'Value'), get(src,'Value')+windowDuration]));
    end

    %% 5) Vertical scaling (optional)
    uicontrol('Style','slider','Parent',figHandle,...
        'Min',0.5,'Max',5,'Value',1,...
        'Units','normalized','Position',[0.82,0.02,0.15,0.03],...
        'Callback',@amplitudeCallback);

    function amplitudeCallback(src,~)
        scaleFactor = get(src,'Value');
        for c = 1:nTotal
            ylim(axAll(c), fixedYLim * scaleFactor);
        end
    end
end
