
% GENERATE_HEATMAP_TIMESERIES_HEATMAP1
%  - Prompts for Project folder (contains code, images, output, figure)
%  - Ensures figure folder exists
%  - Loads stage1To2RawData.mat, stage3To4RawData.mat, stage5To6RawData.mat
%    from <Project>/output
%  - Expects variables T1..T8 and Time_line inside each .mat
%  - Removes first 2 rows from each T (metadata rows)
%  - Plots channels 1..4 (T1..T4) as 1x4 subplots using heatmap1()
%  - Uses time vector to place xticks; default tick interval = 4 hours
%  - Saves PNG and FIG to <Project>/figure
%
% REQUIREMENT: heatmap1 must be on the MATLAB path (user-provided).
%
% Usage: put this file in Project/code and run it.

    close all; clc;
    % Adding the directory of important function that will be used here.
    addpath ('E:\rahmans4\EXPERIMENTAL DATA\NIA_Experiment\Scripts');
    % --- Ask user for project folder
    projectFolder = uigetdir(pwd, 'Select Project folder (contains code, images, output, figure)');
    if isequal(projectFolder, 0)
        disp('No folder selected. Exiting.');
        return;
    end

    outputFolder = fullfile(projectFolder, 'output_4Pxl_NucArea_CellPose');
    figureFolder = fullfile(projectFolder, 'figure');

    if ~exist(figureFolder, 'dir')
        mkdir(figureFolder);
    end

    % Conditions
    conds = {
        '5xY', fullfile(outputFolder, 'stage1To2RawData.mat');
        '5xM', fullfile(outputFolder, 'stage3To4RawData.mat');
        'M'   , fullfile(outputFolder, 'stage5To6RawData.mat')
        };

    % Channel Names
    channelNames = { ...
        'Methoxy-x04', ...
        'LyspH Green', ...
        'A\beta-AF555', ...
        'LysoPrime deep red'};

    % Tick interval
    intervalH = input('Enter tick interval in hours (default = 4): ');
    if isempty(intervalH) || ~isnumeric(intervalH) || intervalH <= 0
        intervalH = 4;
    end

    mycolormap = 'jet';

    for ci = 1:size(conds,1)

        condName = conds{ci,1};
        matFile  = conds{ci,2};

        if ~exist(matFile,'file')
            warning('Missing file: %s', matFile);
            continue;
        end

        S = load(matFile);

        if ~isfield(S,'Time_line')
            warning('Time_line missing in %s', matFile);
            continue;
        end

        tvec = double(S.Time_line(:))./60; % minutes → hours

        % Extract T1-T4
        T = cell(1,4);
        for ch = 1:4
            varname = sprintf('T%d',ch);
            if isfield(S,varname)
                temp = double(S.(varname));
                if size(temp,1) > 2
                    temp = temp(3:end,:); % remove metadata rows
                end
            else
                temp = [];
            end
            T{ch} = temp;
        end

        fh = figure('Units','inches','Position',[0 0 15 8],'Color','w');

        for ch = 1:4

            subplot(1,4,ch);
            data = T{ch};

            if isempty(data)
                axis off
                title(['No Data: ' channelNames{ch}]);
                continue;
            end
            % -------------------------------
            % SORT CELLS BY NUMBER OF NON-ZERO TIMEPOINTS
            % -------------------------------
            
            % Count non-zero elements per cell (row-wise)
            nonZeroCount = sum(data ~= 0);
            
            % % Handle NaNs (treat as zero)
            % nonZeroCount(isnan(nonZeroCount)) = 0;
            
            % Sort descending (most active first)
            [~, sortIdx] = sort(nonZeroCount, 'ascend');
            
            % Reorder matrix
            data = data(:,sortIdx);
            
            % -------------------------------
            % Robust color scaling
            v = data(:);
            v = v(isfinite(v));
            lo = prctile(v,1);
            hi = prctile(v,99);
            if lo == hi
                hi = lo + 1;
            end

            heatmap1((data)',[],[],[],...
                'NaNColor',[1 1 1],...
                'ColorMap',mycolormap,...
                'MinColorValue',lo,...
                'MaxColorValue',hi);

            % Title using biological names
            title(channelNames{ch},...
                'FontSize',12,...
                'FontWeight','normal');

            % ---- X Tick Handling ----
            nT = size(data,2);
            if length(tvec) >= nT
                timeCols = tvec(1:nT);
            else
                timeCols = linspace(tvec(1), tvec(end), nT);
            end

            minT = min(timeCols);
            maxT = max(timeCols);
            tickTimes = minT:intervalH:maxT;

            xtic = zeros(size(tickTimes));
            for k = 1:length(tickTimes)
                [~, idx] = min(abs(timeCols - tickTimes(k)));
                xtic(k) = idx;
            end

            set(gca,'XTick',xtic);
            set(gca,'XTickLabel',...
                arrayfun(@(x) num2str(round(x)), tickTimes,'UniformOutput',false));

            xlabel('Time (hours)');
            ylabel('Cells');

            set(gca,'FontSize',10);

            cb = colorbar('southoutside');
            ylabel(cb,'Intensity');

        end

        sgtitle(['Heatmap Time Series - ' condName],'FontSize',14);

        outBase = fullfile(figureFolder,['Heatmap_' condName]);
        print(fh,[outBase '.png'],'-dpng','-r300');
        savefig(fh,[outBase '.fig']);

        close(fh);
    end

    disp('All heatmaps generated successfully.');
