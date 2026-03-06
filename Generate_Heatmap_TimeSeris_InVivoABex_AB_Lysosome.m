% function generate_heatmap_timeseries()
% GENERATE_HEATMAP_TIMESERIES
%   GUI to choose project folder. Loads mat files from <Project>/output and
%   generates heatmap time-series (T1..T4) for three conditions:
%     5xY -> stage1To2RawData.mat   (stages 1,2)
%     5xM -> stage3To4RawData.mat   (stages 3,4)
%     M   -> stage5To6RawData.mat   (stages 5,6)
%
%   Assumptions:
%     - Each mat file contains variables T1..T8 and Time_line
%     - Each T* is H x W or cells x time matrix (we treat them as cells x time)
%     - First two rows of each T contain metadata and will be removed
%
%   Output:
%     - Figures saved to <Project>/figure as PNG and FIG

    close all;
    clc;

    %% Ask user for project folder
    projectFolder = uigetdir(pwd, 'Select Project folder (contains code, images, output, figure)');
    if isequal(projectFolder,0)
        disp('No folder selected. Exiting.');
        return;
    end

    % Define important sub-folders
    outputFolder = fullfile(projectFolder, 'output_4Pxl_NucArea_CellPose');
    figureFolder = fullfile(projectFolder, 'figure');

    % Create figure folder if not exist
    if ~exist(figureFolder, 'dir')
        mkdir(figureFolder);
    end

    % Map conditions to filenames and names
    conds = {
        '5xY', fullfile(outputFolder, 'stage1To2RawData.mat');
        '5xM', fullfile(outputFolder, 'stage3To4RawData.mat');
        'M'   , fullfile(outputFolder, 'stage5To6RawData.mat')
        };

    % Ask for tick interval in hours (default 4)
    prompt = 'Enter tick interval in hours for x-axis (default = 4): ';
    intervalH = input(prompt);
    if isempty(intervalH) || ~isnumeric(intervalH) || intervalH <= 0
        intervalH = 4;
    end

    % Colors/limits defaults (you may change inside loop)
    defaultMinColor = []; % leave empty to auto compute per-channel
    defaultMaxColor = []; % leave empty to auto compute per-channel

    % Loop over conditions
    for ci = 1:size(conds,1)
        condName = conds{ci,1};
        matFile   = conds{ci,2};

        fprintf('Processing condition: %s\n', condName);

        if ~exist(matFile, 'file')
            warning('File not found: %s. Skipping %s.', matFile, condName);
            continue;
        end

        % Load mat file into struct so we don't pollute workspace
        S = load(matFile);

        % Validate Time_line exists
        if ~isfield(S, 'Time_line')
            warning('Time_line variable missing in %s. Skipping.', matFile);
            continue;
        end

        % Convert time to hours (divide by 60)
        t_raw = S.Time_line;
        % Time_line might be column or row; make a column vector
        if ismatrix(t_raw)
            tvec = t_raw(:) ./ 60; % hours
        else
            tvec = double(t_raw(:)) ./ 60;
        end

        % Prepare channel matrices T1..T4
        T = cell(1,4);
        for ch = 1:4
            fieldname = sprintf('T%d', ch);
            if ~isfield(S, fieldname)
                warning('Variable %s missing in %s. Filling with empty.', fieldname, matFile);
                T{ch} = [];
                continue;
            end
            mat = double(S.(fieldname)); %#ok<UDIM>
            % Expect mat is (rows x timepoints). Remove first two rows if exist
            if size(mat,1) > 2
                mat = mat(3:end, :); % delete first 2 rows
            else
                % If there are <=2 rows, result is empty
                mat = zeros(0, size(mat,2));
            end
            T{ch} = mat;
        end

        % Prepare figure
        fh = figure('Name', ['Heatmaps - ' condName], 'Units','inches', 'Position',[0 0 14 4], 'Color','w');
        % 1x4 subplots
        for ch = 1:4
            subplot(1,4,ch);
            ax = gca;

            mat = T{ch};
            if isempty(mat)
                % show placeholder
                text(0.5,0.5, sprintf('No data: T%d', ch), 'HorizontalAlignment','center');
                axis off;
                title(sprintf('Channel %d (T%d)', ch, ch), 'FontSize',12);
                continue;
            end

            % imagesc expects rows = cells, columns = timepoints
            imagesc(mat);
            axis xy; % so first row at top
            colormap(ax, jet);
            % Compute color limits if not provided
            if isempty(defaultMinColor) || isempty(defaultMaxColor)
                % robust limits: 1st and 99th percentile to avoid outliers
                v = mat(:);
                vfinite = v(isfinite(v));
                if isempty(vfinite)
                    clim = [0 1];
                else
                    lo = prctile(vfinite, 1);
                    hi = prctile(vfinite, 99);
                    if lo==hi
                        lo = min(vfinite);
                        hi = max(vfinite);
                        if lo==hi
                            hi = lo + 1;
                        end
                    end
                    clim = [lo hi];
                end
            else
                clim = [defaultMinColor defaultMaxColor];
            end
            caxis(clim);

            % Labels and title
            title(sprintf('Channel %d', ch), 'FontSize',12, 'FontWeight','normal');
            xlabel('Time (hours)', 'FontSize',10);
            ylabel('Cell index', 'FontSize',10);

            % x ticks: map from time vector to column indices.
            nt = size(mat,2);
            if length(tvec) >= nt
                % find best mapping: assume tvec length >= nt and times correspond to columns
                timeCols = tvec(1:nt);
            elseif length(tvec) == 1
                % single time value (unlikely): create equal spaced
                timeCols = (0:nt-1) * (tvec);
            else
                % fallback: create linear times 1..nt
                timeCols = 1:nt;
            end

            % Compute tick positions based on intervalH
            minT = min(timeCols);
            maxT = max(timeCols);
            tickPositionsH = minT:intervalH:maxT;
            if isempty(tickPositionsH)
                tickPositionsH = [minT, maxT];
            end

            % For each tick time find nearest column index
            xtic = zeros(size(tickPositionsH));
            for k = 1:length(tickPositionsH)
                [~, idxmin] = min(abs(timeCols - tickPositionsH(k)));
                xtic(k) = idxmin;
            end

            % Set xticks if they are inside range
            xtick(xtic); %#ok<*ASGLU>
            xticks(xtic);
            xticklabels(arrayfun(@(x) num2str(x), round(tickPositionsH), 'UniformOutput', false));
            set(gca,'FontSize',10);
            % colorbar for each subplot placed at bottom-ish for compact layout, but show only for first two perhaps
            cb = colorbar('southoutside');
            % shrink colorbar a bit
            cbPos = cb.Position;
            cb.Position = [cbPos(1), cbPos(2)*0.95, cbPos(3), cbPos(4)*0.6];
            ylabel(cb, sprintf('Intensity (T%d)', ch));
        end

        % Save figure: PNG and FIG
        outBase = fullfile(figureFolder, ['Heatmap_' condName]);
        try
            % Increase paper size/resolution if desired
            set(fh, 'PaperPositionMode','auto');
            print(fh, [outBase '.png'], '-dpng', '-r300');
            savefig(fh, [outBase '.fig']);
            fprintf('Saved figures: %s.png and %s.fig\n', outBase, outBase);
        catch ME
            warning('Failed to save figure for %s: %s', condName, ME.message);
        end

        % Close figure to keep environment tidy
        close(fh);
    end

    disp('All conditions processed.');
% end