% function plot_heatmaps_by_Methoxy_posneg_v2()
% PLOT_HEATMAPS_BY_METHOXY_POSNEG (modified so each column has a correct per-column colorbar)
    close all; clc;
    addpath ('E:\rahmans4\EXPERIMENTAL DATA\NIA_Experiment\Scripts');
    %% User selects project folder
    projectFolder = uigetdir(pwd, 'Select Project folder (contains code, images, output, figure)');
    if isequal(projectFolder,0)
        disp('No folder selected. Exiting.');
        return;
    end

    outputFolder = fullfile(projectFolder, 'output_4Pxl_NucArea_CellPose');
    figureFolder = fullfile(projectFolder, 'figure');
    if ~exist(figureFolder, 'dir')
        mkdir(figureFolder);
    end

    %% Files / condition mapping
    conds = { ...
        '5xY', fullfile(outputFolder, 'stage1To2RawData.mat'); ...
        '5xM', fullfile(outputFolder, 'stage3To4RawData.mat'); ...
        'M',   fullfile(outputFolder, 'stage5To6RawData.mat') ...
        };

    % Channel names (for titles)
    channelNames = {'Methoxy-x04', 'LyspH Green', 'A\beta-AF555', 'LysoPrime deep red'};

    % Threshold for Methoxy-x04 basal
    thresholdVal = 400;

    % Tick interval (hours) for x-axis; prompt user (default 4)
    intervalH = input('Enter tick interval in hours for x-axis (default = 4): ');
    if isempty(intervalH) || ~isnumeric(intervalH) || intervalH <= 0
        intervalH = 4;
    end

    mycolormap = 'jet';
    % Colorbar geometry tweaks (normalized units, since figure units are normalized)
    cbHeight = 0.03;    % height of colorbar in normalized units
    cbYOffset = 0.02;   % extra offset below the lower subplot (normalized units)

    %% Loop through conditions
    for ci = 1:size(conds,1)
        condName = conds{ci,1};
        matFile  = conds{ci,2};

        fprintf('Processing %s ...\n', condName);

        if ~exist(matFile, 'file')
            warning('File not found: %s. Skipping %s', matFile, condName);
            continue;
        end

        S = load(matFile);

        % Validate Time_line and T1..T4
        if ~isfield(S, 'Time_line')
            warning('Time_line missing in %s. Skipping.', matFile);
            continue;
        end

        % Convert time_line to hours and make vector
        t_raw = double(S.Time_line(:));
        t_hours = t_raw ./ 60;  % adjust if Time_line is already hours

        % Extract T1..T4 and remove first 2 rows (metadata)
        T = cell(1,4);
        for ch = 1:4
            fname = sprintf('T%d', ch);
            if isfield(S, fname)
                mat = double(S.(fname));
                if size(mat,1) > 2
                    mat = mat(3:end, :);  % remove metadata rows (first 2)
                else
                    mat = zeros(0, size(mat,2));
                end
            else
                mat = zeros(0, 0);
            end
            T{ch} = mat;  % mat is cells x time
        end

        % Check T1 exists
        if isempty(T{1})
            warning('T1 missing or empty in %s. Skipping.', matFile);
            continue;
        end

        % Compute basal Methoxy-x04 = mean of first 2 timepoints per cell (columns 1:2)
        if size(T{1},2) < 2
            warning('Not enough timepoints in T1 of %s. Skipping.', matFile);
            continue;
        end
        basal = mean(T{1}(1:2,:), 1, 'omitnan');  % cells x 1

        % Partition into positive (>= threshold) and negative (< threshold)
        posIdx = basal >= thresholdVal;
        negIdx = basal < thresholdVal;

        % Optionally remove NaNs/infs from both sets
        posIdx(~isfinite(basal)) = false;
        negIdx(~isfinite(basal)) = false;

        nPos = sum(posIdx);
        nNeg = sum(negIdx);
        fprintf('  %s: %d positive, %d negative cells (threshold = %g)\n', condName, nPos, nNeg, thresholdVal);

        % Prepare combined color limits per channel (so top/bottom comparable)
        clim = cell(1,4);
        for ch = 1:4
            dataAll = T{ch};
            if isempty(dataAll)
                clim{ch} = [0 1];
                continue;
            end
            v = dataAll(:);
            v = v(isfinite(v));
            if isempty(v)
                clim{ch} = [0 1];
            else
                p = prctile(v, [1 99]);
                if p(1)==p(2)
                    p = [min(v) max(v)];
                    if p(1)==p(2)
                        p(2) = p(1) + 1;
                    end
                end
                clim{ch} = p;
            end
        end

        %% Create figure with 2 rows x 4 columns (use normalized units for simple geometry)
        fh = figure('Name', ['Heatmaps_posneg_' condName], 'Units', 'normalized', 'Position', [0.05 0.1 0.9 0.75], 'Color', 'w');

        % store axes handles to compute colorbar geometry later
        axTop = gobjects(1,4);
        axBottom = gobjects(1,4);

        for ch = 1:4
            % Extract channel data
            dataAll = T{ch};  % cells x time
            if isempty(dataAll)
                % Create placeholders for both positive and negative panels
                axTop(ch) = subplot(2,4,ch);
                axis off; title(axTop(ch), channelNames{ch}, 'FontSize', 11);
                axBottom(ch) = subplot(2,4,4+ch);
                axis off;
                continue;
            end

            % Positive group data and sorting by non-zero count (desc)
            dataPos = dataAll(:,posIdx);  % rows = selected cells
            if ~isempty(dataPos)
                nzPos = sum(dataPos ~= 0, 2);
                nzPos(isnan(nzPos)) = 0;
                [~, sPos] = sort(nzPos, 'descend');
                dataPos = dataPos(sPos, :);
            end

            % -------------------------------
            % SORT COLUMNS (timepoints) by non-zero count across rows (optional)
            nonZeroCountCols = sum(dataPos ~= 0, 1);
            [~, sortIdxCols] = sort(nonZeroCountCols, 'ascend');
            if ~isempty(sortIdxCols)
                dataPos = dataPos(:, sortIdxCols);
            end
            clear sortIdxCols nonZeroCountCols
            % -------------------------------

            %--------------------------------
            % Convert the zeros to NaN
            dataPos (dataPos == 0) = NaN;
            %--------------------------------


            % Negative group data and sorting by non-zero count (desc)
            dataNeg = dataAll(:,negIdx);
            if ~isempty(dataNeg)
                nzNeg = sum(dataNeg ~= 0, 2);
                nzNeg(isnan(nzNeg)) = 0;
                [~, sNeg] = sort(nzNeg, 'descend');
                dataNeg = dataNeg(sNeg, :);
            end

            % -------------------------------
            % SORT COLUMNS for negative
            nonZeroCountCols = sum(dataNeg ~= 0, 1);
            [~, sortIdxCols] = sort(nonZeroCountCols, 'ascend');
            if ~isempty(sortIdxCols)
                dataNeg = dataNeg(:, sortIdxCols);
            end
            clear sortIdxCols nonZeroCountCols
            % -------------------------------
            %--------------------------------
            % Convert the zeros to NaN
            dataNeg (dataNeg == 0) = NaN;
            %--------------------------------
            % Determine time columns for xticks based on length of timepoints
            nT_pos = size(dataPos,2);
            nT_neg = size(dataNeg,2);
            % prefer using nT from available data (use max)
            nT = max([nT_pos, nT_neg, 1]);

            if length(t_hours) >= nT
                timeCols = t_hours(1:nT);
            else
                % fallback evenly spaced
                if length(t_hours) >= 2
                    timeCols = linspace(t_hours(1), t_hours(end), nT);
                else
                    timeCols = 1:nT;
                end
            end

            % compute xticks based on intervalH
            minT = min(timeCols); maxT = max(timeCols);
            tickTimes = minT:intervalH:maxT;
            if isempty(tickTimes)
                tickTimes = [minT maxT];
            end
            xtic = zeros(size(tickTimes));
            for k = 1:length(tickTimes)
                [~, idxmin] = min(abs(timeCols - tickTimes(k)));
                xtic(k) = idxmin;
            end
            xticklabels = arrayfun(@(x) num2str(round(x)), tickTimes, 'UniformOutput', false);

            % -----------------------
            % Top: positive
            % -----------------------
            axTop(ch) = subplot(2,4,ch);
            if isempty(dataPos)
                text(0.5, 0.5, 'No positive cells','HorizontalAlignment','center','Parent',axTop(ch));
                axis(axTop(ch), 'off');
            else
                try
                    heatmap1((dataPos)', [], [], [], 'NaNColor', [1 1 1], ...
                        'ColorMap', mycolormap, 'MinColorValue', clim{ch}(1), 'MaxColorValue', clim{ch}(2));
                    axCur = ancestor(gca, 'axes');
                    if ~isempty(axCur)
                        axTop(ch) = axCur;
                    end
                catch
                    imagesc(dataPos); axis xy; colormap(mycolormap); caxis(clim{ch});
                    axTop(ch) = gca;
                end
                title(axTop(ch), channelNames{ch}, 'FontSize', 11, 'FontWeight','normal');
                try
                    set(axTop(ch), 'XTick', xtic, 'XTickLabel', xticklabels);
                catch
                end
                ylabel(axTop(ch), 'Cells (pos, sorted)');
            end

            % ensure top axis uses the per-channel clim
            try
                caxis(axTop(ch), clim{ch});
            catch
            end

            % -----------------------
            % Bottom: negative
            % -----------------------
            axBottom(ch) = subplot(2,4,4+ch);
            if isempty(dataNeg)
                text(0.5, 0.5, 'No negative cells','HorizontalAlignment','center','Parent',axBottom(ch));
                axis(axBottom(ch), 'off');
            else
                try
                    heatmap1((dataNeg)', [], [], [], 'NaNColor', [1 1 1], ...
                        'ColorMap', mycolormap, 'MinColorValue', clim{ch}(1), 'MaxColorValue', clim{ch}(2));
                    axCur = ancestor(gca, 'axes');
                    if ~isempty(axCur)
                        axBottom(ch) = axCur;
                    end
                catch
                    imagesc(dataNeg); axis xy; colormap(mycolormap); caxis(clim{ch});
                    axBottom(ch) = gca;
                end
                try
                    set(axBottom(ch), 'XTick', xtic, 'XTickLabel', xticklabels);
                catch
                end
                xlabel(axBottom(ch), 'Time (hours)');
                ylabel(axBottom(ch), 'Cells (neg, sorted)');
            end

            % ensure bottom axis uses the per-channel clim
            try
                caxis(axBottom(ch), clim{ch});
            catch
            end

            % small visual spacing adjustments
            try set(axTop(ch), 'FontSize', 9); end
            try set(axBottom(ch), 'FontSize', 9); end

        end % end channels

        % --- Reduce vertical gap between the two rows ---
        gapShrink = 0.04;   % try 0.02–0.06 (increase to reduce gap more)
        
        for ch = 1:4
            if isgraphics(axTop(ch)) && isgraphics(axBottom(ch))
                posTop = get(axTop(ch), 'Position');
                posBot = get(axBottom(ch), 'Position');
        
                % Move bottom row upward
                posBot(2) = posBot(2) + gapShrink;
        
                % Optionally increase bottom height slightly (keeps proportions nicer)
                posBot(4) = posBot(4) + gapShrink/2;
        
                set(axBottom(ch), 'Position', posBot);
            end
        end

        % --- Per-column colorbars (one colorbar centered under each column) -----
        colormap(fh, mycolormap);

        % gather positions for each column and create a colorbar per column
        posTop = zeros(numel(axTop),4);
        posBottom = zeros(numel(axBottom),4);
        for k = 1:numel(axTop)
            if isgraphics(axTop(k))
                posTop(k,:) = get(axTop(k), 'Position');
            else
                posTop(k,:) = [0 0 0 0];
            end
        end
        for k = 1:numel(axBottom)
            if isgraphics(axBottom(k))
                posBottom(k,:) = get(axBottom(k), 'Position');
            else
                posBottom(k,:) = [0 0 0 0];
            end
        end

        for ch = 1:4
            % compute left/right extents of the column (top and bottom)
            left_col = min(posTop(ch,1), posBottom(ch,1));
            right_col = max(posTop(ch,1)+posTop(ch,3), posBottom(ch,1)+posBottom(ch,3));
            cb_x = left_col;
            cb_w = right_col - left_col;
            % bottom y for this column (lowest y of its bottom axis)
            bottom_y = posBottom(ch,2);
            % compute colorbar y
            cb_h = cbHeight;
            cb_y = bottom_y - cbYOffset - cb_h;
            if cb_y < 0
                cb_y = 0.01;
            end

            % Create colorbar attached to the bottom axis of the column so it shares the same axes limits
            try
                cb = colorbar(axBottom(ch), 'southoutside');
                % force exact position (normalized units)
                cb.Position = [cb_x, 0.65*cb_y, cb_w, 0.5*cb_h];
                cb.Label.String = 'Intensity';
                % ensure the colorbar reflects the axis caxis (no need to set Limits directly)
            catch
                % If attaching fails, fallback to previous method but ensure caxis is set on cbax
                cbax = axes('Position', [cb_x, 0.65*cb_y, cb_w, 0.5*cb_h], 'Visible', 'off', 'Units', 'normalized');
                cb = colorbar(cbax, 'southoutside');
                cb.Position = [cb_x, cb_y, cb_w, cb_h];
                cb.Label.String = 'Intensity';
                try
                    caxis(cbax, clim{ch});
                catch
                end
                set(cbax, 'Visible', 'off');
            end
        end
        
        % --- Set font size of all axes labels and tick labels to 14 ---
        allAxes = findall(fh, 'Type', 'axes');
        
        for k = 1:length(allAxes)
            ax = allAxes(k);
        
            % Skip invisible colorbar helper axes if any
            if strcmp(get(ax,'Visible'),'off')
                continue
            end
        
            % Tick labels
            ax.FontSize = 14;
        
            % X and Y axis labels
            if isprop(ax, 'XLabel') && ~isempty(ax.XLabel)
                ax.XLabel.FontSize = 14;
            end
            if isprop(ax, 'YLabel') && ~isempty(ax.YLabel)
                ax.YLabel.FontSize = 14;
            end
        end

        % supertitle and save
        sgtitle(sprintf('Heatmaps (%s) — Top: Methoxy-x04 >= %g ; Bottom: < %g', condName, thresholdVal, thresholdVal), 'FontSize', 14);

        outBase = fullfile(figureFolder, ['Heatmap_posneg_' condName]);
        try
            set(fh, 'PaperPositionMode', 'auto');
            print(fh, [outBase '.png'], '-dpng', '-r300');
            savefig(fh, [outBase '.fig']);
            fprintf('Saved %s.png and %s.fig\n', outBase, outBase);
        catch ME
            warning('Could not save for %s: %s', condName, ME.message);
        end
        %------------------------------------------------------------------
        % --- Save figure as PDF (vector when possible) and TIFF (high-res), plus PNG and FIG ---
        outBase = fullfile(figureFolder, ['Heatmap_posneg_' condName]);
        try
            set(fh, 'PaperPositionMode', 'auto');
        
            % 1) PDF (vector if exportgraphics available)
            try
                % R2020a+: exports vector PDF when ContentType='vector'
                exportgraphics(fh, [outBase '.pdf'], 'ContentType', 'vector');
            catch
                % fallback for older MATLABs
                print(fh, [outBase '.pdf'], '-dpdf', '-r300');
            end
        
            % 2) TIFF (high-resolution raster)
            try
                % exportgraphics supports resolution for raster outputs
                exportgraphics(fh, [outBase '.tif'], 'Resolution', 300);
            catch
                % fallback: print to TIFF (raster)
                % -dtiffn produces uncompressed TIFF; change to -dtiff if you prefer compression
                print(fh, [outBase '.tif'], '-dtiffn', '-r300');
            end
        
            % 3) keep PNG and FIG saves you already used
            print(fh, [outBase '.png'], '-dpng', '-r300');
            savefig(fh, [outBase '.fig']);
        
            fprintf('Saved %s.pdf, %s.tif, %s.png and %s.fig\n', outBase, outBase, outBase, outBase);
        catch ME
            warning('Could not save for %s: %s', condName, ME.message);
        end

        close(fh);
    end

    disp('All conditions processed.');
% end