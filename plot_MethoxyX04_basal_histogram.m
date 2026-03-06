% function plot_MethoxyX04_basal_histogram()

    close all; clc;

    % ---- Select Project Folder ----
    projectFolder = uigetdir(pwd, 'Select Project folder');
    if isequal(projectFolder,0)
        disp('No folder selected. Exiting.');
        return;
    end

    outputFolder = fullfile(projectFolder,'output_4Pxl_NucArea_CellPose');
    figureFolder = fullfile(projectFolder,'figure');

    if ~exist(figureFolder,'dir')
        mkdir(figureFolder);
    end

    % ---- Condition files ----
    conds = {
        '5xY', fullfile(outputFolder,'stage1To2RawData.mat');
        '5xM', fullfile(outputFolder,'stage3To4RawData.mat');
        'M'   , fullfile(outputFolder,'stage5To6RawData.mat')
        };

    basalData = cell(1,3);

    % ---- Extract basal intensities ----
    for i = 1:3

        condName = conds{i,1};
        matFile  = conds{i,2};

        if ~exist(matFile,'file')
            warning('Missing file: %s',matFile);
            basalData{i} = [];
            continue;
        end

        S = load(matFile);

        if ~isfield(S,'T1')
            warning('T1 not found in %s',matFile);
            basalData{i} = [];
            continue;
        end

        data = double(S.T1);

        % Remove first 2 metadata rows
        if size(data,1) > 2
            data = data(3:end,:);
        else
            basalData{i} = [];
            continue;
        end

        % Ensure at least 2 time points exist
        if size(data,2) < 2
            warning('Not enough time points in %s',matFile);
            basalData{i} = [];
            continue;
        end

        % Compute basal intensity (mean of first 2 timepoints)
        basal = mean(data(1:2,:),1,'omitnan');

        % Remove NaN and Inf
        basal = basal(isfinite(basal));

        basalData{i} = basal;

        fprintf('%s: n = %d cells\n',condName,length(basal));
    end

    % ---- Plot Histogram ----
    fh = figure('Units','inches','Position',[0 0 6 5],'Color','w');

    hold on;

    colors = lines(3);  % distinct colors

    for i = 1:3
        if isempty(basalData{i})
            continue;
        end

        histogram(basalData{i}, ...
            'Normalization','probability', ...
            'DisplayStyle','stairs', ...
            'LineWidth',2, ...
            'EdgeColor',colors(i,:));
    end

    hold off;

    xlabel('Methoxy-x04 Basal Intensity (mean of first 2 time points)','FontSize',12);
    ylabel('Probability','FontSize',12);
    legend({'5xY','5xM','M'},'Location','best');
    title('Basal Methoxy-x04 Intensity Distribution','FontSize',14);

    set(gca,'FontSize',12);
    box on;

    % ---- Save Figure ----
    outFile = fullfile(figureFolder,'Histogram_MethoxyX04_Basal');
    print(fh,[outFile '.png'],'-dpng','-r300');
    savefig(fh,[outFile '.fig']);

    disp('Histogram saved successfully.');

% end