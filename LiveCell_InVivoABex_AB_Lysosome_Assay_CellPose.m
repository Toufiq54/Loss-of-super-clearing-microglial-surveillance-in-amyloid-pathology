%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by: Shah Md Toufiqur Rahman
% Date: 03.04.2026 (updated)
% Research Fellow
% Transcription Systems Dynamics and Biology Unit
% Laboratory of Molecular Biology and Immunology 
% National Institute on Aging, NIH
% National Institute of Health
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------- EXPT. DETAILS ------------------------------- 
% Experiment date: 3/2/2026
% Sample	            5xY	    5xM	     M
% Stage position	    1,2	    3,4	    5,6
% 1. File name invivoAB_exvivoAB_lysosome_R1_2026_03_02__10_33_51 and invivoAB_exvivoAB_lysosome_R2_2026_03_02__18_04_58
% 2. microglia-Replicate-channel
%   •	5xFAD young 5xY
%   •	5xFAD (+) middle-age  5xM
%   •	5xFAD (-) middle-age  M
% 3.	Channels
%   •	Channel 1 – Methoxy-x04
%   •	Channel 2 – LyspH Green (activity)
%   •	Channel 3 – Aβ-AF555
%   •	Channel 4 – LysoPrime deep red (amount)
%   •	Channel 5 – Brightfield
% 
% 4. 2 cycles before adding Aβ-AF555 treatment
% 35 cycles after adding Aβ-AF555 treatment
% 2x2 tiling, 7 min time interval, 2 min down time, 105 s waiting time 
%% Check the cell segmentation 
%--------------------------------------------------------------------------
clc
clear
close all
%--------------------------------------------------------------------------
% Adding the directory of important function that will be used here.
addpath ('E:\rahmans4\EXPERIMENTAL DATA\NIA_Experiment\Scripts');
%--------------------------------------------------------------------------
% -----------READING THE INPUT IMAGE FILE (sStageNumber.tif)---------------
% % Take input from the user : Number of stage position want to analyse
% stageN = input ('Give the number of stage position for the experiment-->');
% Take input from the user how many stage position want to extract
% stageNi = input ('Give the initial stage number:-->');
% stageNf = input ('Give the final stage number:-->');
%==========================================================================
% --------------------Get Image directory----------------------------------
disp('Select project directory (folder that contains "images")');
% Ask user for project directory
projDir = uigetdir(pwd, 'Select project directory (folder that contains "images")');
currentFolder = projDir;
ImgName = '';
AnalysisFolder = 'Tracked_images_CellPose';
% make output directory
cd (currentFolder);
cd ..\.
resultFolder = fullfile(pwd,['\',AnalysisFolder,'_output_CellPose']);
if ~exist(resultFolder, 'dir')
    mkdir(resultFolder);
end
% resultFolder = 'E:\NIA_Experiment\Expt# MEF_GFP-RelA_LPS_100ng_ml\output';
%==========================================================================
cd (currentFolder);
for stageID = 2:2 %stageNi:stageNf
    %-----------------------------------------------------
    stage_name = ['s',num2str(stageID),'-4'];       % Channel 4 – LysoPrime deep red (amount)
    %---------------------------------------------------
    % Make a folder as named by the stage name
    outputFolder = fullfile(pwd, AnalysisFolder,[ImgName,'s',num2str(stageID)]);
    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end
% Reading live-cell image movie file (GFP-image-stack)   
    fname = [stage_name,'.tif'];
    info = imfinfo (fname);
%--------------------------------------------------------------------------
    NumFrame = numel(info);    % or Assian number of Frame
    %NumFrame = 50;       % USER DEFINED IMAGE FRAME NUMBER
%--------------------------------------------------------------------------
%   Segmentation of image sequence. This section will make new folder for
%   each of the stage position and save the segmented image sequence in that
%   folder. 
   
   for imageID = 1:1% NumFrame
        % Read cell segmentation binary image
        imgCell = imread (fname, imageID);
        %----------------------------------------
        % load Cellpose models
        cpCyto2 = cellpose(Model="cyto2");
        cpNuclei = cellpose(Model="nuclei");
        %----------------------------------------
        averageCellDiameter = 40;   % Adjustable parameters
        %----------------------------------------
        labelsDefault = segmentCells2D(cpCyto2,imgCell,ImageCellDiameter=averageCellDiameter);  % using the cpCyto2 model to segment the cells
        loverlayDefault = labeloverlay(imgCell,labelsDefault);
        figure;imshow(loverlayDefault);
        % show the distribution of cellular area
        stat = regionprops (labelsDefault,'Area');
        A = cell2mat(struct2cell(stat));
        figure; subplot (1,2,1);histogram(A);xlabel ('Cellular area (pixels)'); ylabel ('# Cells');
        subplot (1,2,2); imagesc(labelsDefault);
        Arange = input ('Provide Cellular area pixel range as [minA, maxA]:-');
        %--------------------------------------------------
        labelIdx = label2idx(labelsDefault);
        labelsFiltered = zeros(size(labelsDefault),'like',labelsDefault);
        newID = 1;
        for i = 1:numel(labelIdx)
            A = numel(labelIdx{i});
            if A >= Arange(1) && A <= Arange(2)
                labelsFiltered(labelIdx{i}) = newID;
                newID = newID + 1;
            end
        end
        labelsDefault = labelsFiltered;
        %--------------------------------------------------
        stat = regionprops (labelsDefault,'Area');
        A = cell2mat(struct2cell(stat));
        figure, subplot (1,2,1);histogram(A); xlabel ('Cellular Area (pixels)'); ylabel ('# Cells');
        subplot (1,2,2); imagesc(labelsDefault);
    end  % for check 
end
%% Cellular area segmentation and pre-processing for all stage position
%--------------------------------------------------------------------------
% Take input from the user how many stage position want to extract
stageNi = input ('Give the initial stage number:-->');
stageNf = input ('Give the final stage number:-->');
%--------------------------------------------------------------------------
cd (currentFolder);
for stageID = stageNi:stageNf
    disp(['Now processing Cellular segmentation and pre-processing of position #',num2str(stageID)]);
    %-----------------------------------------------------
    stage_name = ['s',num2str(stageID),'-4'];       % Channel 4 – LysoPrime deep red (amount)
    %---------------------------------------------------
    % Make a folder as named by the stage name
    outputFolder = fullfile(pwd, AnalysisFolder,[ImgName,'s',num2str(stageID)]);
    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end
    % Reading live-cell image movie file 
    fname = [stage_name,'.tif'];
    info = imfinfo (fname);
    %----------------------------------------------------------------------
    NumFrame = numel(info);    % or Assian number of Frame
    %NumFrame = 50;       % USER DEFINED IMAGE FRAME NUMBER
    d = numel(num2str(NumFrame));
%--------------------------------------------------------------------------
%   Segmentation of image sequence. This section will make new folder for
%   each of the stage position and save the segmented image sequence in that
%   folder. 
   for imageID = 1:NumFrame
% Read cell segmentation binary image
        imgCell = imread (fname, imageID);
        %----------------------------------------
        % load Cellpose models
        cpCyto2 = cellpose(Model="cyto2");
        cpNuclei = cellpose(Model="nuclei");
        %----------------------------------------
        averageCellDiameter = 40;   % Adjustable parameters
        %----------------------------------------
        labelsDefault = segmentCells2D(cpCyto2,imgCell,ImageCellDiameter=averageCellDiameter);
        loverlayDefault = labeloverlay(imgCell,labelsDefault);
        %--------------------------------------------
        labelIdx = label2idx(labelsDefault);
        labelsFiltered = zeros(size(labelsDefault),'like',labelsDefault);
        newID = 1;
        for i = 1:numel(labelIdx)
            A = numel(labelIdx{i});
            if A >= Arange(1) && A <= Arange(2)
                labelsFiltered(labelIdx{i}) = newID;
                newID = newID + 1;
            end
        end
        labelsDefault = labelsFiltered;
        %--------------------------------------------
        CellMask (:,:,imageID) = labelsDefault;

        % Change the directory to created folder
        cd (outputFolder)
%       Save the labelled, segmented cell 
        imwrite (uint16(labelsDefault),['segmented_image_',sprintf(['%0',num2str(d),'d'],(imageID)),'.tif'])
%       Change the folder directory back to original location
        cd (currentFolder);
   end
    %------------------- Saving as GIF file -------------------------------
    % cd (currentFolder);     
    % % save (['Labels_Nuc','_s',num2str(stageID),'.mat'],'NucMask');
    % save(['Labels_Cell','_s',num2str(stageID),'.mat'],'CellMask','-v7.3')
    % save_3D_matrix_as_gif_v2(['Labels_Cell','_s',num2str(stageID),'.gif'],CellMask)
    %------------------- Tracking cells-------------------------------------
    cd (outputFolder);
    ctrack_dir = 'LineageMapper.exe';
    cd('E:\rahmans4\EXPERIMENTAL DATA\NIA_Experiment\Scripts\LM\Lineage-Mapper-compiled-matlab');
    command = [ctrack_dir,' ','segmented_images_path',' ','"',outputFolder,'"',' ','tracked_images_path',' ','"',outputFolder,'"'];
    status = jsystem (command);
    cd (currentFolder); 
   %-----------------------------------------------------------------------
end
%% %%%%%%%%%%%%%%%%%%%%%%%------ANALYSIS--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for stageID = stageNi:stageNf
    disp(['Now Analysing of position #',num2str(stageID)]);
    cd (currentFolder);
    stage_name = [ImgName,'s',num2str(stageID)];
    
    % Make a folder as named by the stage name
    outputFolder = fullfile(pwd, AnalysisFolder,[ImgName,'s',num2str(stageID)]);
    if ~exist(outputFolder, 'dir')
      mkdir(outputFolder);
    end

    % make output directory
    cd (currentFolder);
    cd ..\.
    resultFolder = fullfile(pwd,['\','output_4Pxl_NucArea_CellPose']);
    if ~exist(resultFolder, 'dir')
        mkdir(resultFolder);
    end
    
    %-------------------------------------------
    % Reading live-cell image movie file
    stage_name1 = ['s',num2str(stageID),'-1'];       %   •	Channel 1 – Methoxy-x04    
    fname1 = [stage_name1,'.tif'];    
    % info1 = imfinfo (fname1);         
    stage_name2 = ['s',num2str(stageID),'-2'];       %   •	Channel 2 – LyspH Green (activity)      
    fname2 = [stage_name2,'.tif'];    
    % info2 = imfinfo (fname2);
    stage_name3 = ['s',num2str(stageID),'-3'];       %   •	Channel 3 – Aβ-AF555     
    fname3 = [stage_name3,'.tif'];    
    % info3 = imfinfo (fname3);
    stage_name4 = ['s',num2str(stageID),'-4'];       %   •	Channel 4 – LysoPrime deep red (amount)     
    fname4 = [stage_name4,'.tif'];    
    % info4 = imfinfo (fname4);
    %--------------------------------------------
    % Read the position.csv file which have all the centroid positions for all
    % the tracks
    cd (outputFolder);
    T = xlsread ('positions.csv');      % Tracking information from the Lineage mapper
    counter = 0;
    d = numel(num2str(NumFrame));

for imageID = 1: NumFrame
    cd (currentFolder);
    %----------------------------------------------------------------------
    % Reading all the channels 
    imgMain1 = imread (fname1, imageID); % Channel 1 – Methoxy-x04
    imgMain1 = imgMain1 - ProcessedBG1(imgMain1);  %imgaussian (imgMain1,500,550);
    imgMain1 = medfilt2 (imgMain1);

    imgMain2 = imread (fname2, imageID); % Channel 2 – LyspH Green (activity)
    imgMain2 = imgMain2 - ProcessedBG1(imgMain2);  %imgaussian (imgMain2,500,550);
    imgMain2 = medfilt2(imgMain2);
    
    imgMain3 = imread (fname3, imageID); % Channel 3 – Aβ-AF555  
    imgMain3 = imgMain3 - ProcessedBG1(imgMain3);  %imgaussian (imgMain2,500,550);
    imgMain3 = medfilt2(imgMain3);

    imgMain4 = imread (fname4, imageID); % Channel 4 – LysoPrime deep red (amount)    
    imgMain4 = imgMain4 - ProcessedBG1(imgMain4);  %imgaussian (imgMain2,500,550);
    imgMain4 = medfilt2(imgMain4);
    %----------------------------------------------------------------------
    cd (outputFolder);   
    mask = imread ([sprintf(['%0',num2str(d),'d'],(imageID)),'.tif']);  % cell mask tracking from LM toolbox
%
%%%% Save the Segmented image overleyed with either GFP-RelA or mScarlet-cRel%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SI = imfuse (imgMain2,bwperim(mask));   % overlay the cell mask in original GFP-RelA image
    % SI = imfuse (SI,bwperim(img2));         % adding overlay of the numclear mask 
    MSI (:,:,:,imageID) = SI;               % Saving the overlayed sequence. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SIT = SI;

    cellN = sum(T(:,2) == imageID);     % number of cell in each frame
    id = find(T(:,2)==imageID);
    TT = T(id(1):id(end),:);

    for cellID = 1: cellN
        counter = counter +1;
        cell = (mask == TT(cellID,1));  % selecting individual cell
        cell_label = TT(cellID,1);
        s = regionprops (cell,'BoundingBox','Centroid');
        if ~isempty(s)
            bb = s.BoundingBox;
            C_mask = extend_im (cell,bb,0);          % crop the tracked image for individual cell mask
            OI_crop1 = extend_im (imgMain1,bb,0);    % channel 1 
            OI_crop2 = extend_im (imgMain2,bb,0);    % channel 2
            OI_crop3 = extend_im (imgMain3,bb,0);    % channel 3 
            OI_crop4 = extend_im (imgMain4,bb,0);    % channel 4
                       
            % Adding cell number to the overlayed image
            cen = s.Centroid;
            [row,col] = find(bwperim(cell),1);
            SIT = insertText (SI,cen,num2str(cell_label),'FontSize',18,'BoxColor',...
                                {'yellow'},'BoxOpacity',0,'TextColor','white');
            SI = SIT;
            
            %%%%%%%%%%%% Generate and save the overlayed image %%%%%%%%%%%%%%%
            % Get the boundary of segmented cell
            B = bwboundaries(cell);
            C_mask = C_mask>0;  % binearize the cellular area
            C_mask = bwareafilt (C_mask,1);
            % Calculate the area of nuclear, and cell 
            Acell = bwarea (C_mask);    % Area of cell
            % Saving the area information 
            ACell (TT(cellID,2),TT(cellID,1)) = Acell;

            % mean fluorescence intensity  calculation 
            fcell1 = mean(mean(nonzeros((double(C_mask)).*(double(OI_crop1)))));
            fcell2 = mean(mean(nonzeros((double(C_mask)).*(double(OI_crop2)))));
            fcell3 = mean(mean(nonzeros((double(C_mask)).*(double(OI_crop3)))));
            fcell4 = mean(mean(nonzeros((double(C_mask)).*(double(OI_crop4)))));
    
            % Saving the mean fluorescence measurement using the regionprops
            MCellF1(TT(cellID,2),TT(cellID,1)) = fcell1;
            MCellF2(TT(cellID,2),TT(cellID,1)) = fcell2;
            MCellF3(TT(cellID,2),TT(cellID,1)) = fcell3;
            MCellF4(TT(cellID,2),TT(cellID,1)) = fcell4;
            
            % Saving the total fluorescence intensities
            Fcell1(TT(cellID,2),TT(cellID,1)) = sum (nonzeros(double(C_mask).*(double(OI_crop1))));
            Fcell2(TT(cellID,2),TT(cellID,1)) = sum (nonzeros(double(C_mask).*(double(OI_crop2))));
            Fcell3(TT(cellID,2),TT(cellID,1)) = sum (nonzeros(double(C_mask).*(double(OI_crop3))));
            Fcell4(TT(cellID,2),TT(cellID,1)) = sum (nonzeros(double(C_mask).*(double(OI_crop4))));

        end

        
    end
    clear TT id cellN
    % SAVE THE SEGMENTED OVERLEY WITH CELL NUMBER
    MSIT (:,:,:,imageID) = SIT;
        
end
clear T counter
cd (resultFolder);
% Saving the data in .xlsx file
xlswrite ([stage_name,'TotalFluorescence.xlsx'],Fcell1, 'Channel_1')
xlswrite ([stage_name,'TotalFluorescence.xlsx'],Fcell2, 'Channel_2')
xlswrite ([stage_name,'TotalFluorescence.xlsx'],Fcell3, 'Channel_3')
xlswrite ([stage_name,'TotalFluorescence.xlsx'],Fcell4, 'Channel_4')
% Saving the data in .mat file
save([stage_name,'TotalFluorescence.mat'],'Fcell1','Fcell2',"Fcell3",'Fcell4');

% Saving the data in .xlsx file
xlswrite ([stage_name,'MeanFluorescence.xlsx'],MCellF1, 'Channel_1')
xlswrite ([stage_name,'MeanFluorescence.xlsx'],MCellF2, 'Channel_2')
xlswrite ([stage_name,'MeanFluorescence.xlsx'],MCellF3, 'Channel_3')
xlswrite ([stage_name,'MeanFluorescence.xlsx'],MCellF4, 'Channel_4')
% Saving the data in .mat file
save([stage_name,'MeanFluorescence.mat'],'MCellF1','MCellF2','MCellF3','MCellF4');

xlswrite ([stage_name,'Area.xlsx'],ACell,'CellArea')
save([stage_name,'Area.mat'],'ACell');

% SAVING THE GIF FILE
save_stack_as_avi ([stage_name,'overlay','.gif'],MSI,4);
save_stack_as_avi ([stage_name,'overlay_withCellNum','.gif'],MSIT,4);

cd(currentFolder);      % Changing directory to main folder again

clear ACell 
clear MCellF1 MCellF2 MCellF3 MCellF4 
clear Fcell1 Fcell2  Fcell3 Fcell4
       
end





