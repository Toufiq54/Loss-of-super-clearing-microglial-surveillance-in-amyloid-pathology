%%
% This script was written to extract the GFP-RelA and mScarlet-cRel data
% from the individual stage position and store in single file for further
% analysis at the same time 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Variable inforation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% T1  --> contains ratio of mean Nuclear/Total cell intensity of GFP-RelA
% T2  --> contains ratio of mean Nuclear/Total cell intensity of  mScarlet-cRel
% T3  --> contains ratio of integrated Nuclear/Total cell intensity of GFP-RelA
% T4  --> contains ratio of integrated Nuclear/Total cell intensity of  mScarlet-cRel
% T5  --> contains mean nuclear intensity of GFP-RelA
% T6  --> contains mean nuclear intensity of mScarlet-cRel
% T7  --> contains mean cellular intensity of GFP-RelA
% T8  --> contains mean cellular intensity of mScarlet-cRel
% T9  --> contains integrated nuclear intensity of GFP-RelA
% T10 --> contains integrated nuclear intensity of mScarlet-cRel
% T11 --> contains integrated cellular intensity of GFP-RelA
% T12 --> contains integrated cellular intensity of mScarlet-cRel
% T13 --> contains Nuclear Area
% T14 --> contains Cellular Area

clc
clear
close all

% Select the current working directory to the output folder
% Ask user for project directory
projDir = uigetdir(pwd, 'Select project directory (folder that contains "output")');
d = projDir;
% d = 'E:\rahmans4\EXPERIMENTAL DATA\NIA_Experiment\CRISPRi Project\01.23.2025_iF_NLS-iRFP-NLS_TRE-dCas9-KRAB_Pooled CROPseq_TNF\Tracked_images_output'; % Insert the directory of output folder
cd(d);  % Change the directory
A1 = [];     % Preallocation of Mean intensity ratio of GFP-RelA
A2 = [];     % Preallocation of Mean intensity ratio of mScarlet-cRel
A3 = [];
A4 = [];
A5 = [];
A6 = [];
A7 = [];
A8 = [];
A9 = [];
A10 = [];
A11 = [];
A12 = [];
A13 = [];
A14 = [];
% Take input from the user how many stage position want to extract
stageNi = input ('Give the initial stage number:-->');
stageNf = input ('Give the final stage number:-->');
ConFrame = input ('Give the number of frame considered:-->');
for stageID = stageNi:stageNf
    %%%%%%%%%%%%%%%%%%%% LOADING DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Loading the Ratio data
    % Loading the mean intensities
    filename2 = ['s',num2str(stageID),'MeanFluorescence','.mat'];
    L2 = load (filename2);
    MCellF1 = L2.MCellF1;
    MCellF2 = L2.MCellF2;
    MCellF3 = L2.MCellF3;
    MCellF4 = L2.MCellF4;
    % Loading the total intensities
    filename3 = ['s',num2str(stageID),'TotalFluorescence','.mat'];
    L3 = load (filename3);
    Fcell1 = L3.Fcell1;
    Fcell2 = L3.Fcell2;
    Fcell3 = L3.Fcell3;
    Fcell4 = L3.Fcell4;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%% PROCESSING THE DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Only considered Time frames are taken, rest of matrix elements are discarded
    
    MCellF1 (ConFrame:end,:) = [];
    MCellF2  (ConFrame:end,:) = [];
    MCellF3 (ConFrame:end,:) = [];
    MCellF4  (ConFrame:end,:) = [];
    
    Fcell1 (ConFrame:end,:) = [];
    Fcell2 (ConFrame:end,:) = [];
    Fcell3  (ConFrame:end,:) = [];
    Fcell4  (ConFrame:end,:) = [];

    
    % Discard all the cell response data that have missing track within the
    % experimental period. Only the cell that stay live and stay within the
    % imaging field are retained.
    % Convert NaN values in reference variable "rNf_cellf1" to zero.
    temp1 = MCellF1;
    temp1(isnan(temp1))=0;
    % Delete column having any zero element
    [~,colsNotZeros] = find (temp1(1,:)~=0);
    % colsWithZeros = any(temp1(1,:)==0);
    temp2 = temp1(:, colsNotZeros); % temp2 contains filtered cell response
    % Tracking data to individial cell
    [~,idx1] = ismember (temp2(1,:),MCellF1(1,:)); % idx1 constains the global CellID number
    
    stgID = ones(size(idx1)).*stageID;
    
    temp3 = [stgID;idx1;temp2];
    T1 ((1:size(temp3,1)), (size(A1,2)+1):size(A1,2)+size(temp3,2)) = temp3;      % r
    A1 = T1; % MCellF1
    clear temp2 temp3

    temp2 = MCellF2 (:,idx1);
    temp3 = [stgID;idx1;temp2];
    T2 ((1:size(temp3,1)), (size(A2,2)+1):size(A2,2)+size(temp3,2)) = temp3; 
    A2 = T2;     % T2 =MCellF2
    clear temp2 temp3
    
    temp2 = MCellF3 (:,idx1);
    temp3 = [stgID;idx1;temp2];
    T3 ((1:size(temp3,1)), (size(A3,2)+1):size(A3,2)+size(temp3,2)) = temp3; 
    A3 = T3;     % T3 = MCellF3
    clear temp2 temp3
    
    temp2 = MCellF4 (:,idx1);
    temp3 = [stgID;idx1;temp2];
    T4 ((1:size(temp3,1)), (size(A4,2)+1):size(A4,2)+size(temp3,2)) = temp3; 
    A4 = T4;     % T4 = MCellF4
    clear temp2 temp3
    
    temp2 = Fcell1 (:,idx1);
    temp3 = [stgID;idx1;temp2];
    T5 ((1:size(temp3,1)), (size(A5,2)+1):size(A5,2)+size(temp3,2)) = temp3; 
    A5 = T5;     % T5 = Fcell1
    clear temp2 temp3
    
    temp2 = Fcell2 (:,idx1);
    temp3 = [stgID;idx1;temp2];
    T6 ((1:size(temp3,1)), (size(A6,2)+1):size(A6,2)+size(temp3,2)) = temp3; 
    A6 = T6;     % T6 = Fcell2
    clear temp2 temp3
    
    temp2 = Fcell3 (:,idx1);
    temp3 = [stgID;idx1;temp2];
    T7 ((1:size(temp3,1)), (size(A7,2)+1):size(A7,2)+size(temp3,2)) = temp3; 
    A7 = T7;     % T7 = Fcell3
    clear temp2 temp3
    
    temp2 = Fcell4 (:,idx1);
    temp3 = [stgID;idx1;temp2];
    T8 ((1:size(temp3,1)), (size(A8,2)+1):size(A8,2)+size(temp3,2)) = temp3; 
    A8 = T8;     % T8 = Fcell4
    clear temp2 temp3
    
    
end
% Automatic detection and filter out the very noisy signal.
% a = T1 (3:end,:);
% NoisyCell = find (max(movstd(a,4))>0.65);
% Noisy_cell_info = [NoisyCell; T1(1:2,NoisyCell)];
% save (['stage',num2str(stageNi),'To',num2str(stageNf),'NoisyCell_info''.mat'],'Noisy_cell_info')
%     T1 (:,NoisyCell) = [];
%     T2 (:,NoisyCell) = [];
%     T3 (:,NoisyCell) = [];
%     T4 (:,NoisyCell) = [];
%     T5 (:,NoisyCell) = [];
%     T6 (:,NoisyCell) = [];
%     T7 (:,NoisyCell) = [];
%     T8 (:,NoisyCell) = [];
%     T9 (:,NoisyCell) = [];
%     T10 (:,NoisyCell) = [];
%     T11 (:,NoisyCell) = [];
%     T12 (:,NoisyCell) = [];
%     T13 (:,NoisyCell) = [];
%     T14 (:,NoisyCell) = [];
% 
% % Get the time series inforatmion based on the stage position
% Stimulation_Frame = 8;%input ('Give Image frame number right after stimulation -->');
% f = 7;%input ('Give the frequency of imaging cycle (in minutes) -->');
% tr = 140;%input ('Give the time (in second) remaining after before grabbing next cycle -->');
% T_stage = 40;%input ('Give the total number of positions -->');
% DownTime = 2+(30/60);%input('Give the down time During Stimulation -->');
% StageID = T1 (1,:); 
% numFrame = size(T1,1) - 2;
% for k = 1: size(StageID,2)
%     Transition_time = (StageID(k) - 1).*(f-(tr/60))/T_stage;
%     Time_line (:,k) = ((1:numFrame)' - Stimulation_Frame+1).*f + Transition_time;
% end
% Time_line (Stimulation_Frame:end,:) = Time_line (Stimulation_Frame:end,:) + DownTime;
% save (['stage',num2str(stageNi),'To',num2str(stageNf),'RawData','.mat'],'Time_line','T1','T2','T3','T4','T5','T6','T7','T8','T9','T10','T11','T12','T13','T14');

%% Manual noise reduction
% Filter the response data based on the manual visualization of data
% Plot the any response data, i.e. plot (T1(3:end,:))
% Go View/Plot Browser
% Manually select the noisy response
% Store the noise cell ID --> the column number

%%%%%%%%%%% Run this 3 lines code when filtering first time %%%%%%%%%%%%%
clear NoisyCell Noisy_cell_info
NoisyCell = []; % insert the column number of noisy response
Noisy_cell_info = [NoisyCell; T1(1:2,NoisyCell)];
save (['stage',num2str(stageNi),'To',num2str(stageNf),'NoisyCell_info''.mat'],'Noisy_cell_info')
    T1 (:,NoisyCell) = [];
    T2 (:,NoisyCell) = [];
    T3 (:,NoisyCell) = [];
    T4 (:,NoisyCell) = [];
    T5 (:,NoisyCell) = [];
    T6 (:,NoisyCell) = [];
    T7 (:,NoisyCell) = [];
    T8 (:,NoisyCell) = [];
%     T9 (:,NoisyCell) = [];
%     T10 (:,NoisyCell) = [];
%     T11 (:,NoisyCell) = [];
%     T12 (:,NoisyCell) = [];
%     T13 (:,NoisyCell) = [];
%     T14 (:,NoisyCell) = [];
% % 
% % 
% %%%%%%%%%%%% Run this 3 line code if you alread have noisy response %%%%%%
% L = load (['stage',num2str(stageNi),'To',num2str(stageNf),'NoisyCell_info''.mat']);
% NoisyCell = L.Noisy_cell_info;
% if length(NoisyCell)~=0
%     NoisyCell(2:3,:)=[];
% 
%     T1 (:,NoisyCell) = [];
%     T2 (:,NoisyCell) = [];
%     T3 (:,NoisyCell) = [];
%     T4 (:,NoisyCell) = [];
%     T5 (:,NoisyCell) = [];
%     T6 (:,NoisyCell) = [];
%     T7 (:,NoisyCell) = [];
%     T8 (:,NoisyCell) = [];
%     T9 (:,NoisyCell) = [];
%     T10 (:,NoisyCell) = [];
%     T11 (:,NoisyCell) = [];
%     T12 (:,NoisyCell) = [];
%     T13 (:,NoisyCell) = [];
%     T14 (:,NoisyCell) = [];
% else
%     T1 = T1;
%     T2 = T2;
%     T3 = T3;
%     T4 = T4;
%     T5 = T5;
%     T6 = T6;
%     T7 = T7;
%     T8 = T8;
%     T9 = T9;
%     T10 = T10;
%     T11 = T11;
%     T12 = T12;
%     T13 = T13;
%     T14 = T14;
% end
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the time series inforatmion based on the stage position
Stimulation_Frame = input ('Give Image frame number right after stimulation -->');
f = input ('Give the frequency of imaging cycle (in minutes) -->');
tr = input ('Give the time (in second) remaining after before grabbing next cycle -->');
T_stage = input ('Give the total number of positions -->');
DownTime = input('Give the down time During Stimulation -->');
StageID = T1 (1,:); 
numFrame = size(T1,1) - 2;
for k = 1: size(StageID,2)
    Transition_time = (StageID(k) - 1).*(f-(tr/60))/T_stage;
    Time_line (:,k) = ((1:numFrame)' - Stimulation_Frame+1).*f + Transition_time;
end
Time_line (Stimulation_Frame:end,:) = Time_line (Stimulation_Frame:end,:) + DownTime;
save (['stage',num2str(stageNi),'To',num2str(stageNf),'RawData','.mat'],'Time_line','T1','T2','T3','T4','T5','T6','T7','T8');

% SMOOTHING DATA (MINIMAL SMOOTHING)
ST1 = T1;
ST2 = T2;
ST3 = T3;
ST4 = T4;
ST5 = T5;
ST6 = T6;
ST7 = T7;
ST8 = T8;
% ST9 = T9;
% ST10 = T10;
% ST11 = T11;
% ST12 = T12;
% ST13 = T13;
% ST14 = T14;

if Stimulation_Frame>5
% preallocation 
ST1 = T1;
ST2 = T2;
ST3 = T3;
ST4 = T4;
ST5 = T5;
ST6 = T6;
ST7 = T7;
ST8 = T8;
% ST9 = T9;
% ST10 = T10;
% ST11 = T11;
% ST12 = T12;
% ST13 = T13;
% ST14 = T14;
sp = 0.04;
for i = 1:size(T1,2)
    ST1(3:end,i) = smooth (T1(3:end,i),sp,'lowess');
    ST2(3:end,i) = smooth (T2(3:end,i),sp,'lowess');
    ST3(3:end,i) = smooth (T3(3:end,i),sp,'lowess');
    ST4(3:end,i) = smooth (T4(3:end,i),sp,'lowess');
    ST5(3:end,i) = smooth (T5(3:end,i),sp,'lowess');
    ST6(3:end,i) = smooth (T6(3:end,i),sp,'lowess');
    ST7(3:end,i) = smooth (T7(3:end,i),sp,'lowess');
    ST8(3:end,i) = smooth (T8(3:end,i),sp,'lowess');
    % ST9(3:end,i) = smooth (T9(3:end,i),sp,'lowess');
    % ST10(3:end,i) = smooth (T10(3:end,i),sp,'lowess');
    % ST11(3:end,i) = smooth (T11(3:end,i),sp,'lowess');
    % ST12(3:end,i) = smooth (T12(3:end,i),sp,'lowess');
    % ST13(3:end,i) = smooth (T13(3:end,i),sp,'lowess');
    % ST14(3:end,i) = smooth (T14(3:end,i),sp,'lowess');

end
else
    % preallocation 

sp = 0.06;
for i = 1:size(T1,2)
    ST1((Stimulation_Frame+1):end,i) = smooth (T1((Stimulation_Frame+1):end,i),sp,'lowess');
    ST2((Stimulation_Frame+1):end,i) = smooth (T2((Stimulation_Frame+1):end,i),sp,'lowess');
    ST3((Stimulation_Frame+1):end,i) = smooth (T3((Stimulation_Frame+1):end,i),sp,'lowess');
    ST4((Stimulation_Frame+1):end,i) = smooth (T4((Stimulation_Frame+1):end,i),sp,'lowess');
    ST5((Stimulation_Frame+1):end,i) = smooth (T5((Stimulation_Frame+1):end,i),sp,'lowess');
    ST6((Stimulation_Frame+1):end,i) = smooth (T6((Stimulation_Frame+1):end,i),sp,'lowess');
    ST7((Stimulation_Frame+1):end,i) = smooth (T7((Stimulation_Frame+1):end,i),sp,'lowess');
    ST8((Stimulation_Frame+1):end,i) = smooth (T8((Stimulation_Frame+1):end,i),sp,'lowess');
    % ST9((Stimulation_Frame+1):end,i) = smooth (T9((Stimulation_Frame+1):end,i),sp,'lowess');
    % ST10((Stimulation_Frame+1):end,i) = smooth (T10((Stimulation_Frame+1):end,i),sp,'lowess');
    % ST11((Stimulation_Frame+1):end,i) = smooth (T11((Stimulation_Frame+1):end,i),sp,'lowess');
    % ST12((Stimulation_Frame+1):end,i) = smooth (T12((Stimulation_Frame+1):end,i),sp,'lowess');
    % ST13((Stimulation_Frame+1):end,i) = smooth (T13((Stimulation_Frame+1):end,i),sp,'lowess');
    % ST14((Stimulation_Frame+1):end,i) = smooth (T14((Stimulation_Frame+1):end,i),sp,'lowess');
end
end
save (['stage',num2str(stageNi),'To',num2str(stageNf),'SmoothedData','.mat'],'Time_line','ST1','ST2','ST3','ST4','ST5','ST6','ST7','ST8');
%%










% Parameter Extraction 
% Use function to find out the Amplitude, Duration, AUC, Time to peak 
% Nuclear occupancy1 = area under [nuc:total ratio -1(basal ratio)] curve (unit: ratio*hr, ignore negative area)
% Nuclear occupancy2 = area under [gfp_nuc-gfp_nuc(0)] curve (unit: intensity*hr, ignore negative area) 
% Time-averaged occupancy1 = nuc occ1/[duration of time course]
% Time-averaged occupancy2 = nuc occ2/[duration of time course]
% Magnitude1 = peak height of nuc:total ratio in first cycle 
% Magnitude2 = peak height of gfp_nuc in first cycle
% Duration of first cycle in nuc:total ratio scale (hrs)
% Time to first peak of nuc:total ratio (hrs)

% rNucOcc1 --> Area under [nuc:total ratio of GFP-RelA - basal], negative area ignored
% rNucOcc2 --> Area under [nuc:total ratio of mScarlet-cRel - basal], negative area ignored
% tNucOcc1 --> Area under [GFP-RelA - basal], negative area ignored
% tNucOcc2 --> Area under [mScarlet-cRel - basal], negative area ignored

% rPH1 --> peak height of nuc:total ratio in first cycle of GFP-RelA
% rPH2 --> peak height of nuc:total ratio in first cycle of mScarlet-cRel
% tPH1 --> peak height of gfp_nuc in first cycle of GFP-RelA
% tPH2 --> peak height of gfp_nuc in first cycle of mScarlet-cRel

% rDur1 --> Duration of first cycle in nuc:total ratio scale (hrs) of GFP-RelA
% rDur2 --> Duration of first cycle in nuc:total ratio scale (hrs) of mScarlet-cRel

% rTFP1 --> Time to first peak of nuc:total ratio (hrs) of GFP-RelA
% rTFP2 --> Time to first peak of nuc:total ratio (hrs) of mScarlet-cRel






