%% Load coordinates:
AnimalNo = [1 2 7 8 9 10 19 20 31 32]; % Enter all animals of one group
SliceNo = [2 3 4]; % Enter all levels you want to plot per animal (we work with left and right striatum separately)
DataDir = fullfile(cd,'OriginalData'); % Enter Input folder, must contain 2 files per animal, level and side: 
    % 1 txt with outline coordinates ("AnimalNo_treatment_lesiontype_brainregion_level_hemisphere_Outline.txt")
    % 1 csv with cell coordinates ("AnimalNo_treatment_lesiontype_brainregion_level_hemisphere_Cells.csv")
Lesion = "OHDA";    % Select lesion type.
Drug = "LD";        % Select drug treatment.
Side = "L"          % Select hemisphere.

fnames = cell(numel(AnimalNo),numel(SliceNo));
fnames_outline = cell(numel(AnimalNo),numel(SliceNo));
clear T T_outlines
for iAnimal = 1:numel(AnimalNo)
    for iSlice = 1:numel(SliceNo)
        % Read coordinates of activated cells:
        fnames{iAnimal,iSlice} = ...
            sprintf('%s%s%d_%s_%s_Striatum_%d_%s_Cells.csv', ...% Enter treatment ("LD24"), surgery type ("6OHDA") and hemisphere ("L") for each group
            DataDir,filesep,AnimalNo(iAnimal),Drug,Lesion,...
            SliceNo(iSlice),Side);
        T{iAnimal,iSlice} = readtable(fnames{iAnimal,iSlice});

        % Read outline:
        fnames_outline{iAnimal,iSlice} = ...
            sprintf('%s%s%d_LD24_6OHDA_Striatum_%d_L_Outline.txt', ... % Enter treatment name, surgery type and hemisphere for each group
            DataDir, filesep, AnimalNo(iAnimal), SliceNo(iSlice));
        fnames_outline{iAnimal,iSlice} = ...
            sprintf('%s%s%d_%s_%s_Striatum_%d_%s_Outline.txt', ... % Enter treatment name, surgery type and hemisphere for each group
            DataDir, filesep, AnimalNo(iAnimal),Drug,Lesion,...
            SliceNo(iSlice),Side);
        T_outlines{iAnimal,iSlice} = readtable(fnames_outline{iAnimal,iSlice});
        T_outlines{iAnimal,iSlice}.Properties.VariableNames = {'X' 'Y'};
    end
end

%% Plot and align striatal outlines (open levels after each other and align them individually if necessary)

% Manual adjustment of rotation (First transformation matrix)
    % Vector for angle input for rotation matrix, enter 1 value per level in degrees, positive is counter-clockwise, no change when 0
    v = [0 0 0 0 0 0 0 -5 -10 5     5 -5 0 -10 0 0 0 -5 -10 5       0 -15 15 -10 0 0 0 -10 -10 -5];        
    vr = deg2rad(v); % Converts angle measure to radian measure

% Manual adjustment of scaling (Second transformation matrix, factor for x and y are the same to preserve proportionality)
    % Vector for scaling factor for x/y coordinates, one value per level, no change when 1
    scaling = [1 1 0.9 1.02 0.98 1.02 1 1.02 1 1     0.98 1.02 0.98 1 1.05 1.02 1 1.02 1 1   1 1 1 1.1 1.1 1.05 1 1.02 1 1.02];  

% Manual adjustment of XY position (Third transformation matrix, one value for x, one for y each)
    shift_x = [0 -150 0 -200 0 -280 50 200 150 -150    -280 -500 100 -250 -250 -350 -250 -150 -200 -350     -300 -150 -500 -300 -250 -200 -350 -250 -200 -500];
    shift_y = [0 -50 -150 -200 0 -120 -200 -200 -300 0     230 120 -150 100 100 -100 100 0 0 150      100 -80 300 -100 100 0 50 0 50 50];

% Plotting
figure(1)
clf
cmap = lines(30); % Select color map
hold on
AllDots = cell.empty; % Creates an empty cell where the new coordinates for each level are saved
for iFile = 1:numel(T)
    % Define rotation matrix:
    r = [cos(vr(iFile)) -sin(vr(iFile)); ...
         sin(vr(iFile)) cos(vr(iFile))];
    % Define scaling matrix:
    s = [scaling(iFile) 0; 0 scaling(iFile)];
    % Define shifting matrix:
    sh = [shift_x(iFile); shift_y(iFile)];
    % Calculate center of mass: (used to center all outlines)
    xm = mean(T_outlines{iFile}.X);
    ym = mean(T_outlines{iFile}.Y);
    % Plot outlines:
    x = T_outlines{iFile}.X-xm;
    y = -(T_outlines{iFile}.Y-ym);
    XY = s*(r*[x';y']+sh); % ' transposes rows to columns, first rotation then shift then scaling
    plot(XY(1,:), XY(2,:), 'Color', cmap(iFile,:))
    % Plot active brain cells:
    x = T{iFile}.X-xm;
    y = -(T{iFile}.Y-ym);
    XY = s*(r*[x';y']+sh); % First rotation, then shift, then scaling
    plot(XY(1,:), XY(2,:), '.', 'Color', cmap(iFile,:))
    AllDots{iFile} = XY; % Saves new coordinates per level in a cell
end
hold off
axis equal

%save(fullfile(cd,'AlignedData',sprintf('%s_%s_%s.mat',Lesion,Drug,Side)),'AllDots'); % Saves all adjusted coordinates as MATLAB file for further analysis.

% Plot heatmap for the whole group
figure(2)
clf
XY = horzcat(AllDots{:})'; % Combines all saved x and all saved y values from the different levels into one pool
h = histogram2(XY(:,1),XY(:,2),150,'DisplayStyle','tile','ShowEmptyBins','on','FaceColor','flat');
axis equal

disp('maximum is: '+ string(max(h.BinCounts, [], 'All')))
disp('minimum is: '+ string(min(h.BinCounts, [], 'All')))
