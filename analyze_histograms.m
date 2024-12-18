% Load adjusted coordinates for active brain cells as prepared by
% align_sections.m:

DataDir = 'AlignedData'; % Define input folder, should contain one MATLAB 
                         % file per group and hemisphere as output from 
                         % preparation script ("lesiontype_treatment_hemisphere.mat")
Lesion = categorical(["OHDA";"SHAM"]);                  % Lesion types
Drug = categorical(["LD";"LD_ROP";"ROP";"saline"]);     % Treatments
Side = categorical(["L";"R"]);                          % Hemisphere
ix = fullfact([numel(Lesion) numel(Drug) numel(Side)]); % All combinations.
T = table(Lesion(ix(:,1)), Drug(ix(:,2)), Side(ix(:,3)), ...
    'VariableNames', {'Lesion' 'Drug' 'Side'});             % Table with all combinations.
fnames = compose("%s_%s_%s.mat", T.Lesion, T.Drug, T.Side); % All filenames.
sampleId = []; % Table with group IDs.
data = [];     % Cell array with XY coordinates.
for iFile = 1:numel(fnames)
    fname = fullfile(cd, DataDir, fnames(iFile));
    if exist(fname,"file")
        fprintf('Loading %s...', fname)
        L = load(fname);
        data = [data; L.AllDots']; % L.AllDots is a cell array where each 
                                   % cell contains data from one section.
        sampleId = [sampleId; repmat(T(iFile,:),numel(L.AllDots),1)];
        fprintf('\n')
    end
end

% Create 2D histograms:
res = 80;                   % Bin width in micrometers.
xedges = -3000:res:3000;    % Horizontal bin edges.
yedges = -3000:res:3000;    % Vertical bin edges.
rows = numel(xedges)-1;     % Number of horizontal bins.
columns = numel(yedges)-1;  % Number of vertical bins.
H = cellfun(@(x) histcounts2(x(1,:)',x(2,:)',xedges,yedges), data, ...
    'UniformOutput',false); % Create one histogram for each section.
Hv = cellfun(@(x) x(:), H, 'UniformOutput', false); % Vectorize the histograms...
Hv = horzcat(Hv{:});                                % and create a matrix...
Hv = Hv';                                           % where each row is one histogram.
fprintf('Number of bins with at least one cell count: %d\n', sum(sum(Hv)>0));

%% Perform PCA (SVD method):
[coeff,score,~,~,explained] = pca(Hv); % For 'coeff', Col1=PC1. For 'score', 
                                       % Col1=PC1 and Row1=Sample1.

% Save table with loadings for principal components:
T = [sampleId table(score(:,1), score(:,2), score(:,3),...
    'VariableNames', {'PC1' 'PC2' 'PC3'})];
writetable(T, 'PC_loadings.csv')

% Plot PCA
fname_cmap = fullfile(cd,'Colormap BlueRed.mat'); % Load blue-red colormap.
dummy = load(fname_cmap);
cmap = dummy.BlueRed;
h_fig = figure(1);
clf
tiledlayout(h_fig, 'flow')
CLim = 0.06; % Color limits.
for iPc = 1:20 % Plotting the first 20 PCs...
    nexttile
    imagesc(reshape(coeff(:,iPc),rows,columns),[-1 1]*CLim)
    colormap(cmap)
    title(sprintf('PC%d',iPc))
    axis square
end

% Plot heatmap of groupwise loadings for the first 20 PCs, sorted on PC1-PC3:
[G,TID] = findgroups(sampleId(:,{'Lesion' 'Drug'}));
Y = splitapply(@mean,score(:,1:20),G);
h_fig = figure(2);
clf
tiledlayout(h_fig,1,3)
for iPc = 1:3
    [~,ix]=sort(Y(:,iPc)); % Sort the heatmap on PC 1, 2, or 3.
    CLim = 30;
    nexttile;
    h_im = imagesc(Y(ix,:),[-1 1]*CLim);
    labels = compose('%s %s',string(TID.Lesion),string(TID.Drug));
    labels = strrep(labels,'_','+');
    h_im.Parent.YTickLabel = labels(ix);
    xlabel('PC')
    title(sprintf('Sorted on PC%d',iPc))
end

% Plot relative activation density images (Group2 minus Group1):
clear Cond
Cond{1} = {["SHAM" "saline"]; ... % Group 1
           ["SHAM" "LD"]};        % Group 2
Cond{2} = {["SHAM" "saline"]; ... % Group 1
           ["SHAM" "ROP"]};       % Group 2
Cond{3} = {["SHAM" "saline"]; ... % Group 1
           ["SHAM" "LD_ROP"]};    % Group 2
Cond{4} = {["SHAM" "LD"]; ...     % Group 1
           ["SHAM" "ROP"]};       % Group 2
Cond{5} = {["SHAM" "LD"]; ...     % Group 1
           ["SHAM" "LD_ROP"]};    % Group 2
Cond{6} = {["SHAM" "ROP"]; ...    % Group 1
           ["SHAM" "LD_ROP"]};    % Group 2
Cond{7} = {["OHDA" "saline"]; ... % Group 1
           ["OHDA" "LD"]};        % Group 2
Cond{8} = {["OHDA" "saline"]; ... % Group 1
           ["OHDA" "ROP"]};       % Group 2
Cond{9} = {["OHDA" "saline"]; ... % Group 1
           ["OHDA" "LD_ROP"]};    % Group 2
Cond{10} = {["OHDA" "LD"]; ...    % Group 1
            ["OHDA" "ROP"]};      % Group 2
Cond{11} = {["OHDA" "LD"]; ...    % Group 1
            ["OHDA" "LD_ROP"]};   % Group 2
Cond{12} = {["OHDA" "ROP"]; ...   % Group 1
            ["OHDA" "LD_ROP"]};   % Group 2

for iCond = 1:numel(Cond)   % Loop over pairwise comparisons.
    % Group 1:
    bSel1 = sampleId{:,1}==Cond{iCond}{1}(1) & ...
        sampleId{:,2}==Cond{iCond}{1}(2); % True for group 1.
    H1 = cat(3,H{bSel1});   % Concatenate all sections in one matrix.
    H1 = H1/(res/1000)^2;   % Convert to density (count/mm2).
    H1m = mean(H1,3);       % Average across sections.
    % Group 2:
    bSel2 = sampleId{:,1}==Cond{iCond}{2}(1) & ...
        sampleId{:,2}==Cond{iCond}{2}(2); % True for group 2.
    H2 = cat(3,H{bSel2});   % Concatenate all sections in one matrix.
    H2 = H2/(res/1000)^2;   % Convert to density (count/mm2).
    H2m = mean(H2,3);       % Average across sections.
    % Plot difference:
    figure;
    clf
    CLim = 270;
    imagesc(H2m-H1m, [-1 1]*CLim)
    axis square
    title(sprintf('%s,%s - %s,%s', Cond{iCond}{2}, Cond{iCond}{1}), 'Interpreter', 'none')
    colormap(cmap)
    colorbar
    % Statistical non-parametric mapping (Mann-Whitney U test between each
    % pixel pair):
    p = NaN(size(H1m));
    for iRow = 1:size(p,1)
        for iCol = 1:size(p,2)
            p(iRow,iCol) = ranksum(squeeze(H1(iRow,iCol,:)),squeeze(H2(iRow,iCol,:)));
        end
    end
    p((H1m-H2m)==0) = 1;                % Fix NaNs.
    cmap2 = summer;                     % Select colormap.
    cmap2(size(cmap2,1),:) = [1 1 1];   % Set the last color to white.
    % Plot statistical maps:
    figure
    clf
    imagesc(p, [0 0.01]) % Alpha defined as 0.01 according to Benjamini-Hochberg multiple comparison correction.
    axis square
    title(sprintf('%s,%s - %s,%s', Cond{iCond}{2}, Cond{iCond}{1}), 'Interpreter', 'none')
    colormap(cmap2)
    colorbar
end