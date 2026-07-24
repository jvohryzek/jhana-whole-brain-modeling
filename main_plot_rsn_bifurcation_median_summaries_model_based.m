%% Model-based median summaries
% Run this script to generate the four model-based median-summary figures.
% === Median per condition (Counting vs Formless Jhanas J7–J8 or J5–J8) across RSNs ===
repoRoot = fileparts(mfilename('fullpath'));
model = load(fullfile(repoRoot, 'data', 'model_bifurcation_parameters.mat'), 'data');
data = model.data;
clear model;
[~, ~, nrsn] = size(data);

rsnNames = {'Subcortical','Visual','Somato-Motor','Dorsal Attention', ...
            'Salience/Ventral Attention','Limbic','Control','Default'};

idxCounting = 2;
colBlue     = [0.20 0.40 0.85];   % Counting
colJhana    = [0.98 0.75 0.10];   % saturated yellow
colFormless = 0.6*colJhana + 0.4*[1 1 1]; % lighter yellow

%% Counting vs J7-J8
idxFormless = 9:10;

% --- Medians per RSN ---
medCounting  = nan(1,nrsn);
medFormless  = nan(1,nrsn);

for r = 1:nrsn
    % Counting: 40 values → median
    valsC = squeeze(data(:, idxCounting, r));
    medCounting(r) = median(valsC,'omitnan');

    % Formless: append all subj × J values (40×2=80 or 40×4=160)
    valsFL = reshape(data(:, idxFormless, r),[],1);
    medFormless(r) = median(valsFL,'omitnan');
end

% Arrange for grouped barplot: [Counting (left), J7-J8 (right)]
M = [medCounting(:), medFormless(:)];

% --- Figure ---
figure('Color','w','Position',[100 100 1200 520]);

% Shaded band between -0.02 and 0.02
yLimits   = [-1 1]*max(abs(M(:)))*1.1;
bandColor = [0.75 0.90 0.75];
hPatch = patch([0 nrsn+1 nrsn+1 0],[-0.02 -0.02 0.02 0.02],bandColor, ...
      'FaceAlpha',0.3,'EdgeColor','none'); hold on;
hTop = yline(0.02,'--','Color',bandColor*0.9,'LineWidth',1.5);
hBot = yline(-0.02,'--','Color',bandColor*0.9,'LineWidth',1.5);
uistack(hPatch,'bottom'); uistack(hTop,'bottom'); uistack(hBot,'bottom');

% Bars
b = bar(M,'grouped'); hold on;
for k = 1:2, b(k).FaceColor = 'flat'; b(k).EdgeColor = 'none'; end
for r = 1:nrsn
    b(1).CData(r,:) = colBlue;     % Counting
    b(2).CData(r,:) = colFormless; % Formless
end

% Guide line at zero
yline(0,'-','LineWidth',1.2,'Color',[0.5 0.5 0.5]);

% Axes, labels, legend
set(gca,'XTick',1:nrsn,'XTickLabel',rsnNames,'XTickLabelRotation',30, ...
        'FontSize',13,'LineWidth',1.2);
ylabel('Median bifurcation parameter','FontSize',14,'FontWeight','bold');
lgd = legend([b(1), b(2)],{'Counting','J7-J8'},'Location','northeast');
set(lgd,'FontSize',12,'Box','off');

% Symmetric y-limits
ylim(yLimits);

% Value labels on each bar
yr = diff(get(gca,'YLim')); dy = 0.015*yr;
for r = 1:nrsn
    xC = b(1).XEndPoints(r); vC = M(r,1);
    text(xC, vC + sign(vC)*dy, sprintf('%.4f', vC), ...
        'HorizontalAlignment','center','VerticalAlignment', ...
        ternary(vC>=0,'bottom','top'),'FontSize',11,'FontWeight','bold');
    xF = b(2).XEndPoints(r); vF = M(r,2);
    text(xF, vF + sign(vF)*dy, sprintf('%.4f', vF), ...
        'HorizontalAlignment','center','VerticalAlignment', ...
        ternary(vF>=0,'bottom','top'),'FontSize',11,'FontWeight','bold');
end

% Keep statistics in the script workspace
statsMedianJ78 = table((1:nrsn).', rsnNames(:), medCounting.', medFormless.', ...
    'VariableNames',{'RSN','Label','Median_Counting','Median_Formless'});

% Distances to critical point (absolute values)
distCounting = abs(medCounting(:));
distFL       = abs(medFormless(:));

rsnNames_v2 = {'Subc','VIS','SMN','DAN','VAN','LIM','CN','DMN'};

% Means per condition
m1 = mean(distCounting,'omitnan');
m2 = mean(distFL,'omitnan');

% Figure
figure('Color','w','Position',[100 100 700 520]); hold on;

% Bars
b = bar([1 2],[m1 m2],0.6,'FaceColor','flat','EdgeColor','none');
b.CData(1,:) = colBlue;
b.CData(2,:) = colFormless;

% Overlay RSN points + connecting lines
for r = 1:nrsn
    if isfinite(distCounting(r)) && isfinite(distFL(r))
        plot([1 2],[distCounting(r) distFL(r)],'-','Color',[0.6 0.6 0.6 0.5],'LineWidth',2);
        scatter(1,distCounting(r),40,'filled','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5]);
        scatter(2,distFL(r),40,'filled','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5]);
        text(1-0.05,distCounting(r),rsnNames_v2{r}, ...
             'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',10,'FontWeight','bold');
    end
end

% Axes
set(gca,'XTick',[1 2],'XTickLabel',{'Counting','J7-J8'},'FontSize',13,'LineWidth',1.2);
ylabel('Distance to bifurcation / critical point','FontSize',14,'FontWeight','bold');
xlim([0.5 2.5]); box on;

%%
% === Median per condition (Counting vs Formless Jhanas J5–J8) across RSNs ===
idxFormless = 7:10;

% --- Medians per RSN ---
medCounting  = nan(1,nrsn);
medFormless  = nan(1,nrsn);

for r = 1:nrsn
    % Counting: 40 values → median
    valsC = squeeze(data(:, idxCounting, r));
    medCounting(r) = median(valsC,'omitnan');

    % Formless: append all subj × J values (40×2=80 or 40×4=160)
    valsFL = reshape(data(:, idxFormless, r),[],1);
    medFormless(r) = median(valsFL,'omitnan');
end

% Arrange for grouped barplot: [Counting (left), Formless (right)]
M = [medCounting(:), medFormless(:)];

% --- Figure ---
figure('Color','w','Position',[100 100 1200 520]);

% Shaded band between -0.02 and 0.02
yLimits   = [-1 1]*max(abs(M(:)))*1.1;
bandColor = [0.75 0.90 0.75];
hPatch = patch([0 nrsn+1 nrsn+1 0],[-0.02 -0.02 0.02 0.02],bandColor, ...
      'FaceAlpha',0.3,'EdgeColor','none'); hold on;
hTop = yline(0.02,'--','Color',bandColor*0.9,'LineWidth',1.5);
hBot = yline(-0.02,'--','Color',bandColor*0.9,'LineWidth',1.5);
uistack(hPatch,'bottom'); uistack(hTop,'bottom'); uistack(hBot,'bottom');

% Bars
b = bar(M,'grouped'); hold on;
for k = 1:2, b(k).FaceColor = 'flat'; b(k).EdgeColor = 'none'; end
for r = 1:nrsn
    b(1).CData(r,:) = colBlue;     % Counting
    b(2).CData(r,:) = colFormless; % Formless
end

% Guide line at zero
yline(0,'-','LineWidth',1.2,'Color',[0.5 0.5 0.5]);

% Axes, labels, legend
set(gca,'XTick',1:nrsn,'XTickLabel',rsnNames,'XTickLabelRotation',30, ...
        'FontSize',13,'LineWidth',1.2);
ylabel('Median bifurcation parameter','FontSize',14,'FontWeight','bold');
lgd = legend([b(1), b(2)],{'Counting','J5-J8'},'Location','northeast');
set(lgd,'FontSize',12,'Box','off');

% Symmetric y-limits
ylim(yLimits);

% Value labels on each bar
yr = diff(get(gca,'YLim')); dy = 0.015*yr;
for r = 1:nrsn
    xC = b(1).XEndPoints(r); vC = M(r,1);
    text(xC, vC + sign(vC)*dy, sprintf('%.4f', vC), ...
        'HorizontalAlignment','center','VerticalAlignment', ...
        ternary(vC>=0,'bottom','top'),'FontSize',11,'FontWeight','bold');
    xF = b(2).XEndPoints(r); vF = M(r,2);
    text(xF, vF + sign(vF)*dy, sprintf('%.4f', vF), ...
        'HorizontalAlignment','center','VerticalAlignment', ...
        ternary(vF>=0,'bottom','top'),'FontSize',11,'FontWeight','bold');
end

% Keep statistics in the script workspace
statsMedianJ58 = table((1:nrsn).', rsnNames(:), medCounting.', medFormless.', ...
    'VariableNames',{'RSN','Label','Median_Counting','Median_J5-J8'});

% Distances to critical point (absolute values)
distCounting = abs(medCounting(:));
distFL       = abs(medFormless(:));

rsnNames_v2 = {'Subc','VIS','SMN','DAN','VAN','LIM','CN','DMN'};

% Means per condition
m1 = mean(distCounting,'omitnan');
m2 = mean(distFL,'omitnan');

% Figure
figure('Color','w','Position',[100 100 700 520]); hold on;

% Bars
b = bar([1 2],[m1 m2],0.6,'FaceColor','flat','EdgeColor','none');
b.CData(1,:) = colBlue;
b.CData(2,:) = colFormless;

% Overlay RSN points + connecting lines
for r = 1:nrsn
    if isfinite(distCounting(r)) && isfinite(distFL(r))
        plot([1 2],[distCounting(r) distFL(r)],'-','Color',[0.6 0.6 0.6 0.5],'LineWidth',2);
        scatter(1,distCounting(r),40,'filled','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5]);
        scatter(2,distFL(r),40,'filled','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5]);
        text(1-0.05,distCounting(r),rsnNames_v2{r}, ...
             'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',10,'FontWeight','bold');
    end
end

% Axes
set(gca,'XTick',[1 2],'XTickLabel',{'Counting','J5-J8'},'FontSize',13,'LineWidth',1.2);
ylabel('Distance to bifurcation / critical point','FontSize',14,'FontWeight','bold');
xlim([0.5 2.5]); box on;

%% --- small helper (inline ternary for text alignment)
function out = ternary(cond,a,b)
    if cond, out=a; else, out=b; end
end
