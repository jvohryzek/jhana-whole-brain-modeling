%% Model-based bifurcation distributions
% Run this script to generate the three model-based distribution figures.
repoRoot = fileparts(mfilename('fullpath'));
addpath(fullfile(repoRoot, 'third_party', 'violinplot'));
model = load(fullfile(repoRoot, 'data', 'model_bifurcation_parameters.mat'), 'data');
data = model.data;
clear model;

[~, ~, nrsn] = size(data);
rsnNames = {'Subcortical','Visual','Somato-Motor','Dorsal Attention', ...
            'Salience/Ventral Attention','Limbic','Control','Default'};

idxCounting = 2;
idxForm = 3:6;       % J1-J4
idxFormless = 7:10;  % J5-J8
idxJ78 = 9:10;       % J7-J8

colBlue       = [0.20 0.40 0.85];   % Counting
colJhanaBase  = [0.98 0.75 0.10];   % base yellow
intensityJ8   = 1 - (8-1)*(0.65/(8-1)); % = 0.35
colJ8         = colJhanaBase*intensityJ8 + (1-intensityJ8)*[1 1 1];
fade = @(k) 1 - (k-1)*(0.65/(8-1));
intForm = mean(fade(1:4));
intFormless = mean(fade(5:8));
colForm = colJhanaBase*intForm + (1-intForm)*[1 1 1];
colFormless = colJhanaBase*intFormless + (1-intFormless)*[1 1 1];

%% Counting vs J7-J8 distributions
condNamesPair = {'Counting','J7–J8'};
figure('Color','w','Position',[100 100 1400 800]);

nCounting = zeros(nrsn, 1);
nJ78 = zeros(nrsn, 1);
cohensD = nan(nrsn, 1);

for r = 1:nrsn
    subplot(2,4,r); hold on;

    % === Gather values ===
    y1 = squeeze(data(:, idxCounting, r));      % 40×1
    y2 = reshape(data(:, idxJ78, r), [], 1);    % 80×1

    % --- Half violins ---
    [f1,xi1] = ksdensity(y1,'Bandwidth',[]);
    f1 = f1 / max(f1) * 0.3;
    patch(1 - f1, xi1, colBlue, 'FaceAlpha',0.3,'EdgeColor','none');
    [f2,xi2] = ksdensity(y2,'Bandwidth',[]);
    f2 = f2 / max(f2) * 0.3;
    patch(2 + f2, xi2, colJ8, 'FaceAlpha',0.3,'EdgeColor','none');

    % --- Jittered scatter points ---
    scatter(1 + (rand(size(y1))-0.5)*0.1, y1, 18, ...
        'MarkerFaceColor',colBlue, 'MarkerEdgeColor','none', 'MarkerFaceAlpha',0.6);
    scatter(2 + (rand(size(y2))-0.5)*0.1, y2, 18, ...
        'MarkerFaceColor',colJ8, 'MarkerEdgeColor','none', 'MarkerFaceAlpha',0.6);

    % --- Medians ---
    m1 = median(y1,'omitnan');
    m2 = median(y2,'omitnan');
    plot([0.75 1.25], [m1 m1], '-', 'Color',[0 0 0], 'LineWidth',2);
    plot([1.75 2.25], [m2 m2], '-', 'Color',[0 0 0], 'LineWidth',2);

    % Axis labels
    title(rsnNames{r}, 'FontWeight', 'bold');
    ylabel('Bifurcation parameter');
    xlim([0.3 2.7]); box on;
    set(gca,'XTick',[1 2],'XTickLabel',condNamesPair);

    ymax = 0.1;
    ymin = -0.1;
    ylim([ymin,ymax]);

    % Horizontal line y=0
    yline(0, '-', 'LineWidth', 1, 'Color', [0.6 0.6 0.6]);

    % === Cohen's d (two-sample, pooled SD) ===
    n1 = numel(y1); n2 = numel(y2);
    m1m = mean(y1,'omitnan');
    m2m = mean(y2,'omitnan');
    s1 = var(y1,'omitnan'); s2 = var(y2,'omitnan');
    spooled = sqrt(((n1-1)*s1 + (n2-1)*s2) / (n1+n2-2));
    d = (m1m - m2m) / spooled;

    % Effect size stars
    if abs(d) < 0.2
        stars = 'n.s.';
    elseif abs(d) < 0.5
        stars = '*';
    elseif abs(d) < 0.8
        stars = '**';
    else
        stars = '***';
    end

    txt = sprintf('%s  (d=%.2f)', stars, d);
    text(1.5, 0.08, txt, ...
        'HorizontalAlignment','center', 'VerticalAlignment','bottom', 'FontWeight','bold');

    nCounting(r) = n1;
    nJ78(r) = n2;
    cohensD(r) = d;
end

statsPair = table((1:nrsn).', rsnNames(:), nCounting, nJ78, cohensD, ...
    'VariableNames', {'RSN','LabelRSN','Ncount','Nj7j8','CohensD'});
disp('=== Effect sizes (Cohen''s d): Counting vs J7–J8 per RSN ===');
disp(statsPair);

%% Distance to the critical point
% 1) Median per RSN & condition (omit NaNs)
medC  = nan(1,nrsn);
medF  = nan(1,nrsn);
medFL = nan(1,nrsn);
for r = 1:nrsn
    yC  = squeeze(data(:, idxCounting,  r));
    yF  = squeeze(mean(data(:, idxForm,     r), 2, 'omitnan'));   % subj mean J1–J4
    yFL = squeeze(mean(data(:, idxFormless, r), 2, 'omitnan'));   % subj mean J5–J8
    medC(r)  = median(yC,  'omitnan');
    medF(r)  = median(yF,  'omitnan');
    medFL(r) = median(yFL, 'omitnan');
end

% 2) Distances to critical point (= |median|)
distC  = abs(medC(:));
distF  = abs(medF(:));
distFL = abs(medFL(:));

% 3) Violin plot: each violin has 7 RSN points
figure('Color','w','Position',[100 100 800 540]);
vp = violinplot([distC distF distFL], {'Counting','Form','Formless'}, ...
                'ShowMean', false, 'ShowData', false);
vp(1).ViolinColor = colBlue;
vp(2).ViolinColor = colForm;
vp(3).ViolinColor = colFormless;
for c = 1:3
    vp(c).ViolinAlpha = 1.0;
    if isprop(vp(c),'OutlineColor'), vp(c).OutlineColor = 'none'; end
    if isprop(vp(c),'EdgeColor'),    vp(c).EdgeColor    = 'none'; end
end
hold on;

% 4) Overlay RSN points (light grey, jittered)
randomStream = RandStream('mt19937ar', 'Seed', 0);
jit = 0.08;
x1 = 1 + jit*randn(randomStream, nrsn, 1);
x2 = 2 + jit*randn(randomStream, nrsn, 1);
x3 = 3 + jit*randn(randomStream, nrsn, 1);
scatter(x1, distC,  32, 'filled', 'MarkerFaceColor',[0.3 0.3 0.3], 'MarkerFaceAlpha',0.35, 'MarkerEdgeAlpha',0.35);
scatter(x2, distF,  32, 'filled', 'MarkerFaceColor',[0.3 0.3 0.3], 'MarkerFaceAlpha',0.35, 'MarkerEdgeAlpha',0.35);
scatter(x3, distFL, 32, 'filled', 'MarkerFaceColor',[0.3 0.3 0.3], 'MarkerFaceAlpha',0.35, 'MarkerEdgeAlpha',0.35);

% 5) Thick mean lines per group
mC  = mean(distC,'omitnan');
mF  = mean(distF,'omitnan');
mFL = mean(distFL,'omitnan');
plot([0.75 1.25], [mC  mC ], 'k-', 'LineWidth', 2);
plot([1.75 2.25], [mF  mF ], 'k-', 'LineWidth', 2);
plot([2.75 3.25], [mFL mFL], 'k-', 'LineWidth', 2);

% Axes & labels
set(gca,'FontSize',13,'LineWidth',1.2);
ylabel('Distance to bifurcation / critical point','FontSize',14,'FontWeight','bold');
xlim([0.5 3.5]); box on;

% y-limits with headroom for three brackets
yy = [distC; distF; distFL]; if isempty(yy), yy = 0; end
ymin = min(yy); ymax = max(yy); margin = 0.12*(ymax - ymin + eps);
ylim([max(0, ymin - margin), ymax + 3*margin]);

% 6) Wilcoxon signed-rank tests (paired across RSNs)
pairMask12 = isfinite(distC)  & isfinite(distF);
pairMask13 = isfinite(distC)  & isfinite(distFL);
pairMask23 = isfinite(distF)  & isfinite(distFL);

[p12,h12,s12] = signrank(distC(pairMask12),  distF(pairMask12),  'method','approximate');
[p13,h13,s13] = signrank(distC(pairMask13),  distFL(pairMask13), 'method','approximate');
[p23,h23,s23] = signrank(distF(pairMask23),  distFL(pairMask23), 'method','approximate');

W12 = s12.signedrank; z12 = isfield(s12,'zval')*s12.zval;
W13 = s13.signedrank; z13 = isfield(s13,'zval')*s13.zval;
W23 = s23.signedrank; z23 = isfield(s23,'zval')*s23.zval;

% Brackets + p-value labels
yl = ylim; base = yl(2) - 2.2*margin; step = 0.8*margin; tick = 0.18*margin;
plot([1 1 2 2],[base base+tick base+tick base],'k-','LineWidth',1.2);
text(1.5, base+tick+0.05*margin, sprintf('signrank p=%.3g (C vs Form)', p12), ...
    'HorizontalAlignment','center','VerticalAlignment','bottom','FontWeight','bold','FontSize',12);

plot([1 1 3 3],[base+step base+step+tick base+step+tick base+step],'k-','LineWidth',1.2);
text(2.0, base+step+tick+0.05*margin, sprintf('signrank p=%.3g (C vs Formless)', p13), ...
    'HorizontalAlignment','center','VerticalAlignment','bottom','FontWeight','bold','FontSize',12);

plot([2 2 3 3],[base+2*step base+2*step+tick base+2*step+tick base+2*step],'k-','LineWidth',1.2);
text(2.5, base+2*step+tick+0.05*margin, sprintf('signrank p=%.3g (Form vs Formless)', p23), ...
    'HorizontalAlignment','center','VerticalAlignment','bottom','FontWeight','bold','FontSize',12);

% 7) Keep statistics in the script workspace
statsSignrank_3groups = table( ...
    {'C vs Form'; 'C vs Formless'; 'Form vs Formless'}, ...
    [p12; p13; p23], [h12; h13; h23], ...
    [W12; W13; W23], [z12; z13; z23], ...
    'VariableNames', {'Comparison','p','h','W_signedrank','zval'});

%%
% === Half-violins with jitter: Counting vs Form (J1–J4) vs Formless (J5–J8), per RSN ===
groupNames   = {'Counting','Form','Formless'};

figure('Color','w','Position',[50 50 1500 800]);

medianCounting = nan(nrsn, 1);
medianForm = nan(nrsn, 1);
medianFormless = nan(nrsn, 1);
cohensDCountingForm = nan(nrsn, 1);
cohensDFormFormless = nan(nrsn, 1);

for r = 1:nrsn
    subplot(2,4,r); hold on;

    % === Values (append, not avg) ===
    yC  = squeeze(data(:, idxCounting, r));   % 40×1
    yF  = reshape(data(:, idxForm, r), [], 1);      % 160 values
    yFL = reshape(data(:, idxFormless, r), [], 1);  % 160 values

    % --- Half violins ---
    maxw = 0.3; % max width
    % Counting
    [f1,xi1] = ksdensity(yC); f1 = f1/max(f1)*maxw;
    patch(1 - f1, xi1, colBlue, 'FaceAlpha',0.3,'EdgeColor','none');
    % Form
    [f2,xi2] = ksdensity(yF); f2 = f2/max(f2)*maxw;
    patch(2 + f2, xi2, colForm, 'FaceAlpha',0.3,'EdgeColor','none');
    % Formless
    [f3,xi3] = ksdensity(yFL); f3 = f3/max(f3)*maxw;
    patch(3 - f3, xi3, colFormless, 'FaceAlpha',0.3,'EdgeColor','none');

    % --- Jittered points ---
    scatter(1+(rand(size(yC))-0.5)*0.1, yC, 18, 'MarkerFaceColor',colBlue, ...
        'MarkerEdgeColor','none','MarkerFaceAlpha',0.6);
    scatter(2+(rand(size(yF))-0.5)*0.1, yF, 18, 'MarkerFaceColor',colForm, ...
        'MarkerEdgeColor','none','MarkerFaceAlpha',0.6);
    scatter(3+(rand(size(yFL))-0.5)*0.1, yFL, 18, 'MarkerFaceColor',colFormless, ...
        'MarkerEdgeColor','none','MarkerFaceAlpha',0.6);

    % --- Medians ---
    med = [median(yC,'omitnan'), median(yF,'omitnan'), median(yFL,'omitnan')];
    plot([0.75 1.25],[med(1) med(1)],'k-','LineWidth',2);
    plot([1.75 2.25],[med(2) med(2)],'k-','LineWidth',2);
    plot([2.75 3.25],[med(3) med(3)],'k-','LineWidth',2);
    plot(1:3,med,'k-','LineWidth',1.2);

    % Axis
    title(rsnNames{r},'FontWeight','bold');
    ylabel('Bifurcation parameter');
    xlim([0.5 3.5]); set(gca,'XTick',1:3,'XTickLabel',groupNames); box on;
    ylim([-0.1 0.1]);  % fixed y-limits
    yline(0,'-','Color',[0.6 0.6 0.6]);

    % === Cohen’s d (two-sample, pooled SD) ===
    cohd = @(a,b) (mean(a,'omitnan')-mean(b,'omitnan'))/ ...
                  sqrt(((numel(a)-1)*var(a,'omitnan') + (numel(b)-1)*var(b,'omitnan'))/(numel(a)+numel(b)-2));
    d_C_F  = cohd(yC,yF);
    d_F_FL = cohd(yF,yFL);

    % --- Stars at y=0.08 (only C–Form and Form–Formless) ---
    text(1.5, 0.08, sprintf('%s', starfun(d_C_F)), ...
        'HorizontalAlignment','center','FontWeight','bold');
    text(2.5, 0.08, sprintf('%s', starfun(d_F_FL)), ...
        'HorizontalAlignment','center','FontWeight','bold');

    medianCounting(r) = med(1);
    medianForm(r) = med(2);
    medianFormless(r) = med(3);
    cohensDCountingForm(r) = d_C_F;
    cohensDFormFormless(r) = d_F_FL;
end

statsTbl = table((1:nrsn).', rsnNames(:), medianCounting, medianForm, ...
    medianFormless, cohensDCountingForm, cohensDFormFormless, ...
    'VariableNames', {'RSN','LabelRSN','Median_Counting','Median_Form', ...
    'Median_Formless','CohensD_C_Form','CohensD_Form_Formless'});

% --- helper for stars
function s = starfun(d)
    if abs(d) < 0.2
        s = 'n.s.';
    elseif abs(d) < 0.5
        s = '*';
    elseif abs(d) < 0.8
        s = '**';
    else
        s = '***';
    end
end
