%% Model-free global-connectivity distributions
% Run this script to generate the two model-free connectivity figures.
%
% The input is the released subject-level global-connectivity table. This
% script reproduces the Python analysis in MATLAB: condition aggregation,
% complete-case selection, Tukey 1.5-IQR trimming, paired t-tests, exact
% Wilcoxon signed-rank tests, Cohen's dz, BH/FDR correction, and plotting.

repoRoot = fileparts(mfilename('fullpath'));
addpath(fullfile(repoRoot, 'third_party', 'violinplot'));
subjectData = readtable(fullfile(repoRoot, 'data', ...
    'subject_global_connectivity_aligned.csv'), 'TextType', 'string');

groupNames = [ ...
    "Subcortical"
    "Visual"
    "Somato-Motor"
    "Dorsal Attention"
    "Salience/Ventral Attention"
    "Limbic"
    "Control"
    "Default"
];

[pairData, pairStats] = buildCountingJ78(subjectData, groupNames);
[triadData, triadStats] = buildCountingEarlyLate(subjectData, groupNames);

pairFigure = plotCountingJ78(pairData, pairStats, groupNames);
triadFigure = plotCountingEarlyLate(triadData, triadStats, groupNames);


function [plotData, statsTable] = buildCountingJ78(subjectData, groupNames)
nGroups = numel(groupNames);
plotData = cell(nGroups, 1);
stats = repmat(emptyPairStats(), nGroups, 1);

for groupIdx = 1:nGroups
    group = groupNames(groupIdx);
    groupData = subjectData(subjectData.group == group, :);
    subjects = sort(unique(groupData.subject));

    rows = strings(0, 1);
    counting = zeros(0, 1);
    meanJ78 = zeros(0, 1);
    for subject = subjects.'
        subjectRows = groupData(groupData.subject == subject, :);
        countingValues = subjectRows.value(subjectRows.task == "counting");
        j78Values = subjectRows.value(ismember(subjectRows.task, ["j7", "j8"]));
        if numel(countingValues) == 1 && ~isempty(j78Values)
            rows(end + 1, 1) = subject; %#ok<AGROW>
            counting(end + 1, 1) = countingValues; %#ok<AGROW>
            meanJ78(end + 1, 1) = mean(j78Values, 'omitnan'); %#ok<AGROW>
        end
    end

    remove = tukeyOutlierMask(counting) | tukeyOutlierMask(meanJ78);
    rows(remove) = [];
    counting(remove) = [];
    meanJ78(remove) = [];

    plotData{groupIdx} = table(rows, counting, meanJ78, ...
        'VariableNames', {'subject', 'counting', 'mean_j7_j8'});
    stats(groupIdx) = pairedStats(group, "counting", "j7_j8_mean", ...
        counting, meanJ78, sum(remove));
end

pValues = [stats.p_value].';
wValues = [stats.p_wilcoxon].';
qValues = fdrBH(pValues);
qWValues = fdrBH(wValues);
for groupIdx = 1:nGroups
    stats(groupIdx).q_value = qValues(groupIdx);
    stats(groupIdx).q_wilcoxon = qWValues(groupIdx);
end
statsTable = struct2table(stats);
end


function [plotData, statsTable] = buildCountingEarlyLate(subjectData, groupNames)
nGroups = numel(groupNames);
plotData = cell(nGroups, 1);
stats = repmat(emptyPairStats(), nGroups * 2, 1);

for groupIdx = 1:nGroups
    group = groupNames(groupIdx);
    groupData = subjectData(subjectData.group == group, :);
    subjects = sort(unique(groupData.subject));

    rows = strings(0, 1);
    counting = zeros(0, 1);
    early = zeros(0, 1);
    late = zeros(0, 1);
    for subject = subjects.'
        subjectRows = groupData(groupData.subject == subject, :);
        countingValues = subjectRows.value(subjectRows.task == "counting");
        earlyValues = subjectRows.value(ismember(subjectRows.task, ["j1", "j2", "j3", "j4"]));
        lateValues = subjectRows.value(ismember(subjectRows.task, ["j5", "j6", "j7", "j8"]));
        if numel(countingValues) == 1 && ~isempty(earlyValues) && ~isempty(lateValues)
            rows(end + 1, 1) = subject; %#ok<AGROW>
            counting(end + 1, 1) = countingValues; %#ok<AGROW>
            early(end + 1, 1) = mean(earlyValues, 'omitnan'); %#ok<AGROW>
            late(end + 1, 1) = mean(lateValues, 'omitnan'); %#ok<AGROW>
        end
    end

    remove = tukeyOutlierMask(counting) | ...
        tukeyOutlierMask(early) | tukeyOutlierMask(late);
    rows(remove) = [];
    counting(remove) = [];
    early(remove) = [];
    late(remove) = [];

    plotData{groupIdx} = table(rows, counting, early, late, ...
        'VariableNames', {'subject', 'counting', 'early_jhana', 'late_jhana'});

    firstIdx = groupIdx;
    secondIdx = nGroups + groupIdx;
    stats(firstIdx) = pairedStats(group, "counting", "early_jhana", ...
        counting, early, sum(remove));
    stats(secondIdx) = pairedStats(group, "early_jhana", "late_jhana", ...
        early, late, sum(remove));
end

for contrastIdx = 1:2
    indices = (contrastIdx - 1) * nGroups + (1:nGroups);
    qValues = fdrBH([stats(indices).p_value].');
    qWValues = fdrBH([stats(indices).p_wilcoxon].');
    for localIdx = 1:nGroups
        stats(indices(localIdx)).q_value = qValues(localIdx);
        stats(indices(localIdx)).q_wilcoxon = qWValues(localIdx);
    end
end
statsTable = struct2table(stats);
end


function stats = pairedStats(group, taskA, taskB, x, y, nRemoved)
stats = emptyPairStats();
stats.group = group;
stats.task_a = taskA;
stats.task_b = taskB;
stats.n_subjects = numel(x);
stats.n_removed_outliers = nRemoved;
stats.mean_a = mean(x);
stats.mean_b = mean(y);
diffValues = y - x;
stats.mean_diff_b_minus_a = mean(diffValues);
stats.sd_diff = std(diffValues, 0);
stats.t_value = mean(diffValues) / (stats.sd_diff / sqrt(numel(diffValues)));
stats.p_value = 2 * tcdf(-abs(stats.t_value), numel(diffValues) - 1);
stats.cohen_dz = stats.mean_diff_b_minus_a / stats.sd_diff;
[stats.wilcoxon_statistic, stats.p_wilcoxon] = exactWilcoxon(diffValues);
end


function stats = emptyPairStats()
stats = struct( ...
    'group', "", ...
    'task_a', "", ...
    'task_b', "", ...
    'n_subjects', 0, ...
    'n_removed_outliers', 0, ...
    'mean_a', NaN, ...
    'mean_b', NaN, ...
    'mean_diff_b_minus_a', NaN, ...
    'sd_diff', NaN, ...
    't_value', NaN, ...
    'p_value', NaN, ...
    'q_value', NaN, ...
    'cohen_dz', NaN, ...
    'wilcoxon_statistic', NaN, ...
    'p_wilcoxon', NaN, ...
    'q_wilcoxon', NaN);
end


function mask = tukeyOutlierMask(values)
values = values(:);
if numel(values) < 4
    mask = false(size(values));
    return;
end
sortedValues = sort(values);
q1 = numpyQuantile(sortedValues, 0.25);
q3 = numpyQuantile(sortedValues, 0.75);
iqrValue = q3 - q1;
if ismembertol(iqrValue, 0, 1e-12)
    mask = false(size(values));
    return;
end
lower = q1 - 1.5 * iqrValue;
upper = q3 + 1.5 * iqrValue;
mask = values < lower | values > upper;
end


function value = numpyQuantile(sortedValues, probability)
position = (numel(sortedValues) - 1) * probability;
lowerIdx = floor(position) + 1;
upperIdx = ceil(position) + 1;
weight = position - floor(position);
value = sortedValues(lowerIdx) * (1 - weight) + sortedValues(upperIdx) * weight;
end


function [statistic, pValue] = exactWilcoxon(diffValues)
diffValues = diffValues(:);
diffValues = diffValues(abs(diffValues) > 1e-14);
if isempty(diffValues)
    statistic = 0;
    pValue = 1;
    return;
end

ranks = tiedrank(abs(diffValues));
positiveRank = sum(ranks(diffValues > 0));
negativeRank = sum(ranks(diffValues < 0));
statistic = min(positiveRank, negativeRank);

scaledRanks = round(2 * ranks);
totalRank = sum(scaledRanks);
counts = zeros(1, totalRank + 1);
counts(1) = 1;
for rankValue = scaledRanks.'
    counts(rankValue + 1:end) = ...
        counts(rankValue + 1:end) + counts(1:end - rankValue);
end
scaledStatistic = round(2 * statistic);
oneSided = sum(counts(1:scaledStatistic + 1)) / (2 ^ numel(ranks));
pValue = min(1, 2 * oneSided);
end


function adjusted = fdrBH(pValues)
pValues = pValues(:);
adjusted = nan(size(pValues));
valid = isfinite(pValues);
if ~any(valid)
    return;
end
[sortedP, order] = sort(pValues(valid));
nValues = numel(sortedP);
sortedAdjusted = sortedP .* nValues ./ (1:nValues).';
for idx = nValues - 1:-1:1
    sortedAdjusted(idx) = min(sortedAdjusted(idx), sortedAdjusted(idx + 1));
end
sortedAdjusted = min(max(sortedAdjusted, 0), 1);
restored = zeros(size(sortedAdjusted));
restored(order) = sortedAdjusted;
adjusted(valid) = restored;
end


function fig = plotCountingJ78(plotData, statsTable, groupNames)
countingColor = [0.20 0.40 0.85];
j78Color = [0.98 0.75 0.10];

fig = figure('Color', 'w', 'Units', 'inches', ...
    'Position', [0.5 0.5 13.7136 7.63], ...
    'Renderer', 'painters');
layout = tiledlayout(fig, 2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

for groupIdx = 1:numel(groupNames)
    ax = nexttile(layout);
    hold(ax, 'on');
    values = plotData{groupIdx};
    series = [values.counting, values.mean_j7_j8];
    drawViolinSeries(ax, series, {countingColor, j78Color}, [-1 1]);
    drawPairedSubjects(ax, series, {countingColor, j78Color});
    drawMedians(ax, series, false);

    stats = statsTable(statsTable.group == groupNames(groupIdx), :);
    annotation = sprintf([ ...
        't = %.2f\np_t = %.3f\nq_t = %.3f\ndz = %.2f\n' ...
        'W = %.1f\np_w = %.3f\nq_w = %.3f\nrm = %d'], ...
        stats.t_value, stats.p_value, stats.q_value, stats.cohen_dz, ...
        stats.wilcoxon_statistic, stats.p_wilcoxon, stats.q_wilcoxon, ...
        stats.n_removed_outliers);

    stylePanel(ax, groupNames(groupIdx), series, 0.22);
    xticks(ax, [1 2]);
    xticklabels(ax, {'Counting', 'Mean J7-J8'});
    if ismember(groupIdx, [1 5])
        ylabel(ax, 'Global network connectivity (r)', 'FontSize', 13);
    end
    text(ax, 0.03, 0.97, annotation, 'Units', 'normalized', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
        'FontSize', 8.5, 'BackgroundColor', 'w', 'Margin', 2, ...
        'Interpreter', 'none');
end

end


function fig = plotCountingEarlyLate(plotData, statsTable, groupNames)
countingColor = [0.20 0.40 0.85];
earlyColor = [0.96 0.61 0.12];
lateColor = [0.98 0.75 0.10];
colors = {countingColor, earlyColor, lateColor};

fig = figure('Color', 'w', 'Units', 'inches', ...
    'Position', [0.5 0.5 14.5136 8.00], ...
    'Renderer', 'painters');
layout = tiledlayout(fig, 2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

for groupIdx = 1:numel(groupNames)
    ax = nexttile(layout);
    hold(ax, 'on');
    values = plotData{groupIdx};
    series = [values.counting, values.early_jhana, values.late_jhana];
    drawViolinSeries(ax, series, colors, [-1 1 -1]);
    drawPairedSubjects(ax, series, colors);
    drawMedians(ax, series, true);

    groupStats = statsTable(statsTable.group == groupNames(groupIdx), :);
    first = groupStats(groupStats.task_a == "counting", :);
    second = groupStats(groupStats.task_a == "early_jhana", :);
    annotation = sprintf([ ...
        'C-E: t=%.2f, p_t=%.3f, q_t=%.3f, dz=%.2f\n' ...
        '     W=%.1f, p_w=%.3f, q_w=%.3f\n' ...
        'E-L: t=%.2f, p_t=%.3f, q_t=%.3f, dz=%.2f\n' ...
        '     W=%.1f, p_w=%.3f, q_w=%.3f\nrm = %d'], ...
        first.t_value, first.p_value, first.q_value, first.cohen_dz, ...
        first.wilcoxon_statistic, first.p_wilcoxon, first.q_wilcoxon, ...
        second.t_value, second.p_value, second.q_value, second.cohen_dz, ...
        second.wilcoxon_statistic, second.p_wilcoxon, second.q_wilcoxon, ...
        first.n_removed_outliers);

    stylePanel(ax, groupNames(groupIdx), series, 0.24);
    xticks(ax, [1 2 3]);
    xticklabels(ax, {'Counting', 'Early Jhana\newline(J1-J4)', ...
        'Late Jhana\newline(J5-J8)'});
    if ismember(groupIdx, [1 5])
        ylabel(ax, 'Global network connectivity (r)', 'FontSize', 13);
    end
    text(ax, 0.03, 0.97, annotation, 'Units', 'normalized', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
        'FontSize', 7.1, 'BackgroundColor', 'w', 'Margin', 2, ...
        'Interpreter', 'none');
end

end


function drawViolinSeries(ax, series, colors, sides)
labels = arrayfun(@(idx) sprintf('%d', idx), 1:size(series, 2), ...
    'UniformOutput', false);
violins = violinplot(series, labels, ...
    'ShowData', false, 'ShowMean', false, 'ShowNotches', false, ...
    'Width', 0.30);

for idx = 1:numel(violins)
    violins(idx).ViolinColor = colors{idx};
    violins(idx).ViolinAlpha = 0.22;
    if isgraphics(violins(idx).ViolinPlot)
        violins(idx).ViolinPlot.EdgeColor = 'none';
        xData = violins(idx).ViolinPlot.XData;
        if sides(idx) < 0
            xData(xData > idx) = idx;
        else
            xData(xData < idx) = idx;
        end
        violins(idx).ViolinPlot.XData = xData;
    end
    hideGraphic(violins(idx).BoxPlot);
    hideGraphic(violins(idx).WhiskerPlot);
    hideGraphic(violins(idx).MedianPlot);
    hideGraphic(violins(idx).NotchPlots);
    hideGraphic(violins(idx).MeanPlot);
    hideGraphic(violins(idx).ScatterPlot);
end
hold(ax, 'on');
end


function hideGraphic(graphicHandle)
if ~isempty(graphicHandle)
    valid = isgraphics(graphicHandle);
    set(graphicHandle(valid), 'Visible', 'off');
end
end


function drawPairedSubjects(ax, series, colors)
randomStream = RandStream('mt19937ar', 'Seed', 42);
jitter = (rand(randomStream, size(series, 1), 1) - 0.5) * 0.08;
xPositions = (1:size(series, 2)) + jitter;
for seriesIdx = 1:size(series, 2)
    scatter(ax, xPositions(:, seriesIdx), series(:, seriesIdx), 26, ...
        'MarkerFaceColor', colors{seriesIdx}, ...
        'MarkerEdgeColor', colors{seriesIdx} * 0.82, ...
        'MarkerFaceAlpha', 0.50, 'MarkerEdgeAlpha', 0.80);
end
end


function drawMedians(ax, series, connectMedians)
medians = median(series, 1, 'omitnan');
if connectMedians
    plot(ax, 1:numel(medians), medians, '-', ...
        'Color', [0.30 0.30 0.30], 'LineWidth', 1.2);
end
for idx = 1:numel(medians)
    plot(ax, [idx - 0.20 idx + 0.20], [medians(idx) medians(idx)], ...
        'k-', 'LineWidth', 2.2);
end
end


function stylePanel(ax, groupName, series, padFraction)
allValues = series(:);
yMin = min(allValues);
yMax = max(allValues);
if yMax > yMin
    pad = (yMax - yMin) * padFraction;
else
    pad = 0.05;
end
ylim(ax, [yMin - pad * 0.25, yMax + pad]);
xlim(ax, [0.55, size(series, 2) + 0.45]);
title(ax, groupName, 'FontSize', 16, 'FontWeight', 'bold');
set(ax, 'FontSize', 11.5, 'LineWidth', 1.0, ...
    'XTickLabelRotation', 0, 'Layer', 'top', 'Box', 'on');
grid(ax, 'on');
ax.XGrid = 'off';
ax.GridAlpha = 0.14;
ax.TickDir = 'out';
end
