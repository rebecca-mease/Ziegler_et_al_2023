%% PlotFigS2c


% Primary Somatosensory Cortex Bidirectionally Modulates Sensory Gain and
% Nociceptive Behavior in a Layer-Specific Manner - Ziegler et al.,2023

% Plots Optotagging data to generate figures for Supplementary Figure 2c
%%
load(fullfile(basedir, filesep, 'L5_Depths.mat'));
load(fullfile(basedir, filesep, 'L6_Depths.mat'));
L5ind = zeros(length(L5_Depths),1);
L6ind = ones(length(L6_Depths),1);
figure('Color', 'white');
boxplot([L5_Depths; L6_Depths], [L5ind; L6ind], 'Notch', 'on', 'Labels', {'L5-ChR2', 'L6-ChR2'}, 'Colors','rb', 'PlotStyle','traditional')
title('Depth Distributions of Opto-tagged Units');
ax = gca;
ax.FontName = 'Arial';
ax.FontSize = 15;
box off

