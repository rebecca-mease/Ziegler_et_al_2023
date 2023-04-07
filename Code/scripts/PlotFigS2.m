%% PlotFigS2


% Primary Somatosensory Cortex Bidirectionally Modulates Sensory Gain and
% Nociceptive Behavior in a Layer-Specific Manner - Ziegler et al.,2023

% Plots Optotagging data to generate figures for Supplementary Figure 2a
% and b

%% Prepping Data

% Loading depths from directory
tbl = [];
tbl{1,1} = 'L5'; tbl{1,2} = readtable(fullfile(basedir, filesep, 'L5_Table.csv'));
tbl{2,1} = 'L6'; tbl{2,2} = readtable(fullfile(basedir, filesep, 'L6_Table.csv'));

l5_sd = 3.5;
l6_sd = 2;
sd_cutoff = [l5_sd, l6_sd];


% Iterate through tables to plot data
for ct = 1 : size(tbl,1)
    current_tbl = tbl{ct,2};
    depth = current_tbl.depth;
    mn = current_tbl.mn; % mean latency
    sd = current_tbl.sd; % standard deviation
    IDs = current_tbl.IDs; % unit IDs
    fd = current_tbl.fd; % Response fidelity - not used here. Calculated probability of spiking within 25ms per laser pulse.
    spontFR = current_tbl.spontFR; % spontaneous firing rate (outside of trigger response window)

    % Removing and highlighting artefact units from plot
    lowestLat = 2;
    lowestSD = 0.1;

    mn(mn<lowestLat) = NaN;
    sd(sd<lowestSD) = NaN;
    naN = isnan(mn) | isnan(sd);

    mn(naN) = NaN;
    sd(naN) = NaN;



    % Eliminating taggedFS units from plots
    fsTagged = current_tbl.fsTagged == 1;
    sd(fsTagged) = NaN;
    mn(fsTagged) = NaN;
    naN = isnan(mn) | isnan(sd);

    tagged = current_tbl.isTagged == 1;
    TaggedIDs = IDs(tagged);
    nontagged = ~tagged & ~naN;
    minDepth = min(depth(tagged));
    maxDepth = max(depth(tagged));

    name = tbl{ct,1};
    pulseWidth = 10;
   


   

 %% Scatter Plots


    fig = figure('Name', [name, '_scatter'], 'Color', 'White');


    scatter(mn(~tagged),sd(~tagged),30,depth(~tagged), 'o');
    hold on
    scatter(mn(tagged),sd(tagged),45,depth(tagged),'o','filled'); 
    hold off

    hot = hot;
    n_colors = 256;
    blue_colormap = [zeros(n_colors,2), linspace(1,0,n_colors)'];
    red = hot(10:240,:);
    blue = blue_colormap(10:140,:);
    cmap = [blue; red];
    c_axis_limits = [-1400, 0];
    colormap(cmap);
    caxis(c_axis_limits);
    colorbar;
    ylim([0, 10]);
    xlim([2, 17]);
    xlabel('First-Spike Latency [ms]');
    ylabel('StdDev [ms]', 'Interpreter', 'tex')
    ax = gca;
    hold on
    laser = plot([0:pulseWidth], (ax.YLim(2))*ones(pulseWidth + 1,1), 'Color', [0 1 1 0.5], 'LineWidth', 3);
    ax = gca;
    colorbar(ax);
    ax.FontName = 'Arial';
    ax.FontSize = 20;
    ax.XTick = [2:2:16];
    ax.YTick = [1:9];
    lgd = legend;
    lgd.String{1} = ['Non-tagged Units (n=', num2str(sum(nontagged)), ')'];
    lgd.String{2} = ['Opto-tagged Units (n=', num2str(sum(tagged)), ')'];
    lgd.String{3} = ['Laser Pulse (',num2str(pulseWidth), 'ms)'];
    rectangle(ax, 'Position', [ax.XLim(1), ax.YLim(1), pulseWidth-ax.XLim(1), ax.YLim(2) - ax.YLim(1)], ...
        'LineStyle', 'none', 'FaceColor', [0 1 1 0.1])
    lgd.FontSize = 10;


    %%
end
