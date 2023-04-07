%% PlotFig1B


% Primary Somatosensory Cortex Bidirectionally Modulates Sensory Gain and
% Nociceptive Behavior in a Layer-Specific Manner - Ziegler et al.,2023

% Plots Optotagging data to generate figures for Figure 1b

%% Prepping Data
% Loading data tables from directories
tbl = [];
tbl{1,1} = 'L6'; tbl{1,2} = readtable(fullfile(basedir, filesep, 'L6_Table.csv'));

l6_sd = 2;
sd_cutoff = [l6_sd];


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
   


   

    %% Depth vs Latency Plots


    fig = figure('Name', [name, '_depthlatency'], 'Color', 'White');

    
    errorbar(mn(tagged),depth(tagged),sd(tagged), 'horizontal', 'LineStyle', 'none', 'Marker', 'd', 'Color',[0, 0.5, 1], 'MarkerSize', 2.5, 'LineWidth', 0.01, 'CapSize', 0, 'MarkerFaceColor',[0,0.5,1]);

    ylim([-1400, 0]);
    yticks([round(round(min(depth),-2)/2,-2)*2-200:200:0]);
    xlim([0 25]);
    xticks(0:5:25);
    xlabel('First-Spike Latency [ms]');
    ylabel('Unit Depth [\mum]', 'Interpreter', 'tex')
    ax = gca;
    hold on
    laser = plot([0:pulseWidth], ax.YLim(2)*ones(pulseWidth + 1,1), 'Color', [0 1 1 0.5], 'LineWidth', 3);
    ax = gca;
    ax.FontName = 'Arial';
    ax.FontSize = 10;
    lgd = legend;
    lgd.String{1} = (['Unit Mean +/- SD (n=', num2str(sum(nontagged)), ')']);
    lgd.String{2} = (['Opto-tagged Units (n=', num2str(sum(tagged)), ')']);
    lgd.String{3} = ['Laser Pulse (',num2str(pulseWidth), 'ms)'];

    rectangle(ax, 'Position', [ax.XLim(1), ax.YLim(1), pulseWidth, ax.YLim(2) - ax.YLim(1)], ...
        'LineStyle', 'none', 'FaceColor', [0 1 1 0.1])
    lgd.FontSize = 10;


    %%
end
