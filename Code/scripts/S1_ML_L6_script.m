load S1_data_tables_ML
% Are evoked spiking rates R_{M}, R_{ML} and R_{L} different?   Compare counts
% per trial per neuron for M and ML conditions.
for i=1:numel(these_cells)
    Cm=TBL{M}{i}.Ctest;   %mech trials
    Cml=TBL{ML}{i}.Ctest; %mech + light trials
    Cl=TBL{L}{i}.Ctest; %light trials
    % does ML differ significantly from M  evoked responses?  This tests per
    % unit
    %check for uneven trial length
    n=cellfun(@(x,y) min([numel(x) numel(y)]),Cm,Cml,'UniformOutput',false);
    Cm=cellfun(@(x,y) x(1:y),Cm,n,'UniformOutput',false);
    Cml=cellfun(@(x,y) x(1:y),Cml,n,'UniformOutput',false);
    P_MxML{i} = cell2mat(cellfun(@(x,y) signrank(x, y),Cm,Cml,'UniformOutput',false));
end

%% fraction of responsive S1 units for M and ML

INDICES={};
fig=figure;
t=tiledlayout(numel(these_cells),1,'TileSpacing','Compact');
%collect indices of sig. modulated for each condition.
for i=1:numel(these_cells)  %this is holdover loop so that it can be used for cortex too.

    Indices={};
    for J=1:3  %conditions
        tbl=Summary_stats{J};
        depths=TBL{J}{i}.depths;
        borderzone=[750 1050]*10; %for screening, not currently used
        inborder=find(depths>=borderzone(1) & depths<=borderzone(2));
        cellcount=tbl.cellcount(i)-numel(inborder);
        Indices{J}=setdiff(unique([find(TBL{J}{i}.dblZetaP<=.05);find(TBL{J}{i}.P_rs<=.05)]),inborder); %sig_changed
    end
    INDICES{i}=Indices;

end


%% Figure 4c, + Data for Venns in Figure S8e
VennData={};
Sig_any={};
f=figure;
t=tiledlayout(1,4,"TileSpacing",'compact');
for i=1:numel(these_cells)

    m=INDICES{i}{M};l=(INDICES{i}{L});ml=(INDICES{i}{ML});
    m_and_l=intersect(m,l);
    m_and_ml=intersect(m,ml);
    l_and_ml=intersect(l,ml);
    m_and_ml_and_l=intersect(m_and_ml,l_and_ml);  %only

    m_only=setdiff(setdiff(m,ml),l); % correct
    l_only=setdiff(setdiff(l,m),ml); % correct
    ml_only=setdiff(setdiff(ml,m),l);%%  %correct


    m_and_l_only=setdiff(m_and_l,m_and_ml_and_l);
    m_and_ml_only=setdiff(m_and_ml,m_and_ml_and_l);
    l_and_ml_only=setdiff(l_and_ml,m_and_ml_and_l);
    sig_any=unique([m; l; ml]);  %screen here for significance
    Sig_any{i}=sig_any;
    venndata=[numel(m_only)           
        numel(m_and_ml_only)   
        numel(ml_only)         
        numel(l_and_ml_only)         
        numel(l_only)               
        numel(m_and_l_only)
        numel(m_and_ml_and_l)];
    
    VennData{i}=venndata;

    cc=[[numel(find(TBL{L}{i}.dR(l)>0))  numel(find(TBL{L}{i}.dR(l)<0))  numel(find(TBL{L}{i}.dR(l)==0))];...
        [numel(find(TBL{M}{i}.dR(m)>0))  numel(find(TBL{M}{i}.dR(m)<0))  numel(find(TBL{M}{i}.dR(m)==0))];
        [numel(find(TBL{ML}{i}.dR(ml)>0))  numel(find(TBL{ML}{i}.dR(ml)<0))  numel(find(TBL{ML}{i}.dR(ml)==0))]];

    P=[[numel(find(TBL{L}{i}.dR(l)>0))  numel(find(TBL{L}{i}.dR(l)<0))  numel(find(TBL{L}{i}.dR(l)==0))]/numel(sig_any);...
        [numel(find(TBL{M}{i}.dR(m)>0))  numel(find(TBL{M}{i}.dR(m)<0))  numel(find(TBL{M}{i}.dR(m)==0))]/numel(sig_any);...
        [numel(find(TBL{ML}{i}.dR(ml)>0))  numel(find(TBL{ML}{i}.dR(ml)<0))  numel(find(TBL{ML}{i}.dR(ml)==0))]/numel(sig_any)];

    layer=tbl.these_cells(i);
    layer=cell2mat(string(layer));
    cc_new=[cc(:,1)+cc(:,3) cc(:,2)]';
    c=[[cellcount-sum(cc_new)];cc_new];


    nexttile
    bar([1:3],P(:,1));hold on; bar([1:3],-P(:,2));hold on
    ylim([-.8 1])
    box off
    layer=tbl.these_cells(i);
    layer=cell2mat(string(layer));
    title([layer ' n='  num2str(numel(sig_any))])
    %legend(['cat counts: ' [num2str([sum(cc) F(1) F(3)])]])
    %legend("Position",[0.5844,0.88905,0.33274,0.046032])
    ylabel 'response fraction'
    if i==numel(these_cells)
        set(gca,'xticklabel',{"L" 'M' 'ML'})
    end
    set(gca,'yminortick','on')

    %add text
    strs={};
    for j=1:3
        strs{j}=num2str(cc(j,1));
    end
    text([1:3]',P(:,1)+.025,strs')
    %add text
    strs={};
    for j=1:3
        strs{j}=num2str(cc(j,2));
    end
    text([1:3]',-P(:,2)-.025,strs')

    %add text
    strs={};
    for j=1:3
        strs{j}=num2str(cc(j,3));
    end
    text([1:3]',[0 0 0],strs','color',[1 1 1])

    %set(gcf,'position',[ 900.3333  479.0000  333.3333  420.0000])

end


% Supplementary Figure8e  Venn diagram figures (actual figure was made in Affinity Designer using the
%same data. Save VennData structure for breakdown by cell type
for J=1:numel(VennData)
    [~, fig]=vennX(VennData{J}, 0.1);
    title(these_cells(J))
end


%% Figure 4e modulation by depth

colors=[[0 0 0];[0 1 0];[1 0 0];[0 0 1]];

figure
scount=0;
for J=[L M ML]
    scount=scount+1;
    subplot(3,1,scount);
    Tbl=TBL{J};
    m=[];err=[];m_y=[];err_y=[];
    RATES={};

    for i=1:numel(Tbl)
        dilute=[.4 .4 .4];
        cl=colors(i,:);
        thislayer=Tbl{i};
        rng(i);
        thisjitter=randn(numel(thislayer.id),1)*8;
        sen=unique([find(thislayer.dblZetaP<=.05);find(thislayer.P_rs<=.05)]);
        xs=thislayer.depths(sen);
        ys=thislayer.mod_index(sen);
        r=TBL{1}{i}.m_rate_base(sen);
        xs=xs+thisjitter(sen);
        dilute=[.4 .4 .4];
        cl=cl+dilute;
        cl(find(cl>1))=1;

        s= stem(xs,ys,'Color',[.8 .8 .8],'marker','.','linewidth',.5,'markeredgecolor','none');

        hold on
        colormap hot
        brighten(-.25);
        scatter(xs,ys,r*5+.1,r*5+.1,'filled','markeredgecolor','k')
        RATES{i}=xs;
        colorbar
        caxis([0 40])
    end
    box off
    %set(gca,'xdir','reverse')
    %set(gca,'view',[90 -90])
    xlim([200 1300])
    ylim([-1.2 1.2])
    %xlabel depth
    ylabel MI
    if J<3,
        ax=gca;ax.XAxis.Visible = 'off';
    end
    ylabel([conds{J} ' MI'])
end



%%  Supplementary Figure 8e (scatter plots and dR histograms), compare MECH TO MECH LIGHT
dR={};
for i=1:numel(these_cells)
    thislayer=TBL{ML}{i};  %mech light only
    cl=M_color;
    pval=P_MxML{i};
    sig=find(pval<=.05); %filter for significantly different
    sen=unique([INDICES{i}{L}; INDICES{i}{ML}; INDICES{i}{M}]);
    sig=intersect(sig,sen);
    nsig=setdiff(sen,sig);
    %plot insen
    ys_s=TBL{ML}{i}.m_rate(sig);
    xs_s=TBL{M}{i}.m_rate(sig);
    ys_ns=TBL{ML}{i}.m_rate(nsig);
    xs_ns=TBL{M}{i}.m_rate(nsig);

    indices=find(ys_s<50); %exclude very rare outliers

    [S figs h dR{i}] = makePopulationScatter(xs_s(indices),ys_s(indices),xs_ns,ys_ns,colors(i,:));

    layer=thislayer.layers(i)
    layer=cell2mat(string(layer));
    figure(figs(1))
    set(gcf,'Position', [836.3333 815 344 312])
    n=numel(xs_s);
    title([ layer ' n=' num2str(n)])
    ylabel 'mech+L rate [Hz]'

    axis square
    box off

    figure(figs(2))
    subplot(1,2,1)
    title 'ML to M'
    xlim([.5 2.5])
    subplot(1,2,2)
    title '(ML-M) dR vs. M rate'
    close

    figure(figs(3))
    xlabel('FR')
    ylabel('Unit count')
    set(gcf,'Position', [1.2943e+03 983 255.3333 236.6667])
    box off

end

%% Fig. 4d  Modulation by layer and condition
violindata={};

for J=[L M ML]
    Tbl=TBL{J};
    for i=1:numel(Tbl)
        thislayer=Tbl{i};
        sen=Sig_any{i};
        ys=thislayer.mod_index(sen);
        violindata{J}{i}=ys;
    end
end

colors=[L_color;[.4 .4 .4];[0 .2 .7];];
fig=figure;
t=tiledlayout(1,numel(these_cells),'TileSpacing','Compact');
for i=1:numel(these_cells)
    nexttile
    plot([.5 3.5],[0 0],'k');hold on
    v=violinplot([violindata{L}{i} violindata{M}{i} violindata{ML}{i}],{'L','M','ML'},'ShowBox',true,'ShowData',true,'MarkerSize',10,'ShowWhiskers',false,'BoxColor',[0 0 0],'ViolinColor',colors);
    for j=1:numel(v)
        v(j).MedianPlot.SizeData=100;
        v(j).MedianPlot.MarkerFaceColor=[0 0 0];
    end
    ylim([-1.1 1.1])
    box off
    ylabel MI
    xlim([.5 3.5])
end


%%  Fig. 4f Scatter plot of MI by layer

PVALS={};
figure
t=tiledlayout(3,1,'TileSpacing','Compact');
for J=[L M ML]
    nexttile
    data=violindata{J};
    for i=1:4
        v=swarmchart(ones(size(data{i}))*i,data{i},'xjitterwidth',.75);
        v.SizeData=10;
        hold on
        q=quickstats(data{i});
        errs=[q.med-q.quart(1), q.quart(2)-q.med];
        errorbar(i,q.med,errs(1),errs(2),'Marker','.','LineWidth',2,'MarkerSize',30,'color','k')
        grid on
    end

    pvals=nan(4);
    for i=1:numel(data)
        for j=1:numel(data)
            [pvals(i,j),H] = ranksum(data{i},data{j});
        end
    end
    PVALS{J}=pvals;
    ylim([-1.2 1.2])
    xlim([0 5])
    set(gca,'xdir','reverse')
    set(gca,'view',[90 -90])
    box off
end



