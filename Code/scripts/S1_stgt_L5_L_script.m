load S1_data_tables_stgt_L5_L

%% depth resolved modulation
tbl=S1_L{1};
these_exp=unique(tbl.exp);

fig=figure
for i=1:numel(these_exp)
    tbl=S1_L{1};
    tbl=tbl(find(tbl.exp==these_exp(i)),:);
    sen=unique([find(tbl.dblZetaP<=.05);find(tbl.P_rs<=.05)]);
    tbl=tbl(sen,:);
    depths=tbl.depths;
    r=tbl.m_rate_base;
    mi=tbl.mod_index;

    %plot mech sen
    %[S figs] = makePopulationScatter(xs_s,ys_s,xs_s,ys_s,[.1 .1 .1])
    thisjitter=randn(numel(depths),1)*5;
    depths=depths+thisjitter;
    %stem(depths,mi,'o','linewidth',1)

    stem(depths,mi,'Color',[.8 .8 .8],'marker','.','linewidth',1,'markeredgecolor','none')
    hold on
    colormap hot
    brighten(-.25);
    scatter(depths,mi,r*4+.1,r*4+.1,'filled','markeredgecolor','k')

    hold on
end

set(gca,'xdir','reverse')
set(gca,'view',[90 -90])

xlim([200 1300])
ylim([-1.2 1.2])

%xlim([200 1300])
%set(gca,'ydir','reverse');
title 'L5 stgt inhibition'
ylabel MI
xlabel depth
colorbar
caxis([0 40])


%%  violin of MI by layer
violindata={};
for J=[1]
    Tbl=TBL{J}(2:end);%ignore L1 neurons
    m=[];err=[];m_y=[];err_y=[];
    RATES={}
    for i=1:numel(Tbl)
        thislayer=Tbl{i};
        sen=1:size(thislayer,1);  %default all neurons
        sen=unique([find(thislayer.dblZetaP<=.05);find(thislayer.P_rs<=.05)]);
        ys=thislayer.mod_index(sen);
        violindata{J}{i}=ys;
    end
end

%%

fig=figure
t=tiledlayout(1,numel(these_cells),'TileSpacing','Compact');
data=violindata{1};

for i=1:4
v=swarmchart(ones(size(data{i}))*i,data{i},'xjitterwidth',.5)
v.SizeData=20
hold on
q=quickstats(data{i})
errs=[q.med-q.quart(1), q.quart(2)-q.med];
errorbar(i,q.med,errs(1),errs(2),'Marker','.','LineWidth',2,'MarkerSize',30,'color','k')
end

set(gca,'view',[90 -90])

%xlim([200 1300])
set(gca,'xdir','reverse')
ylim([-1.2 1.2])