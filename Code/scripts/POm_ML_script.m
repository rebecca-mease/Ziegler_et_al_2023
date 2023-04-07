load POm_data_tables_ML
%  Indices of reponsive units
Indices={};  %collect indices of sig. modulated for each condition.
for i=1%:numel(these_cells)  %this is holdover loop so that it can be used for cortex too.
    F=[];
    for J=[M L ML]  %conditions
        tbl=Summary_stats{J};
        cellcount=tbl.cellcount(i);
        F(J)=tbl.sig_either(i);   %cell counts per condition from summary table.
        %indices of those units
        Indices{J}=unique([find(TBL{J}{i}.dblZetaP<=.05);find(TBL{J}{i}.P_rs<=.05)]);

    end
end

% Find group membership and put in Venn diagram

m=Indices{M};l=(Indices{L});ml=(Indices{ML});
m_and_l=intersect(m,l);
m_and_ml=intersect(m,ml);
l_and_ml=intersect(l,ml);
m_and_ml_and_l=intersect(m_and_ml,l_and_ml);

m_only=setdiff(setdiff(m,ml),l);
l_only=setdiff(setdiff(l,m),ml);
ml_only=setdiff(setdiff(ml,m),l);%

m_and_l_only=setdiff(m_and_l,m_and_ml_and_l);
m_and_ml_only=setdiff(m_and_ml,m_and_ml_and_l);
l_and_ml_only=setdiff(l_and_ml,m_and_ml_and_l);
sig_any=unique([m; l; ml]);  %screen here for significantly modulated under any circumstance

venndata=[numel(m_only)
    numel(m_and_ml_only)
    numel(ml_only)
    numel(l_and_ml_only)
    numel(l_only)
    numel(m_and_l_only)
    numel(m_and_ml_and_l)];

[~, fig]=vennX(venndata, 0.1);
title(['Response categories VPL:   M:' num2str(numel(m)) ' L:' num2str(numel(l)) '  ML:' num2str(numel(ml))]);

cc=[[numel(find(TBL{M}{i}.dR(m)>0))  numel(find(TBL{M}{i}.dR(m)<0))  numel(find(TBL{M}{i}.dR(m)==0))];...
    [numel(find(TBL{L}{i}.dR(l)>0))  numel(find(TBL{L}{i}.dR(l)<0))  numel(find(TBL{L}{i}.dR(l)==0))];
    [numel(find(TBL{ML}{i}.dR(ml)>0))  numel(find(TBL{ML}{i}.dR(ml)<0))  numel(find(TBL{ML}{i}.dR(ml)==0))]];

P=[[numel(find(TBL{M}{i}.dR(m)>0))  numel(find(TBL{M}{i}.dR(m)<0))  numel(find(TBL{M}{i}.dR(m)==0))]/numel(sig_any);...
    [numel(find(TBL{L}{i}.dR(l)>0))  numel(find(TBL{L}{i}.dR(l)<0))  numel(find(TBL{L}{i}.dR(l)==0))]/numel(sig_any);...
    [numel(find(TBL{ML}{i}.dR(ml)>0))  numel(find(TBL{ML}{i}.dR(ml)<0))  numel(find(TBL{ML}{i}.dR(ml)==0))]/numel(sig_any)];

cc_new=[cc(:,1)+cc(:,3) cc(:,2)]';
c=[[cellcount-sum(cc_new)];cc_new];

fig=figure
bar([1:3],P(:,1));hold on; bar([1:3],-P(:,2));hold on
ylim([-.4 1])
box off
layer=tbl.these_cells(i);
layer=cell2mat(string(layer));
title([layer ' n='  num2str(cellcount)])
legend(['cat counts: ' [num2str([F(M) F(L) F(ML)])]]);
legend("Position",[0.5844,0.88905,0.33274,0.046032])
ylabel 'response fraction'
set(gca,'xticklabel',{conds{1} conds{2} conds{3}})
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
text([1:3]',[0 0 0],strs','color',[1 1 1]);
set(gcf,'position',[ 900.3333  479.0000  333.3333  420.0000]);


% Compare SPONTANEOUS TO LIGHT ONLY

thislayer=TBL{L}{i};
sen=l;
insen=setdiff(1:numel(thislayer.dblZetaP), sen);

%plot insen
ys_ns=thislayer.m_rate(insen);
xs_ns=thislayer.m_rate_base(insen);

%plot  sen
ys_s=thislayer.m_rate(sen);
xs_s=thislayer.m_rate_base(sen);
[S figs] = makePopulationScatter(xs_s,ys_s,xs_ns,ys_ns,L_color);

layer=thislayer.layers(i)
layer=cell2mat(string(layer));
figure(figs(1))
n=numel(xs_s);
title([ layer ' n=' num2str(n)]);
xlabel 'spontaneous rate [Hz]'
ylabel 'light-evoked rate [Hz]'
set(gca,'xscale','log')
set(gca,'yscale','log')
set(gca,'fontsize',14)
box off
axis square
box off
legend ns sig

legend("Position",[0.65181,0.41172,0.16422,0.11492]);

figure(figs(2))
subplot(1,2,1)
title 'spont to L'
xlim([.5 2.5])

subplot(1,2,2)
title '(M-spont) dR vs. spont rate'
close

figure(figs(3))
close


%Figure S12d Comparison of POm and VPL

pom=load('ModIndex_by_cond_pom.csv');
vpl=load('ModIndex_by_cond_vpl.csv');

vpl=vpl(:,2:end);
pom=pom(:,2:end);

figure
for J=1:3
    data={vpl(:,J) pom(:,J)};
    subplot(1,3,J)
    for i=1:2
        v=swarmchart(ones(size(data{i}))*i,data{i},'xjitterwidth',.5);
        v.SizeData=20;
        hold on;
        q=quickstats(data{i});
        errs=[q.med-q.quart(1), q.quart(2)-q.med];
        errorbar(i,q.med,errs(1),errs(2),'Marker','.','LineWidth',2,'MarkerSize',30,'color','k');
        ylim([-1.2 1.2]);
        xlim([.5 2.5]);
        set(gca,'xtick',[1 2],'xticklabel',{"VPL", "POm"});
    end
    [P,H] = ranksum(data{1},data{2},'tail','right');
    title(P);
end