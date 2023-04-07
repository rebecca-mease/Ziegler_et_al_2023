load VPL_data_tables_ML M_color ML_color L_color M L ML Summary_stats TBL conds these_cells;

close all
Indices={};  %collect indices of sig. modulated for each condition.
F=[];
i=1;
for J=1:3  % stimulus conditions
    tbl=Summary_stats{J};
    cellcount=tbl.cellcount(1);
    F(J)=tbl.sig_either(1);   %cell counts per condition from summary table.
    %indices of those units
    Indices{J}=unique([find(TBL{J}{1}.dblZetaP<=.05);find(TBL{J}{i}.P_rs<=.05)]);
end

m=Indices{M};l=Indices{L};ml=Indices{ML};
m_and_l=intersect(m,l);
m_and_ml=intersect(m,ml);
l_and_ml=intersect(l,ml);
m_and_ml_and_l=intersect(m_and_ml,l_and_ml);

m_only=setdiff(setdiff(m,ml),l);
l_only=setdiff(setdiff(l,m),ml);
ml_only=setdiff(setdiff(ml,m),l);

m_and_l_only=setdiff(m_and_l,m_and_ml_and_l);
m_and_ml_only=setdiff(m_and_ml,m_and_ml_and_l);
l_and_ml_only=setdiff(l_and_ml,m_and_ml_and_l);

sig_any=unique([m; l; ml]);  %screen here for significantly modulated under any circumstance

venndata=[numel(m_only)            % m_only
    numel(m_and_ml_only)    %
    numel(ml_only)         % ml_only
    numel(l_and_ml_only)          %
    numel(l_only)                % l_only
    numel(m_and_l_only)
    numel(m_and_ml_and_l)];   %correct %


% FIGURE 3B  Venn diagram figure (actual figure was made in Affinity Designer using the
%same data.
[~, fig]=vennX(venndata, 0.1);
title(['Response categories VPL:   M:' num2str(numel(m)) ' L:' num2str(numel(l)) '  ML:' num2str(numel(ml))]);

i=1;

cc=[[numel(find(TBL{L}{i}.dR(l)>0))  numel(find(TBL{L}{i}.dR(l)<0))  numel(find(TBL{L}{i}.dR(l)==0))];...
    [numel(find(TBL{M}{i}.dR(m)>0))  numel(find(TBL{M}{i}.dR(m)<0))  numel(find(TBL{M}{i}.dR(m)==0))];
    [numel(find(TBL{ML}{i}.dR(ml)>0))  numel(find(TBL{ML}{i}.dR(ml)<0))  numel(find(TBL{ML}{i}.dR(ml)==0))]];

P=[[numel(find(TBL{L}{i}.dR(l)>0))  numel(find(TBL{L}{i}.dR(l)<0))  numel(find(TBL{L}{i}.dR(l)==0))]/numel(sig_any);...
    [numel(find(TBL{M}{i}.dR(m)>0))  numel(find(TBL{M}{i}.dR(m)<0))  numel(find(TBL{M}{i}.dR(m)==0))]/numel(sig_any);...
    [numel(find(TBL{ML}{i}.dR(ml)>0))  numel(find(TBL{ML}{i}.dR(ml)<0))  numel(find(TBL{ML}{i}.dR(ml)==0))]/numel(sig_any)];

% FIGURE 3C
fig=figure;
bar([1:3],P(:,1));hold on; bar([1:3],-P(:,2));hold on
ylim([-.2 1])
box off
layer=tbl.these_cells(i);
layer=cell2mat(string(layer));
title([layer ' n='  num2str(cellcount)]);
legend(['cat counts: ' [num2str([F(L) F(M) F(ML)])]]);
legend("Position",[0.5844,0.88905,0.33274,0.046032]);
ylabel 'response fraction'
set(gca,'xticklabel',{conds{L} conds{M} conds{ML}});
set(gca,'yminortick','on');

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
text([1:3]',-P(:,2)-.025,strs');

%add text
strs={};
for j=1:3
    strs{j}=num2str(cc(j,3));
end
text([1:3]',[0 0 0],strs','color',[1 1 1]);
set(gcf,'position',[ 900.3333  479.0000  333.3333  420.0000]);


% Supplementary Figure8 C1 Compare SPONTANEOUS TO LIGHT ONLY


thislayer=TBL{L}{i};
sen=l;
insen=setdiff(1:numel(thislayer.dblZetaP), sen);

%plot insen
ys_ns=thislayer.m_rate(insen);
xs_ns=thislayer.m_rate_base(insen);

%plot mech sen
ys_s=thislayer.m_rate(sen);
xs_s=thislayer.m_rate_base(sen);
[S figs] = makePopulationScatter(xs_s,ys_s,xs_ns,ys_ns,L_color);

layer=thislayer.layers(i);
layer=cell2mat(string(layer));
figure(figs(1))
n=numel(xs_s);
title([ layer ' n=' num2str(n)])
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

% Supplementary Figure8 C2    %Compare SPONTANEOUS TO MECH ONLY

mech_sen=m;
mech_insen=setdiff(1:numel(thislayer.dblZetaP), mech_sen);

%plot insen
ys_ns=thislayer.m_rate(mech_insen);
xs_ns=thislayer.m_rate_base(mech_insen);

%plot mech sen
ys_s=thislayer.m_rate(mech_sen);
xs_s=thislayer.m_rate_base(mech_sen);
[S figs] = makePopulationScatter(xs_s,ys_s,xs_ns,ys_ns,M_color);

layer=thislayer.layers(i);
layer=cell2mat(string(layer));
figure(figs(1));
n=numel(xs_s);
title([ layer ' n=' num2str(n)]);
xlabel 'spontaneous rate [Hz]'
ylabel 'mech-evoked rate [Hz]'
set(gca,'xscale','log')
set(gca,'yscale','log')
set(gca,'fontsize',14)
axis square
box off
legend ns sig

legend("Position",[0.68746,0.485,0.16422,0.11492]);
figure(figs(2))
subplot(1,2,1)
title 'spont to M'
xlim([.5 2.5])

subplot(1,2,2)
title '(M-spont) dR vs. spont rate'
figure(figs(3))
close
figure(figs(2))
close

% Supplementary Figure8 C3    %Compare SPONTANEOUS TO ML

thislayer=TBL{L}{i};
sen=ml;
insen=setdiff(1:numel(thislayer.dblZetaP), sen);

%plot insen
ys_ns=thislayer.m_rate(insen);
xs_ns=thislayer.m_rate_base(insen);

%plot mech sen
ys_s=thislayer.m_rate(sen);
xs_s=thislayer.m_rate_base(sen);
[S figs] = makePopulationScatter(xs_s,ys_s,xs_ns,ys_ns,ML_color);

layer=thislayer.layers(i);
layer=cell2mat(string(layer));
figure(figs(1))
n=numel(xs_s);
title([ layer ' n=' num2str(n)]);
xlabel 'spontaneous rate [Hz]'
ylabel 'ML-evoked rate [Hz]'
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
figure(figs(3))
close
figure(figs(2))
close

% Figure 3e Compare MECH TO MECH LIGHT

i=1;
thislayer=TBL{L}{i};
%sen=unique([m; ml])
sen=sig_any;
insen=setdiff(1:numel(thislayer.dblZetaP), sen);

ys_s=TBL{ML}{i}.m_rate(sen);
xs_s=TBL{M}{i}.m_rate(sen);
light_only=TBL{L}{i}.m_rate(sen);
[S figs] = makePopulationScatter(xs_s,ys_s,[],[],ML_color);

layer=thislayer.layers(i);
layer=cell2mat(string(layer));
figure(figs(1))
n=numel(xs_s);
title([ layer ' n=' num2str(n)])
ylabel 'mech+L rate [Hz]'
xlabel 'mech-evoked rate [Hz]'
set(gca,'xscale','log')
set(gca,'yscale','log')
axis square
box off
set(gca,'fontsize',14)
%legend ns sig

figure(figs(2))
subplot(1,2,1)
title 'M to ML'
xlim([.5 2.5])

subplot(1,2,2)
title '(ML-M) dR vs. M rate'
figure(figs(3))
close
figure(figs(2))
close


%  Figure 3h (burst probability)    Change in burst probability between M and ML

sen=sig_any;
C_m=TBL{M}{i}.EventcTest;
C_ml=TBL{ML}{i}.EventcTest;
C_l=TBL{L}{i}.EventcTest;

counts_m=cellfun(@(x) numel(x),C_m,'UniformOutput',true);
bc_m=cellfun(@(x) numel(find(x>1)),C_m,'UniformOutput',true);
singlesm=counts_m-bc_m;

counts_ml=cellfun(@(x) numel(x),C_ml,'UniformOutput',true);
bc_ml=cellfun(@(x) numel(find(x>1)),C_ml,'UniformOutput',true);
singlesml=counts_ml-bc_ml;

counts_l=cellfun(@(x) numel(x),C_l,'UniformOutput',true);
bc_l=cellfun(@(x) numel(find(x>1)),C_l,'UniformOutput',true);
singlesl=counts_l-bc_l;

pval_m_ml_b=[];
for j=1:numel(C_m)
    ctable=[[bc_m(j) singlesm(j)];[bc_ml(j) singlesml(j)]];
    pval_m_ml_b(j)=mcnemar(ctable);
end

sen=intersect(sig_any,find(pval_m_ml_b<=1));
insen=setdiff(sig_any, sen);

BP_m=cellfun(@(x) numel(find(x>1))/numel(x),C_m,'UniformOutput',true);
BP_ml=cellfun(@(x) numel(find(x>1))/numel(x),C_ml,'UniformOutput',true);
BP_l=cellfun(@(x) numel(find(x>1))/numel(x),C_l,'UniformOutput',true);
[S figs] = makePopulationScatter(BP_m(sen),BP_ml(sen),BP_m(insen),BP_ml(insen),'k');
figure(figs(1))
box off
axis square
xlabel BP_m
ylabel BP_{ml}
ylim([0 1])
xlim([0 1])
% set(gca,'yscale','log')
% set(gca,'xscale','log')

figure(figs(2))
subplot(1,2,1)
ylabel 'BP'
title 'M vs. ML'
xlim([.5 2.5])

subplot(1,2,2)
title '(ML-M) vs. M'
ylabel '\Delta BP'
xlabel BP_m

figure(figs(3))
close
figure(figs(2))
close

% Supplementary Figure 8d  (change in evoked modulation index vs. mech response, colored by bursting

thislayer=TBL{M}{i};
sen=unique([ml;m])
sen=sig_any;
%sen=m;
insen=setdiff(1:numel(thislayer.dblZetaP), sen);


% xlim([.5 2.5])

sen=sig_any;
%sen=m;
insen=setdiff(1:numel(thislayer.dblZetaP), sen);

% %plot mech sen units for ML and M conditions
ys_s=TBL{ML}{i}.mod_index(sen);
xs_s=TBL{M}{i}.mod_index(sen);
light_only_mi=TBL{L}{i}.mod_index(sen);

% fig=figure
% colors=[L_color; [.4 .4 .4];[0 .2 .7];]
% v=violinplot([light_only_mi xs_s ys_s ],{'L','M','ML'},'ShowBox',true,'ShowData',true,'MarkerSize',10,'ShowWhiskers',false,'BoxColor',[0 0 0],'ViolinColor',colors);
% for j=1:numel(v)
%     v(j).MedianPlot.SizeData=100;
%     v(j).MedianPlot.MarkerFaceColor=[0 0 0];
% end
% ylabel 'Modulation Index'
% title(S.stats_x_y.p_sr)
% box off
% set(gca,'fontsize',14)
% ylim([-1 1])
% xlim([0 4])
% plot([0 4],[0 0])

%set(gcf,'Position', [814.3333  352.3333  420.6667  406.0000])

%
sen=sig_any;
C_m=TBL{M}{i}.EventcTest(sen);
C_ml=TBL{ML}{i}.EventcTest(sen);
BP_m=cellfun(@(x) numel(find(x>1))/numel(x),C_m,'UniformOutput',true);
BP_ml=cellfun(@(x) numel(find(x>1))/numel(x),C_ml,'UniformOutput',true);

[BP_m thisorder]=sort(BP_m);
fig=figure;
colormap hot;

s=scatter(xs_s(thisorder),ys_s(thisorder)-xs_s(thisorder),30,BP_m,'filled',...
    'MarkerEdgeColor',[1 1 1]*0,'MarkerFaceAlpha',1);
colorbar
xlabel 'MI_{Mech}'
ylabel '\Delta Mod Index'
title('colormap= M burst probability')
set(gca,'fontsize',14)
title('colormap= M burst probability')

%

% getting ISIs for Fig.3 h

sen=sig_any;
ISIm=TBL{M}{i}.ISIs_test;
ISIml=TBL{ML}{i}.ISIs_test;

p10m=TBL{M}{i}.ISI_Ptest;
p10ml=TBL{ML}{i}.ISI_Ptest;

%when does fraction of small isis significantly change?
isis10m=TBL{M}{i}.ISIs_test;  %all isis
indices_m=TBL{M}{i}.ISI_indices_test; %small isis indices only
isis10ml=TBL{M}{i}.ISIs_test;  %all isis
indices_ml=TBL{M}{i}.ISIs_test;  %small isis indices only

% find number of all isis and number of small isis
n_m=cellfun(@numel,isis10m,'UniformOutput',false);
n_ml=cellfun(@numel,isis10ml,'UniformOutput',false);
x_m=cellfun(@numel,indices_m,'UniformOutput',false);
x_ml=cellfun(@numel,indices_ml,'UniformOutput',false);


sen=sig_any;
C_m=TBL{M}{i}.EventcTest;
C_ml=TBL{ML}{i}.EventcTest;
C_l=TBL{L}{i}.EventcTest;

counts_m=cellfun(@(x) numel(x),C_m,'UniformOutput',true);
bc_m=cellfun(@(x) numel(find(x>1)),C_m,'UniformOutput',true);
singlesm=counts_m-bc_m;

counts_ml=cellfun(@(x) numel(x),C_ml,'UniformOutput',true);
bc_ml=cellfun(@(x) numel(find(x>1)),C_ml,'UniformOutput',true);
singlesml=counts_ml-bc_ml;

counts_l=cellfun(@(x) numel(x),C_l,'UniformOutput',true);
bc_l=cellfun(@(x) numel(find(x>1)),C_l,'UniformOutput',true);
singlesl=counts_l-bc_l;
% 
% pval_m_ml_b=[];
% for j=1:numel(C_m)
%     ctable=[[bc_m(j) singlesm(j)];[bc_ml(j) singlesml(j)]];
%     pval_m_ml_b(j)=mcnemar(ctable);
% end
% 
% 
% %sen=intersect(sig_any,find(pval_m_ml_b<=1));
sen=sig_any;
insen=setdiff(sig_any, sen);

BP_m=cellfun(@(x) numel(find(x>1))/numel(x),C_m,'UniformOutput',true);
BP_ml=cellfun(@(x) numel(find(x>1))/numel(x),C_ml,'UniformOutput',true);
BP_l=cellfun(@(x) numel(find(x>1))/numel(x),C_l,'UniformOutput',true);

% Figure 3i   plotting difference in ISIs for high BP units

%sen=sig_any; % use sen from above!
hiBP=find(BP_m>.1); these_units=intersect(hiBP,sen);
lowBP=setdiff(sen,hiBP);
BW=1; %binwidth in ms
[n edges]=histcounts([cell2mat(ISIm(these_units)); cell2mat(ISIml(these_units))],'BinWidth',BW);
H_m=cellfun(@(x) histcounts(x,edges,'normalization','cdf'), ISIm(these_units),'UniformOutput',false);
H_ml=cellfun(@(x) histcounts(x,edges,'normalization','cdf'), ISIml(these_units),'UniformOutput',false);
isi_t=edges(1:end-1)+BW/2;
ids=TBL{i}{M}.id(these_units);

%
fig=figure;
semilogx(isi_t',median(cell2mat(H_m)),'color',M_color,'linewidth',2)
hold on
semilogx(isi_t',median(cell2mat(H_ml)),'color',ML_color,'linewidth',2)
xlim([0 10000])
xlabel 'Log_{10} ISI [ms]'
ylabel CDF(ISI)

hold on
BW=1;
[n edges]=histcounts([cell2mat(ISIm(lowBP)); cell2mat(ISIml(lowBP))],'BinWidth',BW);
H_m=cellfun(@(x) histcounts(x,edges,'normalization','cdf'), ISIm(lowBP),'UniformOutput',false);
H_ml=cellfun(@(x) histcounts(x,edges,'normalization','cdf'), ISIml(lowBP),'UniformOutput',false);
isi_t=edges(1:end-1)+BW/2;
semilogx(isi_t',median(cell2mat(H_m)),':','color',M_color,'linewidth',2)
hold on
semilogx(isi_t',median(cell2mat(H_ml)),':','color',ML_color,'linewidth',2)

legend 'bursty M' 'bursty ML' 'non-burst M' 'non-burst ML')
legend box off
title(['hiBP ' num2str(numel(these_units)) ' lowBP ' num2str(numel(lowBP))]);
legend("Position",[0.14268,0.66279,0.28241,0.20216])


% Fig. 3f violin plots MI 

sen=sig_any;
%sen=m;
insen=setdiff(1:numel(thislayer.dblZetaP), sen);

fig=figure;
colors=[L_color; [.4 .4 .4];[0 .2 .7];];
v=violinplot([light_only_mi xs_s ys_s ],{'L','M','ML'},'ShowBox',true,'ShowData',true,'MarkerSize',10,'ShowWhiskers',false,'BoxColor',[0 0 0],'ViolinColor',colors);
for j=1:numel(v)
    v(j).MedianPlot.SizeData=100;
    v(j).MedianPlot.MarkerFaceColor=[0 0 0];
end
ylabel 'Modulation Index'
title(S.stats_x_y.p_sr)
box off
set(gca,'fontsize',14)
ylim([-1 1])
xlim([0 4])
plot([0 4],[0 0])
set(gcf,'Position', [814.3333  352.3333  420.6667  406.0000])


% Fig. 3g Response probability per trial, per unit.  

%first convert counts/trial to 0 or 1 and then average

sen=sig_any;
Cbase=TBL{M}{i}.Cbase(sen);
prob_spont=cellfun(@(x) mean(double(x>0)),Cbase);
Ctest=TBL{M}{i}.Ctest(sen);
probs_mech=cellfun(@(x) mean(double(x>0)),Ctest);
Ctest=TBL{ML}{i}.Ctest(sen);
probs_ml=cellfun(@(x) mean(double(x>0)),Ctest);
Ctest=TBL{L}{i}.Ctest(sen);
probs_l=cellfun(@(x) mean(double(x>0)),Ctest);

fig=figure;
colors=[L_color; [.4 .4 .4];[0 .2 .7];];
%v=violinplot([probs_l probs_mech probs_ml ],{'L','M','ML'},'ShowBox',true,'ShowData',true,'MarkerSize',10,'ShowWhiskers',false,'BoxColor',[0 0 0],'ViolinColor',colors);
for j=1:numel(v)
    v(j).MedianPlot.SizeData=100;
    v(j).MedianPlot.MarkerFaceColor=[0 0 0];
end