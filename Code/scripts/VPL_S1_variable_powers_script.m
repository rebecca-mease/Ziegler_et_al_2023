%% figure 1D and Fig. s8a
load S1_data_tables_L_variable_power
PosCount=zeros(1,numel(conds));
for J=1:numel(conds)
    Tbl=TBL{J};
    thislayer=Tbl{end};
    light_sen=unique([find(thislayer.dblZetaP<=.05 & thislayer.dR>0);find(thislayer.P_rs<=.05 & thislayer.dR>0)]);
    cellcount=size(thislayer,1);
    PosCount(J)=numel(light_sen);
end

Tbl=TBL{p5};  %screen for light sens at highest power
thislayer=Tbl{3}; %just Layer 6 neurons by depth
light_sen=unique([find(thislayer.dblZetaP<=.05 & thislayer.dR>0);find(thislayer.P_rs<=.05 & thislayer.dR>0)]);
depths=thislayer.depths(light_sen);

xs=[1 2 5 10 15] %actual laser power
Rates=[]
for i=3  %just layer 6 by depth
    for J=1:numel(conds)
        thislayer=TBL{J}{i};
        Rates(:,J)=thislayer.m_rate(light_sen);
    end
end
mrate=max(Rates');

%screen for appreciable responsive above 1 Hz
normal_L6=(find(mrate>1));
Rates=Rates(normal_L6,:);

fig=figure
yyaxis right
hold on
errorbar(xs,mean(Rates),std(Rates)/sqrt(numel(normal_L6)),'o-k','LineWidth',2);
set(gca,'XTick',xs,'xticklabel',conds)
xlabel 'power'
ylabel 'mean L6 unit rate [Hz]'
xlim([0 max(xs)+1])

yyaxis left
ylabel('proportion L6 units with sig. positive light response')
bar(xs+.5,PosCount/cellcount)

clearvars -except basedir
cd(basedir)
load VPL_data_tables_L_variable_power

% Figure S8a (VPL data)

PosCount=zeros(1,numel(conds));
%
for J=1:numel(conds)
    Tbl=TBL{J};
    thislayer=Tbl{end};
    light_sen=unique([find(thislayer.dblZetaP<=.05 & thislayer.dR>0);find(thislayer.P_rs<=.05 & thislayer.dR>0)]);
    cellcount=size(thislayer,1);
    PosCount(J)=numel(light_sen);
end

Tbl=TBL{p5};  %screen for light sens
thislayer=Tbl{end};%just VPL
light_sen=unique([find(thislayer.dblZetaP<=.05 & thislayer.dR>0);find(thislayer.P_rs<=.05 & thislayer.dR>0)]);

xs=[1 2 5 10 15]
Rates=[]
for i=1  %vpl
    for J=1:numel(conds)
        %cl=colors(i,:);
        thislayer=TBL{J}{i};
        Rates(:,J)=thislayer.m_rate(light_sen);
    end
end
mrate=max(Rates');

fig=figure
yyaxis right
hold on
errorbar(xs,mean(Rates),std(Rates)/sqrt(numel(light_sen)),'o-k','LineWidth',2)
set(gca,'XTick',xs,'xticklabel',conds)
xlabel 'power'
ylabel 'mean L6 unit rate [Hz]'
xlim([0 max(xs)+1])

yyaxis left
ylabel('proportion VPL units with sig. positive light response')
bar(xs+.5,PosCount/cellcount)

% Fig. S8b Change in burst probability
sen=light_sen;
C={};Cb={}; BP={}; BPb={};
for J=1:numel(conds)
    C{J}=TBL{J}{i}.EventcTest;
    counts=cellfun(@(x) numel(x),C{J},'UniformOutput',true);
    bc=cellfun(@(x) numel(find(x>1)),C{J},'UniformOutput',true);
    singles=counts-bc;
    Cb{J}=TBL{J}{i}.Eventcbase;
    counts=cellfun(@(x) numel(x),Cb{J},'UniformOutput',true);
    bcb=cellfun(@(x) numel(find(x>1)),Cb{J},'UniformOutput',true);
    singles=counts-bcb;
    BP{J}=cellfun(@(x) numel(find(x>1))/numel(x),C{J},'UniformOutput',true);
    BPb{J}=cellfun(@(x) numel(find(x>1))/numel(x),Cb{J},'UniformOutput',true);
end


D=((cell2mat(BP))-(cell2mat(BPb)));
s=sum(D');
Dc=D(~isnan(s),:);
kruskalwallis(Dc);
title('Change in burst proportion vs. laser power')
xlabel power
ylabel '\Delta BP'
%v=violinplot((Dc),{'1','2','3','4','5'},'ShowBox',false,'ShowData',false);
%%