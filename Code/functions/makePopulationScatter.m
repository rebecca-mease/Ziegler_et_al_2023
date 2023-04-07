function [S figs h dR] = makePopulationScatter(xs_s,ys_s,xs_ns,ys_ns,cl)


figs(1)=figure
if ~isempty(xs_ns)
    plot(xs_ns,ys_ns,'.','markersize',10,'color',[.7 .7 .7]);
end
hold on
plot(xs_s,ys_s,'o','markersize',5,'color',cl,'linewidth',1);
plot([0:.1:max(xs_s)],[0:.1:max(xs_s)],'linewidth',1.5);
figs(2)=figure;
pairs=[xs_s ys_s];
base=ones(size(pairs));dim=size(base,2); for d=1:dim;base(:,d)=base(:,d)*d;end


indices=unique(find(~isnan(xs_s) & ~isnan(ys_s)));

xs_s=xs_s(indices);
ys_s=ys_s(indices);

stats_xs=quickstats(xs_s);
stats_ys=quickstats(ys_s);

[p] = signrank(xs_s, ys_s);
[r rp]=corrcoef(xs_s,ys_s);
r=r(2);
stats_x_y.p_sr=p;
stats_x_y.r=r;
stats_x_y.rp=rp;

dR=pairs(indices,2)-pairs(indices,1);
subplot(1,2,1)
plot(base',pairs','o-');
xlim([.5 2.5])
xlabel([num2str(sum(dR<=0)) ' ' num2str(sum(dR>0))]);

subplot(1,2,2)
plot(xs_s,dR,'o');
[r rp]=corrcoef(xs_s,dR);
r=r(2);

stats_dR=quickstats(dR);
stats_dR.r=r;
stats_dR.rp=rp;

S.stats_xs=stats_xs;
S.stats_ys=stats_ys;
S.stats_x_y=stats_x_y;
S.stats_dR=stats_dR;

figs(3)=figure

h=histogram(dR);
