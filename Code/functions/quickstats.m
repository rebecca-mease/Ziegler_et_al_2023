function q=quickstats(x)
x=x(~isnan(x));
q.mean=nanmean(x);
q.std=nanstd(x);
q.stder=nanstd(x)/sqrt(sum(~isnan(x))-1);
q.med=nanmedian(x);
q.quart=quantile(x,[.25 .75]);
