function plotEventTimes(ax,evtimes)

ylims = ax.YLim;

fnames = fieldnames(evtimes);
for i = 1:numel(fnames)
    this = evtimes.(fnames{i});
    plot(ax,[this,this],ylims,'k--')
    % xline(ax,evtimes.(fnames{i}),'k--');
end
ylim(ax,ylims);

end