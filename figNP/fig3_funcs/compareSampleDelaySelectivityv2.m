function compareSampleDelaySelectivityv2(meta,obj,allrez_null,rez_null,allrez_potent,rez_potent)
% compares selectivity in late delay to selectvity before that
% compares bootstrap distributions
%%
labels = {'sample_end','delay_end'};
epochs = {'delay','goCue'};
tm = {[-0.1 0], [-0.1 0]}; % in seconds, relative to respective epochs

for i = 1:numel(epochs)
    e1 = mode(obj(1).bp.ev.(epochs{i})) + tm{i}(1) - 2.5;
    e2 = mode(obj(1).bp.ev.(epochs{i})) + tm{i}(2) - 2.5;
    times.(labels{i}) = obj(1).time>e1 & obj(1).time<e2;
end

cdix = 1;

%%
sel = struct();


% Null
sel.hit.null = squeeze(allrez_null.cd_proj(:,cdix,5,:) - allrez_null.cd_proj(:,cdix,6,:));
sel.miss.null = squeeze(allrez_null.cd_proj(:,cdix,7,:) - allrez_null.cd_proj(:,cdix,8,:));

% Potent
sel.hit.potent = squeeze(allrez_potent.cd_proj(:,cdix,5,:) - allrez_potent.cd_proj(:,cdix,6,:));
sel.miss.potent = squeeze(allrez_potent.cd_proj(:,cdix,7,:) - allrez_potent.cd_proj(:,cdix,8,:));


%%

delay.null     = mean(sel.hit.null(times.delay_end,:),1);
delay.nullci   = abs(getCI(delay.null,0));
delay.potent   = mean(sel.hit.potent(times.delay_end,:),1);
delay.potentci = abs(getCI(delay.potent,0));

% loop through time until start of the delay testing period (last 100 ms of
% delay epoch)
% compare selectivity in those time points to the delay testing period
p.null = zeros(size(sel.hit.null,1),1);
p.potent = zeros(size(sel.hit.null,1),1);
f = figure;
hold on;
for i = 1:find(times.delay_end,1,'first')
    hold on;
    cla
    % get 95% CI
    % nci = abs(getCI(sel.hit.null(i,:),0));
    % npci = abs(getCI(sel.hit.potent(i,:),0));
    nn = sel.hit.null(i,:);
    pp = sel.hit.potent(i,:);

    % [~,p.null(i)] = my_ttest(delay.null - nn);
    % [~,p.potent(i)] = my_ttest(delay.potent - pp);
    [~,p.null(i)] = my_ttest(delay.null,nn);
    [~,p.potent(i)] = my_ttest(delay.potent,pp);
    % p.null(i) = sum(nn >= delay.nullci) / numel(nn);
    % p.potent(i) = sum(pp >= delay.potentci) / numel(pp);
    
    cols = getColors;
    histogram(delay.null-nn,30,'facecolor',cols.null);
    histogram(delay.potent-pp,30,'facecolor',cols.potent);
    title([num2str(obj(1).time(i)) ' ' num2str(p.null(i))])
    hold off
    pause
end


%%
close all
cols = getColors;
lw = 1.5;
f = figure;
f.Renderer = 'painters';
ax = gca;
ax = prettifyPlot(ax);
hold on;
plot(obj(1).time,p.null,'color',cols.null,'linewidth',2)
plot(obj(1).time,p.potent,'color',cols.potent,'linewidth',lw)
% plot(obj(1).time,p.null>0.05,'color',cols.null,'linewidth',2)
% plot(obj(1).time,p.potent>0.05,'color',cols.potent,'linewidth',lw)


alph = 0.3;
xlims = [-2.3 0];
tstart = mode(obj(1).bp.ev.bitStart) - mode(obj(1).bp.ev.goCue);
sample = mode(obj(1).bp.ev.sample) - mode(obj(1).bp.ev.goCue);
delay = mode(obj(1).bp.ev.delay) - mode(obj(1).bp.ev.goCue);
gc = 0;

xlim(xlims)
xline(tstart,'k--'); xline(sample,'k--'); xline(delay,'k--'); xline(gc,'k--')
yline(0,'k-')
xlabel('Time from go cue (s)')
ylabel('p-value')
% ylabel('p-value > 0.01')

ylims = ax.YLim;

pat = patch([tm{2}(1) tm{2}(2) tm{2}(2) tm{2}(1)], [ylims(1) ylims(1) ylims(2) ylims(2)], 'k');
pat.EdgeColor = 'none';
pat.FaceColor = 'k';
pat.FaceAlpha = alph;



end