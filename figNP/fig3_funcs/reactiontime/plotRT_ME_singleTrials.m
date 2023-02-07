function plotRT_ME_singleTrials(obj,params,me,rt)

% how correlated is ME at end of delay last 20ms with RT
edges = [-0.05 0];
for i = 1:numel(edges)
    [~,ix(i)] = min(abs(obj(1).time - edges(i)));
end

cond2use = [8 9]; % right,left hits, including early
cols = getColors;
col{1} = cols.rhit;
col{2} = cols.lhit;

f = figure;
t = tiledlayout('flow');
for sessix = 1:numel(obj)
    ax = nexttile;
    hold on;
    for c = 1:numel(cond2use)
        trix = cell2mat(params(sessix).trialid(cond2use(c))');

        thisrt = rt{sessix}(trix);
        thisme = mean(me(sessix).data(ix(1):ix(2),trix),1);
        
        scatter(thisrt,thisme,10,'MarkerEdgeColor','none','MarkerFaceColor',col{c})
    end
    hold off

end
xlabel(t,'reaction time (s)')
ylabel(t,'motion energy')