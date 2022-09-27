function plotAvgJawVelocityDuringStim_singleTrials(meta,obj,dfparams,params,kin,kinfeats,feats2plot,cond2plot,sav)

dfparams.plt.ms = {'.','.','.','.','.','.'};


[~,mask] = patternMatchCellArray(kin(1).featLeg,feats2plot,'any');
featix = find(mask);

times = [0 0.8]; % relative to alignEv (should be delay)
for i = 1:numel(times)
    [~,ix(i)] = min(abs(dfparams.time - times(i)));
end

% get avg feature value across time and trials for each condition
toplot = cell(numel(feats2plot),1);
for k = 1:numel(featix)
    for j = 1:numel(cond2plot)
        toplot{k}{j} = [];
        for i = 1:numel(obj)
            trials = params(i).trialid{cond2plot(j)};
            temp = kinfeats{i}(ix(1):ix(2),trials,featix(k));
            temp = nanmean(temp,1);
            toplot{k}{j} = [toplot{k}{j} ; temp'];
        end
    end
end



%% plot
f = figure; hold on;
% f.Position = [680   205   477   773];
t = tiledlayout('flow');
for k = 1:numel(featix)
    for j = 1:numel(cond2plot)
        ax = nexttile; hold on;
        s = histogram(toplot{k}{j},50,'FaceColor',dfparams.plt.color{cond2plot(j)},'EdgeColor','none');
        xlabel('Motion Energy')
        ylabel('Frequency')
        ax.FontSize = 12;
        xlim([0 50])
    end




    %     mins = min([ax.XLim(1) ax.YLim(1)]);
    %     maxs = max([ax.XLim(2) ax.YLim(2)]);
    %     ax.XLim = [mins maxs];
    %     ax.YLim = [mins maxs];

    %     plot(ax.XLim,ax.YLim,'k--','LineWidth',2)


end


if sav
    pth = [ 'C:\Users\munib\Documents\Economo-Lab\code\uninstructedMovements\mc_stim\figs\'];
    fn = [ 'AvgJawVelStim_NoStim' ];
    mysavefig(f,pth,fn)
end








end