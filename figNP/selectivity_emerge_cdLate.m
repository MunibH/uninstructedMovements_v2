% for each bootstrap iteration, when is selectivity first greater than
% presample selectivity

% save('stElsayed_100boot','meta','params','obj','rez','cd_null','cd_potent','boot','-v7.3')

% edges = [-2.46 -2.2];
edges = [-2.35 -2.21];
for i = 1:numel(edges)
    [~,presampleIX(i)] = min(abs( obj(1).time - edges(i)));
end
presampleIX = presampleIX(1):presampleIX(2);

%% across all bootstrap iterations, build distribution of presample selectivity

sm = 21;
smtype = 'zeropad';
for iboot = 1:boot.iters
    
    cdix = 1;
    cddat = mySmooth(squeeze(cd_null_all.cd_proj(:,:,cdix,iboot)),sm,smtype);
    sel.null = squeeze(cddat(:,1,:) - cddat(:,2,:));


    cddat = mySmooth(squeeze(cd_potent_all.cd_proj(:,:,cdix,iboot)),sm,smtype);
    sel.potent = squeeze(cddat(:,1,:) - cddat(:,2,:));


    sel.presample.null(iboot) = std(sel.null(presampleIX));
    sel.presample.potent(iboot) = std(sel.potent(presampleIX));

end

cols = getColors;

% figure;
% ax = gca;
% hold on;
% nbins = 10;
% histogram(sel.presample.null,nbins,'facecolor',cols.null,'edgecolor','none');
% histogram(sel.presample.potent,nbins,'facecolor',cols.potent,'edgecolor','none')


%% for each bootstrap iteration, find when selectivity emerges, relative to presample selectivity
clear firstNull firstPotent
for iboot = 1:boot.iters
    
    cdix = 1;
    cddat = squeeze(cd_null_all.cd_proj(:,:,cdix,iboot));
    null = squeeze(cddat(:,1,:) - cddat(:,2,:));
    sel.null(:,iboot) = null;


    cddat = squeeze(cd_potent_all.cd_proj(:,:,cdix,iboot));
    potent = squeeze(cddat(:,1,:) - cddat(:,2,:));
    sel.potent(:,iboot) = potent;

    stdmult = 4;
    try
        firstNull(iboot) = find(null(presampleIX(end)+1:end)>(stdmult*(sel.presample.null(iboot))),1,'first' );
    catch
        firstNull(iboot) = nan;
    end
    try
        firstPotent(iboot) = find(potent(presampleIX(end)+1:end)>(stdmult*(sel.presample.potent(iboot))),1,'first' );
    catch
        firstPotent(iboot) = nan;
    end


end

firstNull = firstNull(~isnan(firstNull));
firstPotent = firstPotent(~isnan(firstPotent));

firstNull = obj(1).time(firstNull+presampleIX(end)+1);
firstPotent = obj(1).time(firstPotent+presampleIX(end)+1);

%%

% cis = [5 95];
% 
% n = prctile(sel.null,cis,2);
% n = n(:,1); % lower bound
% p = prctile(sel.potent,cis,2);
% p = p(:,1); % lower bound
% 
% f = figure;
% ax = gca;
% hold on;
% plot(obj(1).time,n,'Color',cols.null,'LineWidth',2)
% plot(obj(1).time,p,'Color',cols.potent,'LineWidth',2)
% xline(mode(obj(1).bp.ev.sample)-2.5,'k--')
% xline(mode(obj(1).bp.ev.delay)-2.5,'k--')
% % figure; 
% % ax = nexttile;
% % plot(sel.null,'Color',cols.null); hold on; plot(sel.potent,'Color',cols.potent)
% % ax = nexttile;
% % plot(sel.potent,'Color',cols.potent); hold on; plot(sel.null,'Color',cols.null)

%%



cols = getColors;

f = figure;
f.Position = [680   773   320   205];
% ax = nexttile;
% hold on;
% histogram(firstNull,20,'facecolor',cols.null,'edgecolor','none','facealpha',0.4)
% histogram(firstPotent,20,'facecolor',cols.potent,'edgecolor','none','facealpha',0.4)



ax = gca;
hold on;
nbins = 20;
n = cdfplot(firstNull);
n.Color = cols.null;
n.LineWidth = 2;
p = cdfplot(firstPotent);
p.Color = cols.potent;
p.LineWidth = 2;
xlabel('Time from go cue')
ylabel(['Cumulative probability of' newline ' selectivity emerging'])
grid off
title('')
ax.FontSize = 10;
xline(mode(obj(1).bp.ev.sample)-2.5,'k--')
% xline(mode(obj(1).bp.ev.delay)-2.5,'k--')
% xlim([-2.3,ax.XLim(2)])

h = kstest2(firstNull,firstPotent)

%%
% clear firstNull firstPotent
% for iboot = 1:boot.iters
%     
%     cdix = 1;
%     cddat = squeeze(cd_null_all.cd_proj(:,:,cdix,iboot));
%     null = squeeze(cddat(:,1,:) - cddat(:,2,:));
%     sel.null(:,iboot) = null;
% 
% 
%     cddat = squeeze(cd_potent_all.cd_proj(:,:,cdix,iboot));
%     potent = squeeze(cddat(:,1,:) - cddat(:,2,:));
%     sel.potent(:,iboot) = potent;
% 
% end
% 
% for itime = 1:size(sel.null,1)
%     h.null(itime) = ttest(sel.presample.null,sel.null(itime,:),'Tail','both');
%     h.potent(itime) = ttest(sel.presample.potent,sel.potent(itime,:),'Tail','both');
% end
% 
% figure; 
% ax = gca;
% hold on;
% plot(obj(1).time,h.null,'Color',cols.null,'LineWidth',2)
% plot(obj(1).time,h.potent,'Color',cols.potent,'LineWidth',2)







