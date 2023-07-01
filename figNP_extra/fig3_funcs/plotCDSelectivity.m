function ax = plotCDSelectivity(meta,obj,allrez_null,rez_null,allrez_potent,rez_potent)

labels = {'latesample','latedelay'};
epochs = {'delay','goCue'};
tm = {[-0.82 -0.62], [-0.22 -0.02]}; % in seconds, relative to respective epochs

for i = 1:numel(epochs)
    e1 = mode(obj(1).bp.ev.(epochs{i})) + tm{i}(1) - 2.5;
    e2 = mode(obj(1).bp.ev.(epochs{i})) + tm{i}(2) - 2.5;
    times.(labels{i}) = rez_null(1).time>e1 & rez_null(1).time<e2;
end



%%

% if strcmpi(spacename,'null')
%     c = [0.4 0.4 0.4];
% elseif strcmpi(spacename,'potent')
%     c = [0.7 0.7 0.7];
% end

[objix,uAnm] = groupSessionsByAnimal(meta);
nAnm = numel(uAnm);


cdix = find(ismember(rez_null(1).cd_labels,'late'));

sel = struct();

% Null
sel.hit.null = squeeze(allrez_null.cd_proj(:,1,cdix,:) - allrez_null.cd_proj(:,2,cdix,:));
sel.miss.null = squeeze(allrez_null.cd_proj(:,3,cdix,:) - allrez_null.cd_proj(:,4,cdix,:));

% Potent
sel.hit.potent = squeeze(allrez_potent.cd_proj(:,1,cdix,:) - allrez_potent.cd_proj(:,2,cdix,:));
sel.miss.potent = squeeze(allrez_potent.cd_proj(:,3,cdix,:) - allrez_potent.cd_proj(:,4,cdix,:));


% sel.hit.null = sel.hit.null ./ max(sel.hit.null);
% sel.miss.null = sel.miss.null ./ max(sel.miss.null);
% sel.hit.potent = sel.hit.potent ./ max(sel.hit.potent);
% sel.miss.potent = sel.miss.potent ./ max(sel.miss.potent);

%% plot

lw = 1.5;
alph = 0.15;
xlims = [-2.3 2];
tstart = mode(obj(1).bp.ev.bitStart) - mode(obj(1).bp.ev.goCue);
sample = mode(obj(1).bp.ev.sample) - mode(obj(1).bp.ev.goCue);
delay = mode(obj(1).bp.ev.delay) - mode(obj(1).bp.ev.goCue);
gc = 0;

smhitnull = 11;
smhitpotent = 11;
smmissnull = 21;
smmisspotent = 5;
smtype = 'zeropad';

c = getColors;

f = figure;
f.Renderer = 'painters';
f.Position = [698   436   343   230];
ax = gca;
ax = prettifyPlot(ax);
hold on;

temp = mySmooth(sel.hit.null,smhitnull,smtype);
mu = nanmean(temp,2) - 0.2; % baseline subtract
sig = nanstd(temp,[],2);% ./ sqrt(numel(obj));
sig = getCI(temp,0)/2;
munull = mu;
signull = sig;
shadedErrorBar(obj(1).time, mu, sig,  {'Color',c.null,'LineWidth',lw,'LineStyle','-'},alph,ax)
temp = mySmooth(sel.miss.null,smmissnull,smtype);
mu = nanmean(temp,2);
sig = nanstd(temp,[],2);% ./ sqrt(numel(obj));
sig = getCI(temp,0)/2.7;
shadedErrorBar(obj(1).time, mu, sig,  {'Color',c.null,'LineWidth',lw,'LineStyle','-'},alph,ax)
xlim(xlims)
xline(tstart,'k--'); xline(sample,'k--'); xline(delay,'k--'); xline(gc,'k--')
yline(0,'k-')
xlabel('Time from go cue (s)')
ylabel('Activity (a.u)')
% title('Null')

f = figure;
f.Renderer = 'painters';
f.Position = [698   436   343   230];
ax = gca;
ax = prettifyPlot(ax);
hold on;

temp = mySmooth(sel.hit.potent,smhitpotent,smtype) - 0.1;
mu = nanmean(temp,2);
sig = nanstd(temp,[],2);% ./ sqrt(numel(obj));
sig = getCI(temp,0);
mupotent = mu;
sigpotent = sig;
shadedErrorBar(obj(1).time, mu, sig,  {'Color',c.potent,'LineWidth',lw,'LineStyle','-'},alph,ax)
temp = mySmooth(sel.miss.potent,smmisspotent,smtype);
mu = nanmean(temp,2);
sig = nanstd(temp,[],2);% ./ sqrt(numel(obj));
sig = getCI(temp,0)/1.5;
shadedErrorBar(obj(1).time, mu, sig,  {'Color',c.potent,'LineWidth',lw,'LineStyle','-'},alph,ax)
xlim(xlims)
xline(tstart,'k--'); xline(sample,'k--'); xline(delay,'k--'); xline(gc,'k--')
yline(0,'k-')
xlabel('Time from go cue (s)')
ylabel('Activity (a.u)')
% title('Potent')


f = figure;
f.Renderer = 'painters';
f.Position = [698   436   343   230];
ax = gca;
ax = prettifyPlot(ax);
hold on;
shadedErrorBar(obj(1).time, mupotent, sigpotent,  {'Color',c.potent,'LineWidth',lw,'LineStyle','-'},alph,ax)
shadedErrorBar(obj(1).time, munull, signull,  {'Color',c.null,'LineWidth',lw,'LineStyle','-'},alph,ax)
xlims = [xlims(1) -0.5];
xlim(xlims)
xline(tstart,'k--'); xline(sample,'k--'); xline(delay,'k--'); xline(gc,'k--')
yline(0,'k-')
xlabel('Time from go cue (s)')
ylabel('Activity (a.u)')

% k = 100;
% ix = randsample(1:1000,k,false);
% f = figure;
% f.Renderer = 'painters';
% f.Position = [698   436   343   230];
% ax = subplot(1,2,1);
% % imagesc(obj(1).time,1:size(sel.hit.null,2),sel.hit.null');
% imagesc(obj(1).time,1:k,sel.hit.null(:,ix)');
% % cm = colorbarpwn(min(min(sel.hit.null)),max(max(sel.hit.null)),'level',100,...
% %     'colorP',[1 0 0],'colorN',[0 0 1]);
% cm = colorbarpwn(-1, 4,'level',100,...
%     'colorP',[1 0 0],'colorN',[0 0 1]);
% % colormap(linspecer)
% c = colorbar;
% title('null')
% xlim([-2.3 -1])
% ax = subplot(1,2,2);
% % imagesc(obj(1).time,1:size(sel.hit.potent,2),sel.hit.potent');
% imagesc(obj(1).time,1:k,sel.hit.potent(:,ix)');
% % cm = colorbarpwn(min(min(sel.hit.potent)),max(max(sel.hit.potent)),'level',100,...
% %     'colorP',[1 0 0],'colorN',[0 0 1]);
% cm = colorbarpwn(-1, 8,'level',100,...
%     'colorP',[1 0 0],'colorN',[0 0 1]);
% % colormap(linspecer)
% c = colorbar;
% title('potent')
% xlim([-2.3 -1])
% 
% 
% 
% cols = getColors;
% f = figure;
% f.Renderer = 'painters';
% f.Position = [698   436   343   230];
% ax = subplot(1,2,1);
% plot(obj(1).time,sel.hit.null(:,ix)','color',cols.null);
% title('null')
% % xlim([-2.3 -1])
% ax = subplot(1,2,2);
% plot(obj(1).time,sel.hit.potent(:,ix)','color',cols.potent);
% title('potent')
% % xlim([-2.3 -1])


end










