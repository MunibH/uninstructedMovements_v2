function ax = plotCDSelectivity(meta,obj,allrez_null,rez_null,allrez_potent,rez_potent)

labels = {'latesample','latedelay'};
epochs = {'delay','goCue'};
tm = {[-0.82 -0.42], [-0.42 -0.02]}; % in seconds, relative to respective epochs

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


%% plot

condfns = fieldnames(sel);
npfns = fieldnames(sel.hit);

alph = 0.2;
lw = 2;

cols = getColors;

xlims = [-2.4 2.5];

sample = mode(obj(1).bp.ev.sample) - mode(obj(1).bp.ev.goCue);
delay = mode(obj(1).bp.ev.delay) - mode(obj(1).bp.ev.goCue);
gc = 0;

sm = 21;

for i = 1:numel(condfns)
    f = figure;
    f.Position = [698   436   343   230];
    title(condfns{i})
    ax = gca;
    hold on;

    if strcmpi(condfns{i},'hit')
        ls = '-';
    else
        ls = '-.';
    end
    for j = 1:numel(npfns)
        if strcmpi(npfns{j},'null')
            c = cols.null;
        else
            c = cols.potent;
        end
        temp = mySmooth(sel.(condfns{i}).(npfns{j}),sm,'zeropad');
        mu = nanmean(temp,2);
%         sig = nanstd(temp,[],2) ./ sqrt(numel(obj));
        sig = getCI(temp);
        shadedErrorBar(obj(1).time, mu, sig,  {'Color',c,'LineWidth',lw,'LineStyle',ls},alph,ax)
%         plot(obj(1).time,temp,'Color',c,'LineWidth',1)
    end

    xlabel('Time from go cue (s)')
    ylabel('Selectivity in CDTrialType (a.u.)')
    xlim(xlims)
    xline(sample,'k--')
    xline(delay,'k--')
    xline(0,'k--')
end



end










