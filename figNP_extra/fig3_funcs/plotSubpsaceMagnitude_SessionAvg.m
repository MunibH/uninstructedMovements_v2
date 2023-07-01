function plotSubpsaceMagnitude_SessionAvg(obj,meta,params,rez,me,cond2use)

close all

cols = getColors;

lw = 2;
alph = 0.2;

for sessix = 1:numel(meta)
    trix = cell2mat(params(sessix).trialid(cond2use)');


    % Null space
    tempdat = rez(sessix).N_null(:,trix,:);
    tempdat = sum(tempdat.^2,3); % sum across dims
    datmean = mean(tempdat,2);   % mean across trials
    nulldat(:,sessix) = datmean;

    % Potent space
    tempdat = rez(sessix).N_potent(:,trix,:);
    tempdat = sum(tempdat.^2,3); % sum across dims
    datmean = mean(tempdat,2);   % mean across trials
    potentdat(:,sessix) = datmean;

    % Motion energy
    tempdat = me(sessix).data(:,trix);
    datmean = mean(tempdat,2);   % mean across trials
    medat(:,sessix) = datmean;


end

f = figure;

yyaxis left
ax = gca;
hold on;

% dat = nulldat;
% mu = mean(dat,2);
% stderr = std(dat,[],2) ./ sqrt(size(dat,2));
% shadedErrorBar(obj(1).time,mu,stderr,{'Color',cols.null,'LineWidth',lw},alph,ax);

dat = potentdat;
mu = mean(dat,2);
stderr = std(dat,[],2) ./ sqrt(size(dat,2));
shadedErrorBar(obj(1).time,mu,stderr,{'Color',cols.potent,'LineWidth',lw,'LineStyle','-'},alph,ax);
ylabel('Activity (a.u.)')

yyaxis right

dat = medat;
mu = mean(dat,2);
stderr = std(dat,[],2) ./ sqrt(size(dat,2));
shadedErrorBar(obj(1).time,mu,stderr,{'Color','k','LineWidth',lw,'LineStyle','-'},alph,ax);


sample = mode(obj(1).bp.ev.sample) - 2.5;
delay = mode(obj(1).bp.ev.delay) - 2.5;
sample = mode(obj(1).bp.ev.sample) - 2.5;
delay = mode(obj(1).bp.ev.delay) - 2.5;
xline(0,'k--')
xline(sample,'k--')
xline(delay,'k--')
xlim([-2.5 2.5])
xlabel('Time (s) from go cue')
ylabel('Motion energy (a.u.)')
ax.FontSize = 10;
axpotent = ax;

end





















