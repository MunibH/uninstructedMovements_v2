%% single trial cross-corr b/w (me,potent) and (me,null)
clear r
clear mu

dt = params(1).dt;

% maxlag = int64( ./ dt);
maxlag = 2 / dt;

tt = [-2 0]; % only use time points before go cue
for i = 1:numel(tt)
    [~,tix(i)] = min(abs(obj(1).time - tt(i)));
end
tix = tix(1):tix(2);

cond2use = [2 3];
for sessix = 1:numel(meta)

    trix = cell2mat(params(sessix).trialid(cond2use)');

    %     thisme = zscore(me(sessix).data(tix,trix));
    thisme = normalize(me(sessix).data(tix,trix));
    %     null = sum(rez(sessix).N_null(:,trix,:).^2,3);
    %     potent = sum(rez(sessix).N_potent(:,trix,:).^2,3);

    for t = 1:numel(trix)

        % null
        ndims = rez(sessix).dPrep;
        for d = 1:ndims
            null = normalize(rez(sessix).N_null(tix,trix(t),d).^2);
            [r.null{sessix}(:,t,d),lagtm] = xcorr(thisme(:,t),null,maxlag,'normalized'); % {session}(time,trial,dim)
        end

        % potent
        ndims = rez(sessix).dMove;
        for d = 1:ndims
            potent = normalize(rez(sessix).N_potent(tix,trix(t),d).^2);
            [r.potent{sessix}(:,t,d),lagtm] = xcorr(thisme(:,t),potent,maxlag,'normalized');
        end
    end


end

mu.null = cellfun(@(x) squeeze(nanmean(nanmean(x,3),2)),r.null,'UniformOutput',false); % mean across trials for each session
mu.potent = cellfun(@(x) squeeze(nanmean(nanmean(x,3),2)),r.potent,'UniformOutput',false);
nullcc = cell2mat(mu.null);
potentcc = cell2mat(mu.potent);

%% plots
col = getColors();
lw = 2;
alph = 0.2;
f = figure;
f.Position = [680   694   383   284];
f.Renderer = "painters";
ax = gca;
ax = prettifyPlot(ax);
hold on;
% shadedErrorBar(lagtm*dt,mean(nullcc,2),std(nullcc,[],2)./sqrt(numel(meta)),{'Color',col.null,'LineWidth',lw},alph,ax);
% shadedErrorBar(lagtm*dt,mean(potentcc,2),std(potentcc,[],2)./sqrt(numel(meta)),{'Color',col.potent,'LineWidth',lw},alph,ax);
shadedErrorBar(lagtm*dt,mean(nullcc,2),getCI(nullcc),{'Color',col.null,'LineWidth',lw},alph,ax);
shadedErrorBar(lagtm*dt,mean(potentcc,2),getCI(potentcc),{'Color',col.potent,'LineWidth',lw},alph,ax);
% xlim([-0.5 0.5])

xlabel('Time lag, movement and subspace activity (s)')
ylabel('Correlation')
% title('st-elsayed')