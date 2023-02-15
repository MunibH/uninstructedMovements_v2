
% find move onset probabilities for each session
% compare to selectivity onset in potent space for each session

%% boostrap to find move onset

close all

cond2use = [8 9];

cols = getColors;
col{1} = cols.rhit;
col{2} = cols.lhit;

k = 1000;
sm = 21;
xlims = [-2.4 2.4];

f = figure;
t = tiledlayout('flow');
for s = 1:numel(obj)

    ax = nexttile;
    hold on;
    title([meta(s).anm ' ' meta(s).date])
    for ic = 1:numel(cond2use)
        clear mu sig onsetTimes dat

        edges = [-0.3 0];
        tt = mode(obj(s).bp.ev.sample) - mode(obj(s).bp.ev.goCue);
        for i = 1:numel(edges)
            [~,ix(i)] = min(abs(time - edges(i) - tt));
        end

        % get data
        trix = cell2mat(params(s).trialid(cond2use(ic))');
        dat = mySmooth(me(s).data(:,trix),sm,'zeropad');
        % normalize (position =0 at gc and max after gc is 1)
        dat = dat - dat(ps,:);
        dat = dat ./ max(dat);

        % remove inf
        [r,c] = find(isinf(dat));
        dat(:,c) = [];

        % cut off first 15 samples b/c of smoothing artifacts
        dat(1:15,:) = [];
        time = obj(1).time(16:end);


        nsamp = floor(size(dat,2) / 4);
        for i = 1:k
            % sample trials
            temp = datasample(dat,nsamp,2,"Replace",true);
            % detrend mean trace-  linearly detrended the mean trace based on its value between Tdelay + 0.6 s and Tdelay + 1.2 s
            % don't know what that means really, so just calculating mean and std
            mu = nanmean(temp,2);
            sig = nanstd(temp,[],2);
            % We identified the time point in which movement exceeds three times the standard deviation of the baseline before the Go cue (100 ms window).
            % going to find time points in which movements exceeds mean + 15%*std
            mubase = mean(mu(ix(1):ix(2)),1);
            sigbase = mean(sig(ix(1):ix(2)),1);

            onsetTimes{i} = time(mu > (mubase + 0.15*sigbase));

            %             afterSampleFlinch = mode(obj(s).bp.ev.sample) + 0.4 - mode(obj(s).bp.ev.goCue);
            %             onsetTimes{i} = time(find(onsetTimes{i} > afterSampleFlinch, 20));
            %             onsetTimes{i} = onsetTimes{i}(onsetTimes{i} > afterSampleFlinch);

        end

        onsetTimes = cell2mat(onsetTimes);


        histogram(onsetTimes,'FaceColor',col{ic},'EdgeColor','none','Normalization','pdf','FaceAlpha',0.4)

    end
    xlim(xlims);

end


xlabel(t,'Time from go cue (s)')
ylabel(t,'Probability of movement onset')




%% boostrap to find selectivity onset in potent space

close all

cond2use_ta = [1 2]; % right and left hits, corresponding to trial-avg projs onto n/p
cond2use_st = [8 9]; % right and left hits, corresponding to single-trial projs onto n/p
modparams.subTrials = 35;

cols = getColors;
col{1} = cols.rhit;
col{2} = cols.lhit;

k = 1000;
sm = 21;
xlims = [-2.4 2.4];

f = figure;
t = tiledlayout('flow');
for s = 1:numel(obj)

    ax = nexttile;
    hold on;
    title([meta(s).anm ' ' meta(s).date])
    clear mu sig onsetTimes dat sel potent selDims pref nonpref

    edges = [mode(obj(s).bp.ev.delay) mode(obj(s).bp.ev.goCue)] - mode(obj(s).bp.ev.goCue);

    % get data

    [potent, selDims, pref, nonpref] = getSortedSelectivityPref(obj(s),params(s),rez(s),edges,cond2use_ta,cond2use_st,modparams,'potent');
    dat = rez(s).N_potent(:,:,selDims);
    pref = pref(selDims);
    nonpref = nonpref(selDims);

    for i = 1:k
        
        for ic = 1:numel(cond2use_st)
            trix = cell2mat(params(s).trialid(cond2use_st(ic))');
            nsamp = floor(numel(trix) / 4);
            trix = datasample(trix,nsamp,"Replace",true);
            tempdat{ic} = dat(:,trix,:);
        end

        % selectivity
        tempdatmu = cellfun(@(x) squeeze(nanmean(x,2)), tempdat, 'UniformOutput',false);
        tempdatmu = cat(3,tempdatmu{1},tempdatmu{2});
        for ii = 1:size(tempdatmu,2) % for each dim
            sel(:,ii) = tempdatmu(:,ii,pref(ii)) - tempdatmu(:,ii,nonpref(ii));     % (time,dims), trial-averaged activity
        end
        sel = mySmooth(sel,sm,'zeropad');

        edges = [-0.3 0];
        tt = mode(obj(s).bp.ev.sample) - mode(obj(s).bp.ev.goCue);
        for iii = 1:numel(edges)
            [~,ix(iii)] = min(abs(obj(s).time - edges(iii) - tt));
        end

        mu = nanmean(sel,2);
        sig = nanstd(sel,[],2);

        mubase = mean(mu(ix(1):ix(2)),1);
        sigbase = mean(sig(ix(1):ix(2)),1);

        onsetTimes{i} = obj(s).time(mu > (mubase + 0.2*sigbase));

    end

    onsetTimes = cell2mat(onsetTimes);

    histogram(onsetTimes,'FaceColor',cols.potent,'EdgeColor','none','Normalization','pdf','FaceAlpha',1)

    xlim(xlims)
end



xlabel(t,'Time from go cue (s)')
ylabel(t,'Probability of selectivity onset')



%% boostrap to find selectivity onset in null space

close all

cond2use_ta = [1 2]; % right and left hits, corresponding to trial-avg projs onto n/p
cond2use_st = [8 9]; % right and left hits, corresponding to single-trial projs onto n/p
modparams.subTrials = 35;

cols = getColors;
col{1} = cols.rhit;
col{2} = cols.lhit;

k = 1000;
sm = 21;
xlims = [-2.4 2.4];

f = figure;
t = tiledlayout('flow');
for s = 1:numel(obj)

    ax = nexttile;
    hold on;
    title([meta(s).anm ' ' meta(s).date])
    clear mu sig onsetTimes dat sel null selDims pref nonpref

    edges = [mode(obj(s).bp.ev.delay) mode(obj(s).bp.ev.goCue)] - mode(obj(s).bp.ev.goCue);

    % get data

    [null, selDims, pref, nonpref] = getSortedSelectivityPref(obj(s),params(s),rez(s),edges,cond2use_ta,cond2use_st,modparams,'null');
    if isempty(selDims)
        continue
    end
    dat = rez(s).N_null(:,:,selDims);
    pref = pref(selDims);
    nonpref = nonpref(selDims);

    for i = 1:k
        
        for ic = 1:numel(cond2use_st)
            trix = cell2mat(params(s).trialid(cond2use_st(ic))');
            nsamp = floor(numel(trix) / 4);
            trix = datasample(trix,nsamp,"Replace",true);
            tempdat{ic} = dat(:,trix,:);
        end

        % selectivity
        tempdatmu = cellfun(@(x) squeeze(nanmean(x,2)), tempdat, 'UniformOutput',false);
        tempdatmu = cat(3,tempdatmu{1},tempdatmu{2});
        for ii = 1:size(tempdatmu,2) % for each dim
            sel(:,ii) = tempdatmu(:,ii,pref(ii)) - tempdatmu(:,ii,nonpref(ii));     % (time,dims), trial-averaged activity
        end
        sel = mySmooth(sel,sm,'zeropad');

        edges = [-0.3 0];
        tt = mode(obj(s).bp.ev.sample) - mode(obj(s).bp.ev.goCue);
        for iii = 1:numel(edges)
            [~,ix(iii)] = min(abs(obj(s).time - edges(iii) - tt));
        end

        mu = nanmean(sel,2);
        sig = nanstd(sel,[],2);

        mubase = mean(mu(ix(1):ix(2)),1);
        sigbase = mean(sig(ix(1):ix(2)),1);

        onsetTimes{i} = obj(s).time(mu > (mubase + 0.2*sigbase));

    end

    onsetTimes = cell2mat(onsetTimes);

    histogram(onsetTimes,'FaceColor',cols.null,'EdgeColor','none','Normalization','pdf','FaceAlpha',1)

    xlim(xlims)
end



xlabel(t,'Time from go cue (s)')
ylabel(t,'Probability of selectivity onset')











