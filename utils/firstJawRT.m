function rt = firstJawRT(obj)
% get reaction time from jaw (and tongue) trajectories as first time jaw
% significantly moved after go cue

rt = cell(numel(obj),1);

% figure
for sessix = 1:numel(obj)
    traj = obj(sessix).traj;
    for trix = 1:obj(sessix).bp.Ntrials
        gc = obj(sessix).bp.ev.goCue(trix);

        % get tongue trajs

        view = 1;
        featNames = traj{view}(trix).featNames;
        [~,ix] = patternMatchCellArray(featNames,{'tongue'},'all');
        ts1 = traj{view}(trix).ts(:,[1 2], ix);

        view = 2;
        featNames = traj{view}(trix).featNames;
        [~,ix] = patternMatchCellArray(featNames,{'tongue'},'all');
        ts2 = traj{view}(trix).ts(:,[1 2], ix);

        tm = traj{view}(trix).frameTimes - 0.5;
        [~,gcix] = min(abs(tm - gc));

        if size(ts1,1) ~= size(ts2,1)
            t = min(size(ts1,1),size(ts2,1));
            ts1 = ts1(1:t,:,:);
            ts2 = ts2(1:t,:,:);
            tm = tm(1:t);
        end
        tongue = cat(3,ts1,ts2);

        % get jaw trajs
        view = 1;
        featNames = traj{view}(trix).featNames;
        [~,ix] = patternMatchCellArray(featNames,{'jaw'},'all');
        jaw = traj{view}(trix).ts(:,[1 2], ix);

        
        % trim data to time from go cue to first tongue visible

        % first tongue after go cue
        tongue_gc = tongue(tm >= gc,:,:);
        tm_gc = tm(tm >= gc);
        [~,gcix_gc] = min(abs(tm_gc - gc));
        rtix = find(~isnan(tongue_gc),1,'first'); 
        if rtix > numel(tm) % this means that no tongue found
            rt{sessix}(trix) = nan;
            continue;
        end

        % jaw after go cue
        jaw_gc = jaw(tm > gc, :, :);

        
        % go cue to rtix
        ix2use = gcix_gc:rtix;


try
        dat2use = jaw_gc(ix2use,2);
catch
    'a'
end
        mu = mean(dat2use);
        sigma = std(dat2use);
        thresh = mu+sigma;

        [~,crossix] = min(abs(dat2use - thresh));
        rtix = crossix + gcix;

        if rtix > numel(tm)
            rtix = [];
        end
        if ~isempty(rtix)
            rt{sessix}(trix) = tm(rtix) - gc;
        else
            rt{sessix}(trix) = nan;
        end



%         hold on
%         plot(tm,squeeze(tongue(:,2,:)),'k');
%         plot(tm,jaw(:,2),'r')
%         xline(tm(rtix),'g--')
%         xline(tm(gcix),'k--')
%         xlim([tm(gcix)-0.5 tm(gcix)+0.5])
%         pause
%         clf
%         hold off

    end

end

%%

end