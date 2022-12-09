function rt = firstTongueRT(obj)
% get reaction time from tongue trajectories as first time tongue visible
% after go cue

rt = cell(numel(obj),1);

% figure
for sessix = 1:numel(obj)
    traj = obj(sessix).traj;
    for trix = 1:obj(sessix).bp.Ntrials
        gc = obj(sessix).bp.ev.goCue(trix);

        view = 1;
        featNames = traj{view}(trix).featNames;
        [~,ix] = patternMatchCellArray(featNames,{'tongue'},'all');
        ts1 = traj{view}(trix).ts(:,[1 2], ix);

        view = 2;
        featNames = traj{view}(trix).featNames;
        [~,ix] = patternMatchCellArray(featNames,{'tongue'},'all');
        ts2 = traj{view}(trix).ts(:,[1 2], ix);

        tm = traj{view}(trix).frameTimes - 0.5;


        if size(ts1,1) ~= size(ts2,1)
            t = min(size(ts1,1),size(ts2,1));
            ts1 = ts1(1:t,:,:);
            ts2 = ts2(1:t,:,:);
            tm = tm(1:t);
        end
        ts = cat(3,ts1,ts2);

        tsaftergc = ts(tm > gc,:,:);

        tmaftergc = tm(tm > gc);

        rtix = find(~isnan(tsaftergc),1,'first');
        if rtix > size(tsaftergc,1)
            rtix = [];
        end
        if ~isempty(rtix)
            rt{sessix}(trix) = tmaftergc(rtix) - gc;
        else
            rt{sessix}(trix) = nan;
        end

%         tm = tm(tm>gc);
%         trials = cell2mat(params(sessix).trialid(2:5)');
%         if rt{sessix}(trix) > 0.15 && ismember(trix,trials)
%             plot(tm,squeeze(ts(:,1,:)),'r');
%             hold on;
%             plot(tm,squeeze(ts(:,2,:)),'b')
%             xline(rt{sessix}(trix)+gc,'k--')
%             xline(gc,'k--')
%             pause
%             hold off
%         end
    end

end

%%

end