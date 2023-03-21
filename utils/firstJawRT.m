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
        if all(isnan(tm))
            temp = (0:numel(tm))./400;
            tm = temp(1:end-1);
        end
        [~,gcix] = min(abs(tm - gc));

        if gcix == numel(tm) % i think these are no response trials
            rt{sessix}(trix) = nan;
            continue;
        end

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

        % find time from go cue to first tongue visible
        % if tongue is visible during go cue presentation - make sure that
        % the tongue is visible for more than a few frames
        x = squeeze(tongue_gc(:,1,:));
        vis.mask = logical(sum(~isnan(x),2)); % just need either visibile in x or y
        if vis.mask(1) % if tongue visible during go cue prez
            [istart , iend] = ZeroOnesCount(vis.mask);
            if (iend(1)-istart(1)) < 5
                tongue_gc(istart(1):iend(1),:,:) = nan;
            end
        end


        rtix = find(~isnan(tongue_gc),1,'first'); 
        if rtix > numel(tm) % this means that no tongue found
            rt{sessix}(trix) = nan;
            continue;
        else % check that the tongue isn't mostly invisible and just outlier blips are giving an rtix
            [istart , iend] = ZeroOnesCount(vis.mask);
            if ~any((iend-istart)>5)
                rt{sessix}(trix) = nan;
                continue;
            end
        end

        % jaw after go cue
        jaw_gc = jaw(tm > gc, :, :);

        
        % go cue to rtix
        ix2use = gcix_gc:rtix;

        % in time from go cue to first tongue visible, use jaw position to
        % identify a better estimate of reaction time
        dat2use = jaw_gc(ix2use,2);

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
%         plot(tm,squeeze(tongue(:,2,1)),'k','LineWidth',2);
%         plot(tm,jaw(:,2),'r','LineWidth',2)
%         xline(tm(rtix),'b--')
%         xline(tm(gcix),'k--')
% %         xlim([tm(gcix)-1 tm(gcix)+1])
%         xlim([tm(1) tm(gcix)+1])
%         pause
%         clf
%         hold off

    end

end

if numel(rt) == 1
    rt = rt{1}; % just return an array rather than a cell if one session
end

%%

end