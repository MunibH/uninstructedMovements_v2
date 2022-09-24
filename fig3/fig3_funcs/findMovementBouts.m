function findMovementBouts(nprez,cd_allrez,cdrez,me,obj)


time = obj(1).time;
dt = mode(diff(time));

boutlen = 0.5; % 500 ms


% find average movement bout length (excluding go cue)
[~,gcix] = min(abs(time - 0));
ix = gcix - 1;

figure;
hold on;
for sessix = 1:numel(me)
    tempdat = me(sessix).move;
    tempp = me(sessix).data;
    tempp = tempp(:);
    figure; plot(tempp(1:10000))
    for trix = 1:size(tempdat,2)

        align = mode(obj(sessix).bp.ev.goCue);
        [~,sampix] = min(abs(time - (obj(sessix).bp.ev.sample(trix) - align)));

        medat = me(sessix).data(:,trix);
        medat(isnan(medat)) = 0;
        medat = lowpass(medat,5,100);
        
%         mevel = gradient( gradient(medat));
%         figure; plot(mevel)
% figure;        plot(medat)


figure; plot(medat);
hold on;plot(150*(medat > me(sessix).moveThresh))

        medat = smoothdata(medat,'movmean',10,'omitnan'); % make this a filter that finds the minimum in every dt bin

        t = 1:numel(medat);

        thresh = me(sessix).moveThresh + (nanstd(medat(1:sampix)) * 0.75);
%         thresh = me(sessix).moveThresh;
        move = medat > thresh;

        % find bouts of non-movement
        [qstart, qend, qlen, qcounts] = ZeroOnesCount(~tempdat(:,trix));

        % find bouts of movement
        [mstart, mend, mlen, mcounts] = ZeroOnesCount(move);
        
        % merge bouts of movement
        temp = mstart;
        for i = 1:numel(mstart)-1
            % if there are not quiet bouts b/w current mstart and next mstart
            % remove next mstart 
            % or if quiet bout very small
            s = mstart(i);
            ss = mstart(i+1);
            xs = mstart(i):mstart(i+1);
            if ~any(ismember(qstart,xs))
                temp(i) = nan;
            end
        end
        mstart = temp;

        % merge bouts of nonmovement
        temp = qstart;
        for i = 1:numel(qstart)-1
            % if there are not quiet bouts b/w current mstart and next mstart
            % remove next mstart 
            % or if quiet bout very small
            s = qstart(i);
            ss = qstart(i+1);
            xs = qstart(i):qstart(i+1);
            if ~any(ismember(mstart,xs))
                temp(i) = nan;
            end
        end
        qstart = temp;

        qend(isnan(qstart)) = [];
        mend(isnan(mstart)) = [];


%         mstart((mlen./100)<0.2) = [];
%         qstart((qlen./100)<0.02) = [];


        figure; clf;
        hold on;
        plot(medat)
        yline(thresh)
        for i = 1:numel(mstart)
            xline(mstart,'g--','LineWidth',2)
            xline(mend,'b--','LineWidth',2)
        end
        for i = 1:numel(qstart)
            xline(qstart,'r--','LineWidth',2)
            xline(qend,'m--','LineWidth',2)
        end

        pause


%         % find bouts of non-movement
%         [qstart, qend, qlen, qcounts] = ZeroOnesCount(~tempdat(:,trix));
%         % find bouts of movement
%         [mstart, mend, mlen, mcounts] = ZeroOnesCount(tempdat(:,trix));
% 
% 
%         %% find movement bouts greater than 50ms where subsequent quiet bout
%         
% 
%         figure;
%         hold on;
%         plot(medat)
%         for i = 1:numel(mstart)
%             xline(mstart(i),'K:','LineWidth',2)
%         end
%         thresh = nanstd(medat) * 0.10;
%         yline(me(sessix).moveThresh,'k--','LineWidth',1)
%         yline(me(sessix).moveThresh + thresh,'r--','LineWidth',1)
%         xline(sampix,'g--')
        %%
    end
end





end