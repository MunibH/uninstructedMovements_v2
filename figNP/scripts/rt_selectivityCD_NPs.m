
% rt = firstTongueRT(obj);
rt = firstJawRT(obj);

%%

close all
clear dat cols trialsBlock trix

nQuantiles = 4;
for sessix = 1:numel(meta)

    t1 = -0.02; t2 = 0;
    [~,ix1] = min(abs(obj(sessix).time - t1));
    [~,ix2] = min(abs(obj(sessix).time - t2));


    cddim = 1;
    cdlatepotent = cd_potent(sessix).trialdat(:,:,cddim);
    cdlatenull = cd_null(sessix).trialdat(:,:,cddim);
    me_ = me(sessix).data;
    thisrt = rt{sessix};
    thisme = nanmean(me(sessix).data(ix1:ix2,:),1);

    % right and left hit trials
    trix.all = cell2mat(params(sessix).trialid(8:9)');
    trix.right = params(sessix).trialid{8};
    trix.left = params(sessix).trialid{9};

    thisrt = thisrt(trix.all);
%     med = median(thisrt);
%     trix.outlierMask = thisrt > 1.1*med | thisrt<=0.0025;
% %     trix.outlierMask = isoutlier(thisrt)';
%     trix.outlierMask = thisrt > 0.13;
%     trix.outlierMask = thisrt < 0.01;
    trix.outlierMask = logical(zeros(size(thisrt)));

    trix.all = trix.all(~trix.outlierMask);
    thisrt = thisrt(~trix.outlierMask);
    

    cdlatepotent = cdlatepotent(:,trix.all);
    cdlatenull = cdlatenull(:,trix.all);
    me_ = me_(:,trix.all);
    thisme = thisme(trix.all);

    
    quantiles = quantile(thisrt,nQuantiles);
    trialsBlock{sessix} = getTrialsByRTBlock(thisrt,quantiles);

%     quantiles = quantile(thisme,nQuantiles);
%     trialsBlock{sessix} = getTrialsByRTBlock(thisme,quantiles);

    % divide each entry of trialsBlock into right and left trials
    % each entry is same size as trix.all (logical)
    for qix = 1:nQuantiles
        right_ = (ismember(trix.all,trix.right)' & trialsBlock{sessix}{qix}); % array of right trial nums size of trix.all for current quantile
        left_ = (ismember(trix.all,trix.left)' & trialsBlock{sessix}{qix}); % array of left trial nums size of trix.all for current quantile
        trials2use{qix}.right = right_;
        trials2use{qix}.left = left_;
    end
    

    for i = 1:nQuantiles
        % make sure trials are correct
        dat(sessix).rt.right{i} = thisrt(trials2use{i}.right);
        dat(sessix).rt.left{i} = thisrt(trials2use{i}.left);

        dat(sessix).me.right{i} = me_(:,trials2use{i}.right);
        dat(sessix).me.left{i} = me_(:,trials2use{i}.left);

        dat(sessix).cdpotent.right{i} = cdlatepotent(:,trials2use{i}.right);
        dat(sessix).cdpotent.left{i} = cdlatepotent(:,trials2use{i}.left);

        dat(sessix).cdnull.right{i} = cdlatenull(:,trials2use{i}.right);
        dat(sessix).cdnull.left{i} = cdlatenull(:,trials2use{i}.left);

    end


end




%%

clear alldat

alldat.rt = cell(nQuantiles,1);
alldat.mesel = cell(nQuantiles,1);
alldat.cdnullsel = cell(nQuantiles,1);
alldat.cdpotentsel = cell(nQuantiles,1);
alldat.cdnull.right = cell(nQuantiles,1); 
alldat.cdpotent.right = cell(nQuantiles,1);
alldat.cdnull.left = cell(nQuantiles,1);
alldat.cdpotent.left = cell(nQuantiles,1);
alldat.me.right = cell(nQuantiles,1);
alldat.me.left = cell(nQuantiles,1);
for qix = 1:nQuantiles
    alldat.rt{qix} = [];
%     alldat.me{qix} = [];
    alldat.mesel{qix} = [];
    alldat.cdnullsel{qix} = [];
    alldat.cdpotentsel{qix} = [];
    alldat.cdnull.right{qix} = [];
    alldat.cdpotent.right{qix} = [];
    alldat.cdnull.left{qix} = [];
    alldat.cdpotent.left{qix} = [];
    alldat.me.right{qix} = [];
    alldat.me.left{qix} = [];
    for sessix = 1:numel(dat)
        alldat.rt{qix} = [ alldat.rt{qix} ; dat(sessix).rt.right{qix}' ; dat(sessix).rt.left{qix}'];
        alldat.cdnullsel{qix} = cat(2, alldat.cdnullsel{qix} , nanmean(dat(sessix).cdnull.right{qix},2) - nanmean(dat(sessix).cdnull.left{qix},2) );
        alldat.cdpotentsel{qix} = cat(2, alldat.cdpotentsel{qix} , nanmean(dat(sessix).cdpotent.right{qix},2) - nanmean(dat(sessix).cdpotent.left{qix},2) );
        alldat.cdnull.right{qix} = cat(2, alldat.cdnull.right{qix} , nanmean(dat(sessix).cdnull.right{qix},2)  );
        alldat.cdpotent.right{qix} = cat(2, alldat.cdpotent.right{qix} , nanmean(dat(sessix).cdpotent.right{qix},2)  );
        alldat.cdnull.left{qix} = cat(2, alldat.cdnull.left{qix} , nanmean(dat(sessix).cdnull.left{qix},2)  );
        alldat.cdpotent.left{qix} = cat(2, alldat.cdpotent.left{qix} , nanmean(dat(sessix).cdpotent.left{qix},2)  );

        %         alldat.me{qix} = cat(2, alldat.me{qix} , cat(2, dat(sessix).me.right{qix}, dat(sessix).me.left{qix}) );
        alldat.me.right{qix} = cat(2, alldat.me.right{qix} , nanmean(dat(sessix).me.right{qix},2)  );
        alldat.me.left{qix} = cat(2, alldat.me.left{qix} , nanmean(dat(sessix).me.left{qix},2)  );
        alldat.mesel{qix} = cat(2, alldat.mesel{qix} , nanmean(dat(sessix).me.right{qix},2) - nanmean(dat(sessix).me.left{qix},2) );
    end
end

%% selectivity plot

sample = mode(obj(1).bp.ev.sample) - 2.5;
tstart = mode(obj(1).bp.ev.bitStart) - 2.5;
delay = mode(obj(1).bp.ev.delay) - 2.5;

sm = 101;

xlims = [tstart 1.5];

c = gray(nQuantiles*4);
c = linspace(0.15,0.8,nQuantiles);
for i = 1:nQuantiles
    cols(i,:) = [c(i) c(i) c(i)];
%     cols(i,:) = c(i*4,:);
end


% cols(end,:) = [0.83 0.83 0.83];

close all
f = figure;
f.Position = [680    84   574   894];
for i = 1:nQuantiles
    ax1 = subplot(2,1,1); hold on;
    plot(obj(1).time,mySmooth(nanmean(alldat.cdnullsel{i},2),sm,'zeropad'),'LineWidth',2,'Color',cols(i,:));
    title('null');
    xline(sample,'k--')
    xline(0,'k--')
    xline(tstart,'k--')
    xline(delay,'k--')
    xlim(xlims)
    ylabel('Selectivity (null)')

    ax2 = subplot(2,1,2); hold on;
    plot(obj(1).time,mySmooth(nanmean(alldat.cdpotentsel{i},2),sm,'zeropad'),'LineWidth',2,'Color',cols(i,:));
    title('potent');
    xline(sample,'k--')
    xline(0,'k--')
    xline(tstart,'k--')
    xline(delay,'k--')
    xlim(xlims)
    ylabel('Selectivity (potent)')

%     ax = subplot(3,1,3); hold on;
%     plot(obj(1).time,mySmooth(nanmean(alldat.mesel{i},2),sm,'zeropad'),'LineWidth',2,'Color',cols(i,:));
%     xline(sample,'k--')
%     xline(0,'k--')
%     xline(tstart,'k--')
%     xline(delay,'k--')
%     xlim(xlims)
%     ylabel('Motion Energy')


%     histogram(alldat.rt{i},'FaceAlpha',0.7,'EdgeColor','none','FaceColor',cols(i,:))
%     xlim([0 0.4])
end

ylims = [ax1.YLim ax2.YLim];
ylims = [min(ylims) max(ylims)];
ax1.YLim = ylims;
ax2.YLim = ylims;
ax1.FontSize = 12;
ax2.FontSize = 12;
ax1.XTick = [];
ax2.XLabel.String = 'Time (s) from go cue';

% linkaxes([ax1,ax2])

% L/R plot

sample = mode(obj(1).bp.ev.sample) - 2.5;
tstart = mode(obj(1).bp.ev.bitStart) - 2.5;
delay = mode(obj(1).bp.ev.delay) - 2.5;


cleft = flipud(redblue(30));
for i = 1:nQuantiles
    lcols(i,:) = cleft(i*3,:);
end

cright = redblue(30);
for i = 1:nQuantiles
    rcols(i,:) = cright(i*3,:);
end

% close all
figure;
for i = 1:nQuantiles
    ax1 = subplot(2,1,1); hold on;
    plot(obj(1).time,mySmooth(nanmean(alldat.cdnull.right{i},2),sm,'zeropad'),'LineWidth',2,'Color',rcols(i,:));
    plot(obj(1).time,mySmooth(nanmean(alldat.cdnull.left{i},2),sm,'zeropad'),'LineWidth',2,'Color',lcols(i,:));
    xline(sample,'k--')
    xline(0,'k--')
    xline(tstart,'k--')
    xline(delay,'k--')
    xlim(xlims)
    ylabel('CD (Null)')

    ax2 = subplot(2,1,2); hold on;
    plot(obj(1).time,mySmooth(nanmean(alldat.cdpotent.right{i},2),sm,'zeropad'),'LineWidth',2,'Color',rcols(i,:));
    plot(obj(1).time,mySmooth(nanmean(alldat.cdpotent.left{i},2),sm,'zeropad'),'LineWidth',2,'Color',lcols(i,:));
    xline(sample,'k--')
    xline(0,'k--')
    xline(tstart,'k--')
    xline(delay,'k--')
    xlim(xlims)
    ylabel('CD (Potent)')

%     ax3 = subplot(3,1,3); hold on;
%     plot(obj(1).time,mySmooth(nanmean(alldat.me.right{i},2),sm,'zeropad'),'LineWidth',2,'Color',rcols(i,:));
%     plot(obj(1).time,mySmooth(nanmean(alldat.me.left{i},2),sm,'zeropad'),'LineWidth',2,'Color',lcols(i,:));
%     xline(sample,'k--')
%     xline(0,'k--')
%     xline(tstart,'k--')
%     xline(delay,'k--')
%     xlim(xlims)
%     ylabel('Motion Energy')

end

ylims = [ax1.YLim ax2.YLim];
ylims = [min(ylims) max(ylims)];
ax1.YLim = ylims;
ax2.YLim = ylims;
ax1.FontSize = 12;
ax2.FontSize = 12;
ax1.XTick = [];
ax2.XLabel.String = 'Time (s) from go cue';
% 
% linkaxes([ax1,ax2])


% f = figure;
% hold on;
% for i = 1:nQuantiles
%     plot(obj(1).time,mySmooth(nanmean(alldat.me.right{i},2),sm,'zeropad'),'LineWidth',2,'Color',rcols(i,:));
%     plot(obj(1).time,mySmooth(nanmean(alldat.me.left{i},2),sm,'zeropad'),'LineWidth',2,'Color',lcols(i,:));
%     xline(sample,'k--')
%     xline(0,'k--')
%     xline(tstart,'k--')
%     xline(delay,'k--')
%     xlim(xlims)
%     ylabel('Motion Energy')
% end
% 
% f = figure;
% hold on;
% for i = 1:nQuantiles
%     plot(obj(1).time,mySmooth(nanmean(alldat.mesel{i},2),sm,'zeropad'),'LineWidth',2,'Color',cols(i,:));
%     xline(sample,'k--')
%     xline(0,'k--')
%     xline(tstart,'k--')
%     xline(delay,'k--')
%     xlim(xlims)
%     ylabel('Motion Energy')
% end
% 
% f = figure;
% for i = 1:nQuantiles
%     ax = subplot(2,1,1);
%     hold on;
%     plot(obj(1).time,mySmooth(nanmean(alldat.me.right{i},2),sm,'zeropad'),'LineWidth',2,'Color',rcols(i,:));
%     xline(sample,'k--')
%     xline(0,'k--')
%     xline(tstart,'k--')
%     xline(delay,'k--')
%     xlim(xlims)
%     ylabel('motion energy')
% 
%     ax = subplot(2,1,2);
%     hold on;
%     plot(obj(1).time,mySmooth(nanmean(alldat.me.left{i},2),sm,'zeropad'),'LineWidth',2,'Color',lcols(i,:));
%     xline(sample,'k--')
%     xline(0,'k--')
%     xline(tstart,'k--')
%     xline(delay,'k--')
%     xlim(xlims)
%     ylabel('motion energy')
% end





