

rt = firstTongueRT(obj);


%%

close all

cols = linspecer(numel(obj));

nBlocks = 10;
for sessix = 1:numel(rez)

    t1 = -0.5; t2 = -0.01;
    [~,ix1] = min(abs(obj(sessix).time - t1));
    [~,ix2] = min(abs(obj(sessix).time - t2));

    %     trix.all = find(obj(sessix).bp.L & obj(sessix).bp.hit);
    trix.all = 1:obj(sessix).bp.Ntrials;
    trix.all = trix.all';
    %     trix.all = sort(cell2mat(params(sessix).trialid(2:3)'));

    potent = rez(sessix).N_potent;
    cdpotent = cd_potent(sessix).trialdat;
    null = rez(sessix).N_null;
    cdnull = cd_null(sessix).trialdat;
    me_ = me(sessix).data;

    % for each block of trials, divide that block into high and low me
    % trials
    nTrials = numel(trix.all);
    trialsPerBlock = floor(nTrials/nBlocks);
    trialsBlock = mat2cell(trix.all,diff([0:trialsPerBlock:nTrials-1,nTrials]));
    if numel(trialsBlock{end}) < trialsPerBlock
        trialsBlock(end) = [];
    end
    trix.high = [];
    trix.low = [];
    for i = 1:numel(trialsBlock)
        tempme = mean(me_(ix1:ix2,trialsBlock{i}),1); % avg late delay me for current block of trials
        [~,isorted] = sort(tempme);
        thresh = prctile(tempme,50);
        trix.high = [trix.high ; trialsBlock{i}(tempme >= thresh)];
        trix.low =  [trix.low ; trialsBlock{i}(tempme < thresh)];
    end
    if numel(trix.high) ~= numel(trix.low)
        nTrials = min(numel(trix.high),numel(trix.low));
        trix.high = randsample(trix.high,nTrials,false);
        trix.low = randsample(trix.low,nTrials,false);
    end

    %     cddim = 2;
    %     dat.potent{1} = mean(cdpotent(ix1:ix2,trix.high,cddim),1); % ssm all dims
    %     dat.potent{2} = mean(cdpotent(ix1:ix2,trix.low,cddim),1);
    %     dat.null{1} = mean(mySmooth(cdnull(ix1:ix2,trix.high,cddim),21),1);
    %     dat.null{2} = mean(mySmooth(cdnull(ix1:ix2,trix.low,cddim),21),1);
    %
    %     cols = linspecer(4,'qualitative');
    %     for i = 1:2
    %         scatter(dat.null{i},dat.potent{i},30,'filled','MarkerFaceColor',cols(i+2,:));
    %     end



    dat.potent{1} = mean(potent(:,trix.high,:).^2,3); % ssm all dims
    dat.potent{2} = mean(potent(:,trix.low,:).^2,3);
    dat.null{1} = mean(null(:,trix.high,:).^2,3);
    dat.null{2} = mean(null(:,trix.low,:).^2,3);



    cols = linspecer(4,'qualitative');
    figure; hold on;
    for i = 1:2
        scatter(mean(dat.null{i}(ix1:ix2,:),1),mean(dat.potent{i}(ix1:ix2,:),1),20,'filled','MarkerFaceColor',cols(i+2,:),'MarkerEdgeColor','w');
%         scatter(mean(dat.null{i}(73:104,:),1),mean(dat.potent{i}(ix1:ix2,:),1),20,'filled','MarkerFaceColor',cols(i+2,:),'MarkerEdgeColor','w');
    end

    %     temp = zeros(5,nTrials);
    %     temp(:,trix.high) = 1;
    %     figure; imagesc(temp); colormap(gray); ax = gca; ax.YTick = []; ax.FontSize = 12; xlabel('Trials'); colorbar;
    %
    %     figure;
    %     imagesc(me_(:,[trix.high;trix.low])')
    %
    %     figure; hold on;
    %     plot(mean(me_(:,trix.high),2),'Color',[0 0 0])
    %     plot(mean(me_(:,trix.low),2),'Color',[0.5 0.5 0.5])
    %
    %     figure; hold on;
    %     plot(mean(sum(potent(:,trix.high,:).^2,3),2),'Color',[0 0 0]);
    %     plot(mean(sum(potent(:,trix.low,:).^2,3),2),'Color',[0.5 0.5 0.5]);
    %
    %     figure; hold on;
    %     plot(mySmooth(mean(sum(null(:,trix.high,:).^2,3),2),15,'reflect'),'Color',[0 0 0]);
    %     plot(mySmooth(mean(sum(null(:,trix.low,:).^2,3),2),15,'reflect'),'Color',[0.5 0.5 0.5]);
    %
    pause
    close all

end



