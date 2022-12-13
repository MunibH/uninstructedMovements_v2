
% rt = firstTongueRT(obj);
rt = firstJawRT(obj);

%%

close all
clear dat cols

nBlocks = 4;
for sessix = 1:numel(meta)

%     t1 = -0.5; t2 = -0.01;
%     [~,ix1] = min(abs(obj(sessix).time - t1));
%     [~,ix2] = min(abs(obj(sessix).time - t2));

    potent = sum(rez(sessix).N_potent.^2,3);
    cdlatepotent = cd_potent(sessix).trialdat(:,:,2);
    null = sum(rez(sessix).N_null.^2,3);
    cdlatenull = cd_null(sessix).trialdat(:,:,2);
    me_ = me(sessix).data;
    thisrt = rt{sessix};

    figure; subplot(2,1,1); histogram(thisrt,40); subplot(2,1,2); plot(thisrt); pause; close all; continue;

    % right and left hit trials
    trix.all = cell2mat(params(sessix).trialid(2:3)');
    [~,sortix] = sort(thisrt(trix.all),'descend');
    trix.all = trix.all(sortix);
    nTrials = numel(trix.all);
    trialsPerBlock = floor(nTrials/nBlocks);
    trialsBlock = mat2cell(trix.all,diff([0:trialsPerBlock:nTrials-1,nTrials]));
    if numel(trialsBlock{end}) < trialsPerBlock
        trialsBlock(end) = [];
    end
    for i = 1:numel(trialsBlock)
        dat.me{i} = me_(:,trialsBlock{i});
        dat.cdpotent{i} = cdlatepotent(:,trialsBlock{i});
        dat.cdnull{i} = cdlatenull(:,trialsBlock{i});
    end

    xlims = [-2.4 2.4];

    c = gray(nBlocks*3);
    cols(1,:) = c(1,:);
    cols(2,:) = c(3,:);
    cols(3,:) = c(5,:);
    cols(4,:) = c(7,:);

    xl.gc = 0;
    xl.sample =  mode(obj(sessix).bp.ev.sample) - mode(obj(sessix).bp.ev.(params(sessix).alignEvent));
    xl.delay =  mode(obj(sessix).bp.ev.delay) - mode(obj(sessix).bp.ev.(params(sessix).alignEvent));
    xl.trialStart =  mode(obj(sessix).bp.ev.bitStart) - mode(obj(sessix).bp.ev.(params(sessix).alignEvent));
    xlfns = fieldnames(xl);

    f=figure;
    f.Position = [680    42   695   954];
    alph = 0.05;
    fns = fieldnames(dat);
    for i = 1:numel(fns) % data type
        ax = subplot(numel(fns),1,i);
        hold on;
        for j = 1:nBlocks
            if strcmpi(fns{i},'cdnull')
                sm = 21;
            else
                sm = 11;
            end
            temp = mySmooth(dat.(fns{i}){j},sm,'zeropad');
            plot(obj(sessix).time,mean(temp,2),'LineWidth',2.5,'Color',cols(j,:));
        end
        xlim(xlims)
        ylabel([fns{i} ' (a.u.)']);
        for ixl = 1:numel(xlfns)
            xline(xl.(xlfns{ixl}),'k--');
        end

        ax.FontSize = 12;
        if i~=numel(fns)
            ax.XTick = [];
        else
            xlabel('Time (s) from go cue')
        end
        sgtitle([meta(sessix).anm ' ' meta(sessix).date])
    end
    
%     pth = 'C:\Users\munib\Documents\Economo-Lab\code\uninstructedMovements_v2\fig3\figs\me_np_blocks_by_RT';
%     fn = [meta(sessix).anm '_' meta(sessix).date];
%     mysavefig(f,pth,fn)




end






