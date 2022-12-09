

rt = firstTongueRT(obj);

%% simple plots
close all

cols = linspecer(numel(obj));
% trials at top of figures are fast rt, bottom are slow rt
for sessix = 1:numel(rez)
    temprez = rez(sessix);
    tempme = me(sessix);
    trix = cell2mat(params(sessix).trialid(2:3)');

    thisrt = rt{sessix}(trix);
    [rt_sorted,ix] = sort(thisrt,'ascend');
    trix = trix(ix);


    f = figure;
    f.Position = [118   211   529   704];
    plotme = tempme.data(:,trix);
    imagesc(obj(sessix).time,1:numel(trix),plotme');
    colormap(linspecer)
    lims = clim;
    clim([lims(1) lims(2) / 1])

    f = figure;
    f.Position = [650   211   529   704];
    potent = temprez.N_potent;
    potent = potent(:,trix,:);
    potent = sum(potent.^2,3);
    imagesc(obj(sessix).time,1:numel(trix),potent');
    colormap(linspecer);
    lims = clim;
    clim([lims(1) lims(2) / 1.5])

    f = figure;
    f.Position = [1182   211   529   704];
    null = temprez.N_null;
    null = null(:,trix,:);
    null = sum(null.^2,3);
    imagesc(obj(sessix).time,1:numel(trix), null');
    colormap(linspecer)
    lims = clim;
    clim([lims(1) lims(2) / 1.5])

    pause
    close all

end

%% simple plots
close all

cols = linspecer(numel(obj));

allme = cell(3,1); % allme{rtgroup}(time,trials)
allpotent = cell(3,1);
allnull = cell(3,1);
for sessix = 1:numel(rez)

    temprez = rez(sessix);
    tempme = me(sessix);
    trix = cell2mat(params(sessix).trialid(2:3)');

    thisrt = rt{sessix}(trix);

    fast = prctile(thisrt,33);
    slow = prctile(thisrt,66);

    trials_by_rt{1} = trix(thisrt<fast); % fast
    trials_by_rt{2} = trix(thisrt>=fast & thisrt<=slow); % med
    trials_by_rt{3} = trix(thisrt>slow); % slow


    for i = 1:3 % for each rt group
        allme{i} = tempme.data(:,trials_by_rt{i});

        potent = temprez.N_potent;
        potent = potent(:,trials_by_rt{i},:);
        allpotent{i} = sum(potent.^2,3);

        null = temprez.N_null;
        null = null(:,trials_by_rt{i},:);
        allnull{i} = sum(null.^2,3);
    end

    alph = 0.05;
    xlims = [-1.5 0];

    f = figure;
    f.Position = [118   511   529   311];
    ax = gca;
    hold on;
    h = 360; s = 0;
    cols(1,:) = [h s 0] ./ 360;
    cols(2,:) = [h s 100] ./ 360;
    cols(3,:) = [h s 200] ./ 360;
    cols = hsv2rgb(cols);
    dat = allme;
    shadedErrorBar(obj(sessix).time,mean(dat{1},2),std(dat{1},[],2) ./ sqrt(size(dat{1},2)),{'Color',cols(1,:),'LineWidth',2},alph,ax)
    shadedErrorBar(obj(sessix).time,mean(dat{2},2),std(dat{2},[],2) ./ sqrt(size(dat{2},2)),{'Color',cols(2,:),'LineWidth',2},alph,ax)
    shadedErrorBar(obj(sessix).time,mean(dat{3},2),std(dat{3},[],2) ./ sqrt(size(dat{3},2)),{'Color',cols(3,:),'LineWidth',2},alph,ax)
    xlim(xlims)

    f = figure;
    f.Position = [650   511   529   311];
    ax = gca;
    hold on;
    h = 330; s = 250;
    cols(1,:) = [h s 150] ./ 360;
    cols(2,:) = [h s 250] ./ 360;
    cols(3,:) = [h s 360] ./ 360;
    cols = hsv2rgb(cols);
    dat = allpotent;
    shadedErrorBar(obj(sessix).time,mean(dat{1},2),std(dat{1},[],2) ./ sqrt(size(dat{1},2)),{'Color',cols(1,:),'LineWidth',2},alph,ax)
    shadedErrorBar(obj(sessix).time,mean(dat{2},2),std(dat{2},[],2) ./ sqrt(size(dat{2},2)),{'Color',cols(2,:),'LineWidth',2},alph,ax)
    shadedErrorBar(obj(sessix).time,mean(dat{3},2),std(dat{3},[],2) ./ sqrt(size(dat{3},2)),{'Color',cols(3,:),'LineWidth',2},alph,ax)
    xlim(xlims)

    f = figure;
    f.Position = [1182   511   529  311];
    ax = gca;
    hold on;
    h = 130; s = 250;
    cols(1,:) = [h s 150] ./ 360;
    cols(2,:) = [h s 250] ./ 360;
    cols(3,:) = [h s 360] ./ 360;
    cols = hsv2rgb(cols);
    dat = allnull;
    shadedErrorBar(obj(sessix).time,mean(dat{1},2),std(dat{1},[],2) ./ sqrt(size(dat{1},2)),{'Color',cols(1,:),'LineWidth',2},alph,ax)
    shadedErrorBar(obj(sessix).time,mean(dat{2},2),std(dat{2},[],2) ./ sqrt(size(dat{2},2)),{'Color',cols(2,:),'LineWidth',2},alph,ax)
    shadedErrorBar(obj(sessix).time,mean(dat{3},2),std(dat{3},[],2) ./ sqrt(size(dat{3},2)),{'Color',cols(3,:),'LineWidth',2},alph,ax)
    xlim(xlims)

    
    h = 100; s = 0;
    cols(1,:) = [h s 0] ./ 360;
    cols(2,:) = [h s 100] ./ 360;
    cols(3,:) = [h s 200] ./ 360;
    cols = hsv2rgb(cols);
    f = figure; hold on;
    for j = 1:3
        t1 = -0.5; t2 = -0.02;
        [~,ix1] = min(abs(obj(sessix).time - t1));
        [~,ix2] = min(abs(obj(sessix).time - t2));
        dat1 = allnull{j}(ix1:ix2);
        dat2 = allpotent{j}(ix1:ix2);
        mdl = fitlm(dat1(:),dat2(:));
        h = plot(mdl);
        h(1).Marker = 'o';
        h(1).MarkerSize = 1;
        h(1).MarkerFaceColor = cols(j,:);
        h(1).MarkerEdgeColor = 'w';
        h(2).Color = cols(j,:);
        h(2).LineWidth = 2;

        for i = 3:4
            h(i).Visible = 'off';
%             h(i).LineStyle = '-';
%             h(i).Color = cols(j,:);
        end
        %     r2(sessix) = mdl.Rsquared.Ordinary;
        %     title(['R^2 = ' num2str(r2(sessix))])
        xlabel('Late delay CD_Late_Null')
        ylabel('Late delay CD_Late_Potent')

        ax = gca;
        ax.Legend.Visible = 'off';
        ax.FontSize = 15;
        %     ylim([-12 ax.YLim(2)])
    end
    hold off
    


    pause
    %     break
    close all

end














