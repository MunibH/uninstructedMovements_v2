

rt = firstTongueRT(obj);


%% simple plots
close all

% feat = 'jaw_ydisp_view2';
feat = 'jaw_yvel_view2';
featix = find(ismember(kin(1).featLeg, feat));


cols = linspecer(numel(obj));
% trials at top of figures are fast rt, bottom are slow rt
for sessix = 1:numel(rez)
    tempme = me(sessix);
    trix{1} = params(sessix).trialid{2};
    trix{2} = params(sessix).trialid{3};

    thisrt{1} = rt{sessix}(trix{1});
    thisrt{2} = rt{sessix}(trix{2});

    t1 = -0.02; t2 = 0;
    [~,ix1] = min(abs(obj(sessix).time - t1));
    [~,ix2] = min(abs(obj(sessix).time - t2));

    ts{1} = nanmean(kin(sessix).dat(ix1:ix2,trix{1},featix),1);
    ts{2} = nanmean(kin(sessix).dat(ix1:ix2,trix{2},featix),1);

    plotme{1} = mean(tempme.data(ix1:ix2,trix{1}),1);
    plotme{2} = mean(tempme.data(ix1:ix2,trix{2}),1);

    f = figure; hold on;    
%     scatter(thisrt{1},plotme{1},10,'b');
%     scatter(thisrt{2},plotme{2},10,'r');
    scatter(thisrt{1},ts{1},10,'b');
    scatter(thisrt{2},ts{2},10,'r');

    title([ num2str( corr(thisrt{1}',ts{1}') ) ' | ' num2str( corr(thisrt{2}',ts{2}') )  ])

    pause
    close all

end













