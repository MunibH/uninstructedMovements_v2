function ax = plotCDProj_v2(allrez,rez,sav,spacename,plotmiss,plotno,rez_2afc,obj)

pptx.newVersions = [1 0 0]; % 1 - creates new version of pptx, 0 - add to existing if already exists, 1 entry for each figure created here


clrs = getColors();
lw = 2;
alph = 0.3;

sample = mode(obj(1).bp.ev.sample - rez_2afc(1).align);
delay = mode(obj(1).bp.ev.delay - rez_2afc(1).align);


if strcmpi(spacename,'null')
    clrs.rhit = [3, 3, 173] ./ 255;
    clrs.lhit = [173, 3, 3] ./ 255;
elseif strcmpi(spacename,'potent')
    clrs.rhit = [115, 169, 245] ./ 255;
    clrs.lhit = [240, 110, 138] ./ 255;
end



for i = 1:numel(rez_2afc(1).cd_labels) % for each coding direction
    f = figure; hold on
    f.Position = [698   436   343   230];
    f.Renderer = 'painters';
    ax(i) = gca;
    ax(i) = prettifyPlot(ax(i));
    %     ax(i) = nexttile; hold on;
    tempdat = squeeze(allrez.cd_proj(:,i,:,:));

    tempmean = nanmean(tempdat,3);
    temperror = nanstd(tempdat,[],3) ./ 2;%./sqrt(numel(rez));
    % temperror = getCI(tempdat,0);

    for j = 1:size(tempdat,2)
        temp_ = squeeze(tempdat(:,j,:));
%         temperror(:,j) = prctile(temp_,2.5,2);
%         temperror(:,j) = getCI(temp_);
    end
      shadedErrorBar(obj(1).time,tempmean(:,1+4),temperror(:,1+4),{'Color',clrs.rhit,'LineWidth',lw},alph, ax(i))
      shadedErrorBar(obj(1).time,tempmean(:,2+4),temperror(:,2+4),{'Color',clrs.lhit,'LineWidth',lw},alph, ax(i))
    if plotmiss
        sm = 51;
        smtype = 'zeropad';
        shadedErrorBar(obj(1).time,mySmooth(tempmean(:,3+4),sm,smtype),mySmooth(temperror(:,3+4),sm,smtype),{'Color',clrs.rhit,'LineWidth',lw,'LineStyle',':'},alph, ax(i));
        shadedErrorBar(obj(1).time,mySmooth(tempmean(:,4+4),sm,smtype),mySmooth(temperror(:,4+4),sm,smtype),{'Color',clrs.lhit,'LineWidth',lw,'LineStyle',':'},alph, ax(i))
    end

    %     xlim([rez(1).time(1);rez(1).time(end)])
    xlim([obj(1).time(5);2])

    xlabel('Time from go cue (s)')
    ylabel('Activity (a.u.)')
    ax(i).FontSize = 10;
    % title([rez(1).cd_labels{i} ' | ' spacename],'FontSize',7)

    yline(0,'k-');
    xline(sample,'k--','LineWidth',1)
    xline(delay,'k--','LineWidth',1)
    xline(0,'k--','LineWidth',1)

    hold off

end


end