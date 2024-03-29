function ax = plotCDProj(allrez,rez,sav,spacename,plotmiss,plotno)

pptx.newVersions = [1 0 0]; % 1 - creates new version of pptx, 0 - add to existing if already exists, 1 entry for each figure created here


clrs = getColors();
lw = 2;
alph = 0.3;

sample = mode(rez(1).ev.sample - rez(1).align);
delay = mode(rez(1).ev.delay - rez(1).align);


if strcmpi(spacename,'null')
    clrs.rhit = [3, 3, 173] ./ 255;
    clrs.lhit = [173, 3, 3] ./ 255;
elseif strcmpi(spacename,'potent')
    clrs.rhit = [115, 169, 245] ./ 255;
    clrs.lhit = [240, 110, 138] ./ 255;
end



for i = 1:numel(rez(1).cd_labels) % for each coding direction
    f = figure; hold on
    f.Position = [698   436   343   230];
    ax(i) = gca;
    %     ax(i) = nexttile; hold on;
    tempdat = squeeze(allrez.cd_proj(:,:,i,:));

    tempmean = nanmean(tempdat,3);
    temperror = nanstd(tempdat,[],3); %./sqrt(numel(rez));

    for j = 1:size(tempdat,2)
        temp_ = squeeze(tempdat(:,j,:));
%         temperror(:,j) = prctile(temp_,2.5,2);
%         temperror(:,j) = getCI(temp_);
    end
%     if ~plotmiss
      shadedErrorBar(rez(1).time,tempmean(:,1),temperror(:,1),{'Color',clrs.rhit,'LineWidth',lw},alph, ax(i))
        shadedErrorBar(rez(1).time,tempmean(:,2),temperror(:,2),{'Color',clrs.lhit,'LineWidth',lw},alph, ax(i))
%     end
    if plotmiss
        sm = 51;
        smtype = 'zeropad';
        shadedErrorBar(rez(1).time,mySmooth(tempmean(:,3),sm,smtype),mySmooth(temperror(:,3),sm,smtype),{'Color',clrs.rhit,'LineWidth',lw,'LineStyle',':'},alph, ax(i));
        shadedErrorBar(rez(1).time,mySmooth(tempmean(:,4),sm,smtype),mySmooth(temperror(:,4),sm,smtype),{'Color',clrs.lhit,'LineWidth',lw,'LineStyle',':'},alph, ax(i))
    end

    %     xlim([rez(1).time(1);rez(1).time(end)])
    xlim([rez(1).time(5);2])

    xlabel('Time from go cue (s)')
    ylabel('Activity (a.u.)')
    ax(i).FontSize = 10;
    title([rez(1).cd_labels{i} ' | ' spacename],'FontSize',7)


    xline(sample,'k--','LineWidth',1)
    xline(delay,'k--','LineWidth',1)
    xline(0,'k--','LineWidth',1)

    if ~strcmpi(rez(1).cd_labels{i},'ramping')
        curmodename = rez(1).cd_labels{i};
        shadetimes = rez(1).time(rez(1).cd_times.(curmodename));
        x = [shadetimes(1)  shadetimes(end) shadetimes(end) shadetimes(1)];
        y = [ax(i).YLim(1) ax(i).YLim(1) ax(i).YLim(2) ax(i).YLim(2)];
        %     y = [-60 -60 50 50];
        fl = fill(x,y,'r','FaceColor',[93, 121, 148]./255);
        fl.FaceAlpha = 0.3;
        fl.EdgeColor = 'none';


        ylim([y(1) y(3)]);
        %     ylim([-60 50]);
    end

    if sav
        %         pptx.figPath = 'C:\Users\munib\Documents\Economo-Lab\code\uninstructedMovements_v2\fig3\figs\st_elsayed\cd'; % path to save pptx file to
        %         pptx.filename = ['st_elsayed_' spacename '_cd']; % name of ppt file, if already exists, will add fig to new slide
        %         pptx.newVersion = pptx.newVersions(i);
        %         pptx.slideTitle = [rez(1).cd_labels{i} ' | ' spacename];
        %         pptx.fig = f;
        %         myExportToPPTX(pptx)

        pth = 'C:\Users\munib\Documents\Economo-Lab\code\uninstructedMovements\fig1_v2\figs\cd1';
        fn = ['cd_' lower(fns{i}(3:end-7))];
        mysavefig(f,pth,fn);
        %                 exportfig(f, fullfile(pth,fn),'Format','eps','Color','rgb')
    end

    hold off

end


end