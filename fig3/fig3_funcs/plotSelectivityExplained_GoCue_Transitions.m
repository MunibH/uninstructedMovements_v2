function plotSelectivityExplained_GoCue_Transitions(nprez,cd_allrez,cdrez,me,time,sav)

% find average movement bout length (excluding go cue)
[~,gcix] = min(abs(time - 0));
ix = gcix - 1;

figure;
hold on;
for sessix = 1:numel(me)
    tempdat = me(sessix).move;
    for trix = 1:size(tempdat,2)
        % find bouts of non-movement
        [qstart, qlen, qcounts] = ZeroOnesCount(~tempdat(:,trix));
        % find bouts of movement
        [mstart, mlen, mcounts] = ZeroOnesCount(tempdat(:,trix));
        
        medat = me(sessix).data(:,trix);
        t = 1:numel(medat);
        
        qix = [];
        for i = 1:qcounts
            qix = [qix qstart(i):(qstart(i)+(qlen(i)-1))];
        end
        qmask = ismember(t,qix); % logical array, size of trial length, 1 where animal is quiet

        mix = [];
        for i = 1:mcounts
            mix = [mix mstart(i):(mstart(i)+(mlen(i)-1))];
        end
        mmask = ismember(t,mix); % logical array, size of trial length, 1 where animal is moving

        % now can use these masks plus lens to find transitions. HOW HOW
        % HOW HOW HOW HOW HOW HOW???

        clf
        toplot = me(sessix).data(:,trix);
        plot(toplot)
        hold on;
        z = toplot;
        z(~ismember(1:numel(toplot),qix)) = nan;
        plot(z,'r')
        z = toplot;
        z(~ismember(1:numel(toplot),mix)) = nan;
        plot(z,'g')
        yline(me(sessix).moveThresh,'k:')
        pause
    end

%     % find bouts of non-movement
%     tempdat = ~me(sessix).move;
%     [qstart, qlen, qcounts] = ZeroOnesCount(tempdat(:));
%     qlen = qlen / 100; % in seconds
%     % -start is index where non-movement bout starts
%     % -len is num indices of non-movement bout
%     % -counts is the number of non-movement bouts found
% 
%     % find bouts of movement
%     tempdat = me(sessix).move;
%     [mstart, mlen, mcounts] = ZeroOnesCount(tempdat(:));
%     mlen = mlen / 100; % in seconds
%     % -start is index where movement bout starts
%     % -len is num indices of movement bout
%     % -counts is the number of movement bouts found
% 
%     % merge any bouts of non-movement that are less than 0.015 seconds into
%     % the surrounding movement bout
%     mask = qlen < 0.015;
%     qstart(mask) = [];
% 
%     mask = mlen < 0.015;
%     mstart(mask) = [];

end
% merge


tt = find(mask,1,'first')
qstart(tt)


tt = [nan diff(mstart)];
% num indicies < 0.015 s


aa = me(1).data;
aa = aa(:);
figure; 
hold on;
plot((1:length(aa))./100,aa)
for i = 1:numel(mstart)
    xline(mstart(i)/100,'r:'); hold on
    xline(qstart(i)/100,'g:'); hold on
    if i > 200
        break
    end
end
yline(me(1).moveThresh)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% temp = linspecer(12,'qualitative');
% cols = temp(1:4,:);
% cols(4,:) = temp(7,:);

cols(1,:) = [50, 191, 83];
cols(2,:) = [255, 128, 48];
cols(3,:) = [232, 53, 226];
cols(4,:) = [53, 226, 232];
cols = cols./255;


sample = mode(rez(1).ev.sample) - mode(rez(1).align);
delay  = mode(rez(1).ev.delay) - mode(rez(1).align);


sel = allrez.selexp;

lw = 3;
alph = 0.7;
f = figure; ax = axes(f); hold on;

for i = 1:size(sel,2)
    tempdat = squeeze(sel(:,i,:));
    tempmean = nanmean(tempdat,2);
    temperr = nanstd(tempdat,[],2) ./ sqrt(numel(rez));
    shadedErrorBar(rez(1).time,tempmean,temperr,...
                   {'Color',cols(i,:),'LineWidth',lw},alph,ax);
end


xline(sample,'k--','LineWidth',2)
xline(delay,'k--','LineWidth',2)
xline(0,'k--','LineWidth',2)

xlabel('Time (s) from go cue')
ylabel('Selectivity Explained')
% legend('Total selectivity','early','late','go','early + late + go')
xlim([rez(1).time(1),2])
ax = gca;
ax.FontSize = 14;

if sav
%     mysavefig()
end


end