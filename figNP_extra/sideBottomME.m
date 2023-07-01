%%

clear,clc,close all

% add paths for data loading scripts, all fig funcs, and utils
utilspth = 'C:\Users\munib\Documents\Economo-Lab\code\uninstructedMovements_v3';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));
rmpath(genpath(fullfile(utilspth,'fig1/')));
rmpath(genpath(fullfile(utilspth,'fig2/')));
rmpath(genpath(fullfile(utilspth,'mc_stim/')));

% add paths for figure specific functions
addpath(genpath(pwd))

clc

%% get file names

datapth = 'C:\Users\munib\Documents\Economo-Lab\data\Side-Bottom ME';

token = {'side','bottom'};

for i = 1:numel(token)
    contents = dir(fullfile(datapth,['*',token{i},'*.mat']));
    fns{:,i} = {contents(:).name}';
end


%%

badsess = [2 3 7 8 9 11];

clear pctMove move pctMissed
for isess = 1:numel(fns{1})
    clear x
    for i = 1:numel(fns)

        if ismember(isess,badsess) && i==2
            ds = 1;
        else
            ds = 2;
        end

        temp = load(fullfile(datapth,fns{i}{isess}));
        dat.(token{i}) = temp.me;

        medat = dat.(token{i}).data;
        % downsample
        for j = 1:numel(medat)
            medat{j} = downsample(medat{j},ds);
            n(j,i) = numel(medat{j});
        end


    end

    for i = 1:numel(fns)

        if ismember(isess,badsess) && i==2
            ds = 1;
        else
            ds = 2;
        end

        temp = load(fullfile(datapth,fns{i}{isess}));
        dat.(token{i}) = temp.me;

        medat = dat.(token{i}).data;
        % downsample
        for j = 1:numel(medat)
            medat{j} = downsample(medat{j},ds);
            if n(j,1)~=n(j,2)
                if n(j,1) > n(j,2)
                    medat{j} = medat{j}(1:n(j,2));
                end
                if n(j,1) < n(j,2)
                    medat{j} = medat{j}(1:n(j,1));
                end
            end
        end

        thresh = dat.(token{i}).moveThresh;



        nSamp(i) = numel(cell2mat(medat'));

        move{i} = cellfun(@(x) x > thresh, medat, 'UniformOutput',false);


        nMove = cellfun(@(x) sum(x > thresh), medat, 'UniformOutput',false);
        nMove = sum([nMove{:}]');

        pctMove(isess,i) = nMove / nSamp(i) * 100;

        x(:,i) = cell2mat(move{i}')';


    end

    % n1(isess,1) = numel(move{1}{1});
    % n1(isess,2) = numel(move{2}{1});

    moveboth = x(:,1) | x(:,2);
    nMoveBoth = sum(moveboth);

    pctMove(isess,3) = nMoveBoth / nSamp(1) * 100;

end

%%

close all 
clear temp

f=figure;
f.Position = [680   715   318   263];
ax = gca;
hold on;

cols{1} = [0.55 0.55 0.55];
cols{2} = [0.3 0.3 0.3];
cols{3} = [0.7 0.7 0.7];

dat = pctMove;
dat(:,2) = pctMove(:,3);
dat(:,3) = pctMove(:,2);

xs = 1:3;
for i = 1:numel(xs)
    temp = dat(:,i);
    mu = mean(temp);
    b(i) = bar(xs(i),mu);
    b(i).FaceColor = cols{i};
    b(i).EdgeColor = 'none';
    b(i).FaceAlpha = 1;
    b(i).BarWidth = 0.7;
    vs(i) = scatter(xs(i)*ones(size(temp)),temp,30,'MarkerFaceColor','k',...
        'MarkerEdgeColor','w','LineWidth',1);
    errorbar(b(i).XEndPoints,mu,nanstd(temp),'LineStyle','none','Color','k','LineWidth',1)

    xs_(:,i) = vs(i).XData';
    ys_(:,i) = vs(i).YData';
end

plot([xs_(:,1) xs_(:,2)]',[ys_(:,1) ys_(:,2)]','k-','LineWidth',0.1)
plot([xs_(:,2) xs_(:,3)]',[ys_(:,2) ys_(:,3)]','k-','LineWidth',0.1)




ax.XTick = xs;
xticklabels(["Side cam" "Both" "Bottom cam" ])
ylabel("Time spent moving in session (%)")
ax = gca;
ax.FontSize = 12;






















































