% subspace id schematic

% load some data using a script first


%% motion energy + select 99th percentile

pth.vid{1} = 'C:\Users\munib\Documents\Economo-Lab\data\Video\JEB15\2022-07-27\Cam0';
pth.vid{2} = 'C:\Users\munib\Documents\Economo-Lab\data\Video\JEB15\2022-07-27\Cam1';
energy = getME(pth);

%%
close all

f = figure;
ax = gca;
f.Renderer = 'painters';

e = energy(:,:,20);
e(140:180,1:135) = 0;
e(200:220,115:130) = 0;
e(160:195,130:165) = 0;
% e(145:152,134:139) = 0;
e = imgaussfilt(e,1.5);
e = int32(e);
im = imagesc(e);
cmap = colormap(linspecer);
RGB = ind2rgb(e, cmap);
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off';

% imwrite(RGB,'me.tif','tif')

%%
close all

temp = energy(:,:,20);

f = figure;
ax = gca;
hold on;
f.Renderer = 'painters';
ax = prettifyPlot(ax);
hh = histogram(temp(:),100,'Normalization','pdf');

mdl = fitgmdist(temp(:),1);
x = hh.BinEdges';
p = pdf(mdl,x);
plot(x,p,'k','linewidth',2)
xlabel('Single frame motion energy')
ylabel('Probability')
delete(hh);
xlim([0 6])

prctile(p,99)

%% threshold motion energy
close all

allme = me.data(:);
f = figure;
ax = gca;
hold on;
f.Renderer = 'painters';
ax = prettifyPlot(ax);
hh = histogram(allme,1000,'Normalization','probability');


mdl = fitgmdist(allme,3);
x = hh.BinEdges';
p = pdf(mdl,x);
plot(x,p,'k','linewidth',2)
xlabel('Motion energy')
ylabel('Probability')
% delete(hh);
xlim([0 50])

%% label move and non-move time points (ME and neural)

plotData_v2(obj,params,me,0.5)


%% covariance matrices

%% subspace id


function energy = getME(pth)
%%
% pth = load('pth.mat');
% pth = pth.pth;

params.view = 1; % side cam
params.pre = 5;  % num frames prior current frame to use for ME calc
params.post = 5; % num frames post current frame to use for ME calc
params.percentile = 99; % percentile ME values to store - this is to avoid ME values being biased by pixels with no animal or nothing happening

% crop out y pixels below jaw
params.xcrop = [0.01 1]; % keep the x pixels from (1/4)x:(3/4)x
params.ycrop = [0.01 1]; % keep the y pixels from (1/4)y:(3/4)y

% get vid files
vidfn = cell(1,numel(pth.vid));
for i = 1:numel(pth.vid)
    contents = dir(fullfile(pth.vid{i}, '*.avi'));
    vidfn{i} = {contents.name}';
    for j = 1:numel(vidfn{i})
        vidfn{i}{j} = fullfile(pth.vid{i},vidfn{i}{j});
    end
end


% motion energy for each trial
% disp(['----- ' pth.anm ' ' pth.dt ' -----']);
% me = getMotionEnergy(vidfn{params.view},params); % (anm,date)

vidfn = vidfn{params.view};

me = cell(numel(vidfn),1); % (trials,1) cell array
for trix = 10:numel(vidfn)
    disp(['Processing ME for trial  ' num2str(trix) '/' num2str(numel(vidfn))])

    % read video file for current trial
    v = VideoReader(vidfn{trix});

    % all videos are same size, so only need to do this once
    xix = round(params.xcrop(1).*v.Width):round(params.xcrop(2).*v.Width);
    yix = round(params.ycrop(1).*v.Height):round(params.ycrop(2).*v.Height);

    cnt = 0;
    mov = zeros(numel(yix), numel(xix), v.Duration);
    while hasFrame(v)
        cnt = cnt+1;
        frame = readFrame(v);
        mov(:,:,cnt) = frame(yix,xix,1);
    end

    % Motion energy for current trial,
    energy = zeros(size(mov));
    for i = 1:size(mov,3)
        if (i > (params.pre+1)) && (i < (size(mov,3)-params.post))
            energy(:,:,i) = computeME(mov,i,params);
        else
            continue
        end
    end

    % for each trial, keep 99th percentile values of motion energy for each
    % frame
    % ereshape = reshape(energy,size(energy,1)*size(energy,2),size(energy,3));
    % me{trix} = prctile(ereshape,params.percentile,1);
    % clear mov energy ereshape

    if trix == 20
        return
    end
end



end


function e = computeME(mov,i,params)
e = abs( median(mov(:,:,i:i+params.post),3) - median(mov(:,:,i-params.pre:i),3)  );
end

function plotData_v2(obj,params,me,dy)
ts = mode(obj.bp.ev.bitStart) - 2.5;
sample = mode(obj.bp.ev.sample) - 2.5;
delay = mode(obj.bp.ev.delay) - 2.5;
gc = 0;

nsm = 5;
msm = 7;
nclu = 5;

cols = getColors;
lw = 2;
alph = 0.23;
alph2 = 0.13;
cscale = 0.8;

% while true
for trix = 11%1:obj.bp.Ntrials
    if obj.bp.early(trix)
        continue
    end
    clu = randsample(size(obj.psth,2),nclu);
    %     trix = randsample(obj.bp.Ntrials,1);
    %     trix = 139; [58 164 170]
    clu(1) = 36;
    clu(2) = 23;
    clu(3) = 22;
    clu(4) = 51;
    clu(5) = 47;

    me = stitchAndPurgeMoveBouts(me,params,0.05,0.01);

    m = mySmooth(normalize(me.data(:,trix),'range',[0 2]),msm);
    move = me.move(:,trix);
    mnull = m;
    mpotent = m;
    mnull(~move,:) = nan;
    mpotent(move,:) = nan;

    f = figure;
    f.Position = [582   646   635   144];
    f.Renderer = "painters";
    ax = gca;
    ax = prettifyPlot(ax);
    title(num2str(trix))
    hold on;

    plot(obj.time,m+dy,'Color','k','LineWidth',lw)
    % plot(obj.time,mpotent,'Color','k','LineWidth',lw) % actually null data


    xline(ts,'k--'); xline(delay,'k--'); xline(sample,'k--'); xline(gc,'k--')
    xlim([-2.35 2])

    y = ax.YLim;
    Y = [ax.YLim(1) ax.YLim(1) ax.YLim(2) ax.YLim(2)];

    [nstart, nend, nlen] = ZeroOnesCount(~move);
    nend = nend + 1;
    try
        if nend(end) > numel(obj.time)
            nend(end) = numel(obj.time);
        end
    catch
        delete(f)
        continue
    end
    ns = obj.time(nstart);
    ne = obj.time(nend);
    for ii = 1:numel(nstart)
        ys = mpotent(nstart(ii):nend(ii));
        xs = linspace(ns(ii),ne(ii),numel(ys));
        ff(ii) = area(xs,ys+dy);
        ff(ii).FaceAlpha = 1;
        ff(ii).EdgeColor = 'none';
        ff(ii).FaceColor = cols.null;
        uistack(ff(ii),'bottom');
    end
    [mstart, mend, mlen] = ZeroOnesCount(move);
    mstart = mstart - 1;
    if mstart(1) < 1
        mstart(1) = 1;
    end
    mend = mend + 1;
    if mend(end) > numel(obj.time)
        mend(end) = numel(obj.time);
    end
    ms = obj.time(mstart);
    me_ = obj.time(mend);
    for ii = 1:numel(mstart)
        ys = m(mstart(ii):mend(ii));
        xs = linspace(ms(ii),me_(ii),numel(ys));
        ff(ii) = area(xs,ys+dy);
        ff(ii).FaceAlpha = 1;
        ff(ii).EdgeColor = 'none';
        ff(ii).FaceColor = cols.potent;
        uistack(ff(ii),'bottom');
    end

    % for i = 1:nclu
    %     n{i} = mySmooth(normalize(squeeze(obj.trialdat(:,clu(i),trix)),'range',[0 1]),nsm);
    %     nnull{i} = n{i};
    %     npotent{i} = n{i};
    %     nnull{i}(~move,:) = nan;
    %     npotent{i}(move,:) = nan;
    % end
    % 
    % f = figure;
    % f.Position = [437   139   498   659];
    % f.Renderer = "painters";
    % ax = gca;
    % ax = prettifyPlot(ax);
    % hold on;
    % 
    % cc = cols.potent * cscale;
    % plot(obj.time,m,'Color',cc,'LineWidth',lw)
    % cc = cols.null * cscale;
    % plot(obj.time,mpotent,'Color',cc,'LineWidth',lw)
    % %     for i = 1:nclu
    % %         plot(obj.time,n{i}-(1.5*i),'Color',cols.potent,'LineWidth',lw)
    % %         plot(obj.time,npotent{i}-(1.5*i),'Color',cols.null,'LineWidth',lw)
    % %     end
    % %     plot(obj.time,m,'Color','k','LineWidth',lw)
    % %     plot(obj.time,mpotent,'Color','k','LineWidth',lw)
    % for i = 1:nclu
    %     plot(obj.time,n{i}-(1.5*i),'Color','k','LineWidth',lw)
    %     plot(obj.time,npotent{i}-(1.5*i),'Color','k','LineWidth',lw)
    % end
    % title(['Trial ' num2str(trix) ' - Clu ' num2str(clu')])
    % xline(ts,'k--'); xline(delay,'k--'); xline(sample,'k--'); xline(gc,'k--')
    % xlim([-2.35 2])
    % 
    % y = ax.YLim;
    % Y = [ax.YLim(1) ax.YLim(1) ax.YLim(2) ax.YLim(2)];
    % 
    % [nstart, nend, nlen] = ZeroOnesCount(~move);
    % nend = nend + 1;
    % try
    %     if nend(end) > numel(obj.time)
    %         nend(end) = nend(end) - 1;
    %     end
    % catch
    %     delete(f);
    %     continue
    % end
    % ns = obj.time(nstart);
    % ne = obj.time(nend);
    % for ii = 1:numel(nstart)
    %     ff(ii) = fill([ns(ii) ne(ii) ne(ii) ns(ii)], Y, cols.null);
    %     ff(ii).FaceAlpha = alph;
    %     ff(ii).EdgeColor = 'none';
    %     uistack(ff(ii),'bottom');
    %     ff2(ii) = fill([ns(ii) ne(ii) ne(ii) ns(ii)], [-0.3 -0.3 ax.YLim(2) ax.YLim(2)], cols.null);
    %     ff2(ii).FaceAlpha = alph2;
    %     ff2(ii).EdgeColor = 'none';
    % end
    % [mstart, mend, mlen] = ZeroOnesCount(move);
    % mend = mend + 1;
    % if mend(end) > numel(obj.time)
    %     mend(end) = mend(end) - 1;
    % end
    % ms = obj.time(mstart);
    % me_ = obj.time(mend);
    % for ii = 1:numel(mstart)
    %     ff(ii) = fill([ms(ii) me_(ii) me_(ii) ms(ii)], Y, cols.potent);
    %     ff(ii).FaceAlpha = alph;
    %     ff(ii).EdgeColor = 'none';
    %     uistack(ff(ii),'bottom');
    %     ff2(ii) = fill([ms(ii) me_(ii) me_(ii) ms(ii)], [-0.3 -0.3 ax.YLim(2) ax.YLim(2)], cols.potent);
    %     ff2(ii).FaceAlpha = alph2;
    %     ff2(ii).EdgeColor = 'none';
    %     %         uistack(ff2(ii),'bottom');
    % end
    % %     uistack(ff,'bottom')
    % ax.Visible = 'off';
    % %     break
    % pause
    % delete(f)
end
end

function plotData(obj,params,me)
ts = mode(obj.bp.ev.bitStart) - 2.5;
sample = mode(obj.bp.ev.sample) - 2.5;
delay = mode(obj.bp.ev.delay) - 2.5;
gc = 0;

nsm = 5;
msm = 7;
nclu = 5;

cols = getColors;
lw = 2;
alph = 0.23;
alph2 = 0.13;
cscale = 0.8;

% while true
for trix = 1:obj.bp.Ntrials
    if obj.bp.early(trix)
        continue
    end
    clu = randsample(size(obj.psth,2),nclu);
    %     trix = randsample(obj.bp.Ntrials,1);
    %     trix = 139; [58 164 170]
    clu(1) = 36;
    clu(2) = 23;
    clu(3) = 22;
    clu(4) = 51;
    clu(5) = 47;

    me = stitchAndPurgeMoveBouts(me,params,0.05,0.01);

    m = mySmooth(normalize(me.data(:,trix),'range',[0 2]),msm);
    move = me.move(:,trix);
    mnull = m;
    mpotent = m;
    mnull(~move,:) = nan;
    mpotent(move,:) = nan;


    for i = 1:nclu
        n{i} = mySmooth(normalize(squeeze(obj.trialdat(:,clu(i),trix)),'range',[0 1]),nsm);
        nnull{i} = n{i};
        npotent{i} = n{i};
        nnull{i}(~move,:) = nan;
        npotent{i}(move,:) = nan;
    end

    f = figure;
    f.Position = [437   139   498   659];
    f.Renderer = "painters";
    ax = gca;
    ax = prettifyPlot(ax);
    hold on;

    cc = cols.potent * cscale;
    plot(obj.time,m,'Color',cc,'LineWidth',lw)
    cc = cols.null * cscale;
    plot(obj.time,mpotent,'Color',cc,'LineWidth',lw)
    %     for i = 1:nclu
    %         plot(obj.time,n{i}-(1.5*i),'Color',cols.potent,'LineWidth',lw)
    %         plot(obj.time,npotent{i}-(1.5*i),'Color',cols.null,'LineWidth',lw)
    %     end
    %     plot(obj.time,m,'Color','k','LineWidth',lw)
    %     plot(obj.time,mpotent,'Color','k','LineWidth',lw)
    for i = 1:nclu
        plot(obj.time,n{i}-(1.5*i),'Color','k','LineWidth',lw)
        plot(obj.time,npotent{i}-(1.5*i),'Color','k','LineWidth',lw)
    end
    title(['Trial ' num2str(trix) ' - Clu ' num2str(clu')])
    xline(ts,'k--'); xline(delay,'k--'); xline(sample,'k--'); xline(gc,'k--')
    xlim([-2.35 2])

    y = ax.YLim;
    Y = [ax.YLim(1) ax.YLim(1) ax.YLim(2) ax.YLim(2)];

    [nstart, nend, nlen] = ZeroOnesCount(~move);
    nend = nend + 1;
    try
        if nend(end) > numel(obj.time)
            nend(end) = nend(end) - 1;
        end
    catch
        delete(f);
        continue
    end
    ns = obj.time(nstart);
    ne = obj.time(nend);
    for ii = 1:numel(nstart)
        ff(ii) = fill([ns(ii) ne(ii) ne(ii) ns(ii)], Y, cols.null);
        ff(ii).FaceAlpha = alph;
        ff(ii).EdgeColor = 'none';
        uistack(ff(ii),'bottom');
        ff2(ii) = fill([ns(ii) ne(ii) ne(ii) ns(ii)], [-0.3 -0.3 ax.YLim(2) ax.YLim(2)], cols.null);
        ff2(ii).FaceAlpha = alph2;
        ff2(ii).EdgeColor = 'none';
    end
    [mstart, mend, mlen] = ZeroOnesCount(move);
    mend = mend + 1;
    if mend(end) > numel(obj.time)
        mend(end) = mend(end) - 1;
    end
    ms = obj.time(mstart);
    me_ = obj.time(mend);
    for ii = 1:numel(mstart)
        ff(ii) = fill([ms(ii) me_(ii) me_(ii) ms(ii)], Y, cols.potent);
        ff(ii).FaceAlpha = alph;
        ff(ii).EdgeColor = 'none';
        uistack(ff(ii),'bottom');
        ff2(ii) = fill([ms(ii) me_(ii) me_(ii) ms(ii)], [-0.3 -0.3 ax.YLim(2) ax.YLim(2)], cols.potent);
        ff2(ii).FaceAlpha = alph2;
        ff2(ii).EdgeColor = 'none';
        %         uistack(ff2(ii),'bottom');
    end
    %     uistack(ff,'bottom')
    ax.Visible = 'off';
    %     break
    pause
    delete(f)
end
end
