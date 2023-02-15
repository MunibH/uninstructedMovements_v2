obj = obj(1);
meta = meta(1);
params = params(1);
me = me(1);

%% single trial method
close all

sample = mode(obj.bp.ev.sample) - 2.5;
delay = mode(obj.bp.ev.delay) - 2.5;

clu = 13;
trix = 4;
sm = 11;

dat.m = mySmooth(me.data(:,trix),sm,'reflect');
dat.n = mySmooth(obj.trialdat(:,clu,trix),sm,'reflect');

thresh = me.moveThresh - 3;

dat.nullix = dat.m < thresh;
dat.potentix = dat.m >= thresh;

dat.mnull = dat.m;
dat.mpotent = dat.m;
dat.mpotent(dat.nullix) = nan;

dat.nnull = dat.n;
dat.npotent = dat.n;
dat.npotent(dat.nullix) = nan;



f = figure;
cols = getColors;
subplot(2,1,1)
plot(obj.time,dat.mnull,'Color',cols.null,'LineWidth',2);
hold on;
plot(obj.time,dat.mpotent,'Color',cols.potent,'LineWidth',2);
xlim([-2.4,2])
xline(sample,'k--')
xline(delay,'k--')
xline(0,'k--')


subplot(2,1,2)
plot(obj.time,dat.nnull,'Color',cols.null,'LineWidth',2);
hold on;
plot(obj.time,dat.npotent,'Color',cols.potent,'LineWidth',2);
xlim([-2.4,2])
xline(sample,'k--')
xline(delay,'k--')
xline(0,'k--')



%% elsayed 2016 method
close all

sample = mode(obj.bp.ev.sample) - 2.5;
delay = mode(obj.bp.ev.delay) - 2.5;

clu = 13;
sm = 5;

dat.n = mySmooth(obj.psth(:,clu,1),sm,'reflect');

dat.nullix   = obj.time >= delay & obj.time <= 0;
dat.potentix = obj.time >= 0 & obj.time <= 1;

dat.nnull = dat.n;
dat.nnull(~dat.nullix) = nan;

dat.npotent = dat.n;
dat.npotent(~dat.potentix) = nan;




f = figure;
cols = getColors;
plot(obj.time,dat.n,'Color','k','LineWidth',2);
hold on;
plot(obj.time,dat.nnull,'Color',cols.null,'LineWidth',2);
plot(obj.time,dat.npotent,'Color',cols.potent,'LineWidth',2);
xlim([-2.4,2])
xline(sample,'k--')
xline(delay,'k--')
xline(0,'k--')





