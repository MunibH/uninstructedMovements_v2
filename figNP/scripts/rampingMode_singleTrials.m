rt = firstJawRT(obj);

%%

close all

% ramping mode single trials

rampix = find(strcmpi(cd_null(1).cd_labels,'ramping'));

f = figure;
for sessix = 1:numel(obj)
    ax = nexttile; hold on;

    edges = [-0.8 , 0];
    for i = 1:numel(edges)
        [~,t(i)]  = min(abs(obj(sessix).time - edges(i)));
    end

    scores = mySmooth(cd_potent(sessix).trialdat(t(1):t(2),:,rampix),21);

%     shadedErrorBar(obj(1).time,mean(scores,2),std(scores,[],2)./sqrt(250),{'Color','k'},0.2,ax)


%     gradscores = gradient(scores',1./params(1).dt);

%     scores = abs(mean(gradient(,1),1));
    
    thisrt = rt{sessix};
% 
%     scatter(thisrt,scores,10,'filled','Color','k')
    temp = repmat(thisrt,30,1);
    toplot = cat(1,scores,temp);

    imagesc(toplot')
    colormap(linspecer)
    colorbar
end