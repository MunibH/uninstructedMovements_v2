% plot cds in null and potent spaces in a heatmap

null = squeeze(cd_null_all.cd_proj(:,[1,2],1,:));
sel.null = squeeze(null(:,1,:) - null(:,2,:));
% for i = 1:size(null,3)
%     temp = null(:,:,i);
%     [~,pref] = max(temp');
%     nonpref = pref;
%     temp = nonpref;
%     nonpref(temp==1) = 2;
%     nonpref(temp==2) = 1;
%     for t = 1:size(null,1)
%         sel.null(t,i) = squeeze(null(t,pref(t),i) - null(t,nonpref(t),i));
%     end
% end

sel.null = sel.null ./ max(sel.null);
% sel.null = sel.null.^2;

potent = squeeze(cd_potent_all.cd_proj(:,[1,2],1,:));
sel.potent = squeeze(potent(:,1,:) - potent(:,2,:));
% for i = 1:size(potent,3)
%     temp = potent(:,:,i);
%     [~,pref] = max(temp');
%     nonpref = pref;
%     temp = nonpref;
%     nonpref(temp==1) = 2;
%     nonpref(temp==2) = 1;
%     for t = 1:size(potent,1)
%         sel.potent(t,i) = squeeze(potent(t,pref(t),i) - potent(t,nonpref(t),i));
%     end
% end

sel.potent = sel.potent ./ max(sel.potent);
% sel.potent = sel.potent.^2;





close all

xlims = [-2.4 2];
d = [-0.5 0];
for i = 1:numel(xlims)
    [~,ix(i)] = min(abs(obj(1).time - d(i)));
end
ix = ix(1):ix(2);

[~,n] = sort(mean(sel.null(ix,:),1),'ascend');
[~,p] = sort(mean(sel.potent(ix,:),1),'ascend');

sel.null = sel.null(:,n);
% sel.null(:,1)  = sel.null(:,1) ./ 1.5;
sel.potent = sel.potent(:,p);
% sel.potent(:,1) = sel.potent(:,1) ./ 2;


f = figure;

ax = nexttile;
hold on;
imagesc(obj(1).time,1:size(sel.null,2),sel.null');
[h,cmap] = colorbarpwn(min(min(sel.null)),max(max(sel.null)));
colormap(cmap);
colorbar;
xlim([xlims])
ylim([1,size(sel.null,2)])

ax = nexttile;
hold on;
imagesc(obj(1).time,1:size(sel.null,2),sel.potent');
% [h,cmap] = colorbarpwn(0,max(max(potent)));
[h,cmap] = colorbarpwn(min(min(sel.null)),max(max(sel.null)));
colormap(cmap);
colorbar;
xlim([xlims])
ylim([1,size(sel.null,2)])






