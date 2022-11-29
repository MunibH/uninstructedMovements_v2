
sessix = 2;

figure;
for ifeat = 1:size(kin(sessix).dat,3)
    clf
    imagesc(kin(sessix).dat(200:1600,:,ifeat)')
    colorbar
    title(kin(sessix).featLeg{ifeat},'Interpreter','none')
    pause
end