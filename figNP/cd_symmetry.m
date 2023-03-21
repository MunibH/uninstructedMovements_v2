close all

cond2use = [2 3];

f = figure;
for sessix = 1:numel(meta)
    for icond = 1:numel(cond2use)
        trix{icond} = params(sessix).trialid{cond2use(icond)};
        me_{icond} = me(sessix).data(:,trix{icond});
    end

    ax = nexttile;
    hold on;
    plot(obj(sessix).time,mean(me_{1},2),"Color",'b')
    plot(obj(sessix).time,mean(me_{2},2),"Color",'r')
    title([meta(sessix).anm ' ' meta(sessix).date])
end
sgtitle('motion energy')

f = figure;
for sessix = 1:numel(meta)
    cd.null = cd_null(sessix).cd_proj(:,[1 2],1);
    
    ax = nexttile;
    hold on;
    plot(obj(sessix).time,cd.null(:,1),"Color",'b')
    plot(obj(sessix).time,cd.null(:,2),"Color",'r')
    title([meta(sessix).anm ' ' meta(sessix).date])
end
sgtitle('null')

f = figure;
for sessix = 1:numel(meta)
    cd.potent = cd_potent(sessix).cd_proj(:,[1 2],1);
    
    ax = nexttile;
    hold on;
    plot(obj(sessix).time,cd.potent(:,1),"Color",'b')
    plot(obj(sessix).time,cd.potent(:,2),"Color",'r')
    title([meta(sessix).anm ' ' meta(sessix).date])
end
sgtitle('potent')



%%
cond2use = [2 3];
for sessix = 1:numel(meta)
    
    for icond = 1:numel(cond2use)
        trix{icond} = params(sessix).trialid{cond2use(icond)};
        me_{icond} = me(sessix).data(:,trix{icond});
    end

    cd.null = cd_null(sessix).cd_proj(:,[1 2],1);
    cd.potent = cd_potent(sessix).cd_proj(:,[1 2],1);

    plotME(obj,me_)


    figure;
    ax = nexttile;
    hold on;
    plot(obj(sessix).time,cd.null(:,1),"Color",'b')
    plot(obj(sessix).time,cd.null(:,2),"Color",'r')
    title('null')

    ax = nexttile;
    hold on;
    plot(obj(sessix).time,cd.potent(:,1),"Color",'b')
    plot(obj(sessix).time,cd.potent(:,2),"Color",'r')
    title('potent')

    ax = nexttile;
    hold on;
    plot(obj(sessix).time,mean(me_{1},2),"Color",'b')
    plot(obj(sessix).time,mean(me_{2},2),"Color",'r')
    title('ME')

    sgtitle([meta(sessix).anm ' ' meta(sessix).date])

    pause
    close all


end