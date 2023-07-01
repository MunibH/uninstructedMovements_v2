%%

% sess = [11 1 16];
sess = 1:numel(obj);
cond2use = [1:4];

cols = getColors;



%% potent

close all

col{1} = cols.rhit;
col{2} = cols.lhit;
col{3} = cols.rhit_aw;
col{4} = cols.lhit_aw;

ls{1} = '-';
ls{2} = '-';
ls{3} = '-';
ls{4} = '-';



f = figure;
for i = 1:numel(sess)
    ax = nexttile;
    hold on;
    temp = squeeze(cd_potent_all.cd_proj(:,cond2use,2,sess(i)));
    for c = 1:numel(cond2use)
        if c > 2
            sm = 11;
        else
            sm = 3;
        end
        plot(obj(1).time,mySmooth(temp(:,c),sm),'linewidth',2,'color',col{c},'linestyle',ls{c});
    end
    xlim([0 1])
    title([meta(sess(i)).anm ' ' meta(sess(i)).date])
end

%% null

close all

col{1} = cols.rhit;
col{2} = cols.lhit;
col{3} = cols.rhit_aw;
col{4} = cols.lhit_aw;

ls{1} = '-';
ls{2} = '-';
ls{3} = '-';
ls{4} = '-';



f = figure;
for i = 1:numel(sess)
    ax = nexttile;
    hold on;
    temp = squeeze(cd_null_all.cd_proj(:,cond2use,2,sess(i)));
    for c = 1:numel(cond2use)
        if c > 2
            sm = 11;
        else
            sm = 3;
        end
        plot(obj(1).time,mySmooth(temp(:,c),sm),'linewidth',2,'color',col{c},'linestyle',ls{c});
    end
    xlim([0 1])
    title([meta(sess(i)).anm ' ' meta(sess(i)).date])
end
