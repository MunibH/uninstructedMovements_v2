function rez = var_exp_NP_recon(input_data,rez,cond2use,params, trials, me)

% reconstruct single trial neural activity from null and potent spaces
% get R2 b/w reconstruction and true data

clear dat

mask = me.move;
mask = mask(:); 

d = size(input_data);
dat.full = reshape(input_data,d(1)*d(2),d(3)); % (time*trials,clu) 

dat.null = dat.full(~mask,:);
dat.potent = dat.full(mask,:);

scores = rez.N_null;
Q = rez.Qnull;
recon.null = tensorprod(scores,Q,3,2); % (time, trials, neurons)
d = size(recon.null);
null_reshape = reshape(recon.null,d(1)*d(2),d(3));

temp = corrcoef(dat.full,null_reshape);
rez.ve_recon.null_total = temp(1,2).^2;


scores = rez.N_potent; % (time,trials,dims)
Q = rez.Qpotent; % (neurons,dims)
recon.potent = tensorprod(scores,Q,3,2); % (time, trials, neurons)
d = size(recon.potent);
potent_reshape = reshape(recon.potent,d(1)*d(2),d(3));

temp = corrcoef(dat.full,potent_reshape);
rez.ve_recon.potent_total = temp(1,2).^2;

% for cluix = 1:size(input_data,3)
%     dat = input_data(:,:,cluix);
%     pred = recon.null(:,:,cluix);
%     rez.ve_recon.null(cluix) = getR2(dat(:),pred(:));
% end
% for cluix = 1:size(input_data,3)
%     dat = input_data(:,:,cluix);
%     pred = recon.potent(:,:,cluix);
%     rez.ve_recon.potent(cluix) = getR2(dat(:),pred(:));
% end

rez.recon.null = recon.null;
rez.recon.potent = recon.potent;

% reconstruct PSTHs
for i = 1:numel(cond2use)
    trix = params.trialid{cond2use(i)};

    rez.recon_psth.null(:,:,i) = squeeze(mean(recon.null(:,trix,:),2));
    rez.recon_psth.potent(:,:,i) = squeeze(mean(recon.potent(:,trix,:),2));
end


% reconstruct single trial neural activity from null and potent spaces
% get R2 b/w reconstruction and true data
clear dat recon scores

mask = me.move;
mask = mask(:); 

d = size(input_data);
dat.full = reshape(input_data,d(1)*d(2),d(3)); % (time*trials,clu) 

dat.null = dat.full(~mask,:);
dat.potent = dat.full(mask,:);

% 
scores.null = rez.N_null;
d = size(scores.null);
scores.null = reshape(scores.null,d(1)*d(2),d(3));

scores.null_prep = scores.null(~mask,:);
scores.null_move = scores.null(mask,:);

scores.potent = rez.N_potent;
d = size(scores.potent);
scores.potent = reshape(scores.potent,d(1)*d(2),d(3));

scores.potent_prep = scores.potent(~mask,:);
scores.potent_move = scores.potent(mask,:);

%

Q = rez.Qnull;
recon.null_prep = scores.null_prep * Q';
recon.null_move = scores.null_move * Q';

Q = rez.Qpotent;
recon.potent_prep = scores.potent_prep * Q';
recon.potent_move = scores.potent_move * Q';

%
temp = corrcoef(dat.null,recon.null_prep);
rez.ve_recon.null_prep = temp(1,2).^2;

temp = corrcoef(dat.potent,recon.null_move);
rez.ve_recon.null_move = temp(1,2).^2;

temp = corrcoef(dat.null,recon.potent_prep);
rez.ve_recon.potent_prep = temp(1,2).^2;

temp = corrcoef(dat.potent,recon.potent_move);
rez.ve_recon.potent_move = temp(1,2).^2;



end
