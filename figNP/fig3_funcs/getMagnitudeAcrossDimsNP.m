function mag = getMagnitudeAcrossDimsNP(rez,sm)


for i = 1:numel(rez) % for each session
    potent = rez(i).N_potent.^2;
    potent = sum(potent,3); % squared activity in each dimension, summed
    potent = mySmooth(potent,sm,'reflect');
%     mag(i).potent = potent ./ max(potent);
%     mag(i).potent = potent ./ mean(potent(115,:),2);
    mag(i).potent = normalize(potent,'zscore');


    null = rez(i).N_null.^2;
    null = sum(null,3);
    null = mySmooth(null,sm+5,'reflect');
%     mag(i).null = null ./ max(null);
%     mag(i).null = null ./ mean(null(115,:),2);
    mag(i).null = normalize(null,'zscore');

end

end