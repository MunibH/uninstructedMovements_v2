exampleData.anm = meta.anm;
exampleData.date = meta.date;
exampleData.time = obj.time;
exampleData.seq = permute(obj.trialdat,[1 3 2]);
exampleData.motionEnergy = me.data;
exampleData.moveMask = me.data >= me.moveThresh;

save("exampleData.mat",'exampleData')

%%