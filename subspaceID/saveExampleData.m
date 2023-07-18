exampleData.time = obj.time;
exampleData.seq = permute(obj.trialdat,[1 3 2]);
exampleData.motionEnergy = me.data;
exampleData.moveMask = me.data >= me.moveThresh;
exampleData.rightCorrectTrials = params.trialid{2};
exampleData.leftCorrectTrials = params.trialid{3};

save("exampleData.mat",'exampleData')