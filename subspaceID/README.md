# Matlab code for null and potent subspace identification 

Version 1.0  July 19, 2023

@ 2023 Munib Hasnain   munibh@bu.edu | Jackie Birnbaum   jackieb1@bu.edu 

## VERSION HISTORY (will go here)


## HOW TO GET STARTED

See `main.m`

1) Single trial neural firing rates should be formatted as `(number of time bins, number of trials, number of neurons)`

2) You need an binary annotation of when the animal is moving (`1`) or stationary (`0`). This matrix will be of time `(number of time bins, number of trials)`

3) The optimization depends on the manopt toolbox. This github repo contains the latest version as of July 19, 2023. Manopt can be downloaded [here](https://www.manopt.org/). 

3) That's all you need! Now you can run the the example script (`main.m`)

4) You can run through `main.m` using the data from an example session in `exampleData.mat` 