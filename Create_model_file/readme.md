Readme file for creating model files needed in SeisTomoPy 
 From S. Durand March 2017

If you want to use SeisTomoPy with your model you must first create the appropriate model file. To do so you need to generate 2891 files named with the following format vs.DEPTH.xyz, DEPTH ranging from 2890 to 0. Every file must contain the model with this format:
lon     lat      dvs/vs (%)
where lon ranges from 0 to 359 with  a step of 1 degree
and lat ranges from -89 to 89 with  a step of 1 degree

1) Copy all your model files vs.DEPTH.xyz in the directory  model

3) Go to directory src/ and run make clean, make and ./main. This will create two model files YOURMODEL_1km.sph and YOURMODEL_5km.sph that you will need in SeisTomoPy. This can take several minutes to be done.

4) Copy these two model files in ../SeisTomoPy/fortran_files/models
