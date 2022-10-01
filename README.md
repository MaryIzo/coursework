# MSU Coursework 2022
### Investigation of the effectiveness of various cluster methods for analyzing the structure of gamma-hadron families

##### Intoduction.    
We use the Pamir experiment as a data source in the work, the result of a joint study of several scientific institutions led by FIAN. It is the study of hadron interactions by the method of large X-ray emulsion chamber.    

##### Experiment.     
From the most general physical considerations follows that it is necessary to study interactions at the highest energies in order to study the details of the internal structure of hadrons.   
There are two ways to research high-energy particles – accelerators (like Large Hadron Collider) and cosmic rays. But the upper limit of the energy in cosmic radiation is very high. So the cosmic rays make it possible to study the hadron-nucleus and nucleus-nucleus interactions at very high energies.   
As I said before there is a method of X-ray emulsion chambers (REC) installed at mountain heights to study interactions at energies $ E_0 $> 100 TeV (larger). With the help of detectors with areas up to 1000 m $ ^ 2 $ (squared), statistics are accumulated throughout the year.   
The main object of research is groups of genetically related particles (gamma quanta and hadrons) with energy larger than 4 TeV, which arise as a result of the interaction of a primary high-energy particle with the nuclei of air atoms (nuclear-electromagnetic cascade). Such groups of particles we call families.  
There are lots of difficulties in the experiments with cosmic rays:
the low intensity of high-energy particles, 
the complex chemical composition of the primary cosmic radiation, 
the presence of cascade processes (electromagnetic and nuclear) in the atmosphere. 
That’s why we need modeling data and comparing the experiment with the results of this model.   
In the experiment, we study the electron beam of the electron-photon cascade, which leads to the appearance of dark spots on the X-ray film. In addition, it is possible to determine the incidence angle of the cascade.   
Also we know from calculations that events in the Pamir experiment are a mixture of secondary particles on average of two to three generations in the nuclear-electromagnetic cascade.    
##### Data format.   
As a data source, we used files of experimental (BC.DAT and BCR.DAT) and model (BC.MC0 and BCR.MC0) data in unformatted form.   
After reading, the data was converted to formats, txt, xlsx and csv.   
BC.DAT contains clustered data information. BCR.DAT contains information about all points, about 1000 events in total. Information includes abscissa (X), ordinate (Y) and energy value (E) for each spot.   
Similarly, for files containing simulated data, BC.MC0 contains information about clustered data, BCR.MC0 contains information about all points, about 1500 events in total.   
The information about point includes the abscissa (X), ordinate (Y), energy value (E), the height of occurrence for each spot (the height of the last interaction) (H), and the number of protons in the original particle (A0).   
The values of the X and Y coordinates of all blackening spots are measured from the energy center of gravity, calculated in accordance with the formula.    
With such initial data, the task is to compare the experimental and model data with each other and, by the cluster structure, find particles that came from the same height (the last interaction).    
Coordinates are measured from the energy center of gravity.   
##### Distributions.     
Since the clustering algorithm is selected on the model data and further applied to the experimental data, we need a verification of the homogeneity of the distributions of the experimental and model data.   
Radius is a Euclidean distance of a particle from energy center of gravity. K is the number of particles in the family.    
Let’s consider the distributions:    
total energy in family    
average radius      
average composition of radius and energy      
plurality    
total threshold energy     
energy    
radius     
composition of radius and energy.     
s – the number of threshold particles which satisfy the inequality:     
the ratio of particle energy to the sum of the energies of all the more preceding energetic particles is higher value 0.04.    
To compare the distributions, we determine the total energy for each family and all families in data is divided into intervals of total energies (100-200, 200-400, 400-700,> 700,> 100, 100-700).   
The homogeneity of the experimental and model data was tested using several tests:     
Smirnov (equality of distribution functions),     
Mann-Whitney (equality of distribution functions),     
Welch's t-student,     
Levene.     
Tests the null hypothesis of equality of means.     
student's criterion is implemented independently in Python.   For another tests was used a ready-made implementation.    
Lent’s consider distributions in which, as a result of the statistical solution, on the significance level α = 0.01 in different energy intervals, we can accept the hypothesis of the homogeneity.     
Let us consider a table with the numerical values of the calculated statistics (left column) and the critical value (right column) for some distributions in energy intervals.    
Distributions in which, as a result of the statistical solution, on the significance level α = 0.01 in different energy intervals, we can accept the hypothesis of the equal averages (left) and equal variances (right).     
this is the probability of how often the statistic of the test, given that the hypothesis H0 is correct, will take on a value (for other samples) greater than that obtained on the given sample.   
##### Clustering.     
In the process of work were considered 7 clustering algorithms with an unknown number of clusters. As a result 3 of them showed the best result.    
For each family, the parameter is chosen as the N-th percentile of the ER distribution (product of energy and radius).   
   
AGGLOMERATIVE CLUSTERING.     
Start with many small clusters and merge them together to create bigger clusters.   
1) Make each data point a cluster   
2) Take the two closest clusters and make them one cluster.  
3) Repeat step 2 until there is only one cluster
4) The linkage criteria refers to how the distance between clusters is calculated.
We tried to use average linkage:
(The distance between clusters is the average distance between each point in one cluster to every point in other cluster)
Similar to gradient descent, you can tweak certain parameters to get drastically different results.
  
DBSCAN    
There are two parameters to the algorithm,  minPoints  and  eps, which define formally what we mean when we say dense.   
eps: specifies how close points should be to each other to be considered a part of a cluster. It means that if the distance between two points is lower or equal to this value (eps), these points are considered neighbors.   
minPoints: the minimum number of points to form a dense region. For example, if we set the minPoints parameter as 4, then we need at least 4 points to form a dense region.    
  
ALGORITHM PAMIR      
… to do description …   

##### Coefficients
There are 4 coefficients of clustering quality for the model families.    
1) Purity P  = (the number of particles from the dominant interaction in the given cluster, arriving from the same height) / (the total number of particles in the considered cluster)     
The purity P is calculated for each selected cluster, and then averaged over all clusters in the given family.    
2) Splitting S = (the number of different dominant interactions contributing to the considered family) / (the number of clusters selected in the family)     
3) Integrity I  = (the number of different dominant interactions that contribute to the family under consideration) / (the total number of recent interactions that contribute to the family under consideration)     
4) Efficiency E  = (the number of particles from the dominant interaction in the given cluster that came from the same height) / (the number of particles in the entire family formed in the dominant interaction under consideration)     
Efficiency is calculated for each selected cluster, and then averaged over all clusters in a given family.     
The value of N varies from 5 to 100 with a step of 5; additionally consider the value 98.   
#####  Conclusions.     
As a result of the work, we investigated the bank of experimental and model data.    
Using the criterion of homogeneity of Smirnov, Mann-Whitney, Student's modification of Welch, and Levene, the homogeneity of experimental and model data. The result was obtained for the uniformity of the distributions of PE, R, ER over the energy intervals and PE, R, Nγ over the entire energy spectrum. Found and explained deviations in homogeneity in the intervals of the total energies of the families 100-200 TeV and > 700 TeV.    
Visualization of distributions (histograms and empirical functions) to compare experiment and model, check ranges, and find errors.    
Visualization of families, divided into groups of primary particles, 2d and 3d in png and html formats for comparing and understanding clustered data.    
A distance metric has been introduced, which is used for further clustering.    
The coefficients of the clustering quality have been determined.   
7 clustering methods were studied, among which 2 algorithms were selected. The selection of optimal parameters for the most efficient clustering of families has been started.    
Comparison clustering algorithms.     
