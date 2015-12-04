# Bayasian Sequential Gaussian Simulator (BSGS) #
Bayesian Sequential Gaussian Simulation purpose is to use a secondary variable $Z$ to guide a sequential simulation of a primary variable $X$. 

## History and Developpement ##
* The Bayasian Sequential Gaussian Simulation technique was introduced by  Doyen & Boer in 1996 to extrapolate non linear data of lithology. 
* Ruggeri developped a two-step approch of the inital BSGS to combine electrical and hydrological conductivity.

![CharFlow_BSGS_resized.png](https://bitbucket.org/repo/gABK6j/images/26574965-CharFlow_BSGS_resized.png) ![BSGS_Simulation_2.png](https://bitbucket.org/repo/gABK6j/images/3694505833-BSGS_Simulation_2.png)

## Theorie of BSGS ##
The overall objectif of BSGS is to generate a realisations of a primary variable using two dataset : (1) a sparce high resolution primary and secondary variable at identical location and  (2) low-resolution secondary variable covering the entire region.

The method estimate all unsampled location as follow :
1. Randomly choose a point
2. Using kriging technique find the mean and standard deviation determining the prior estimatet
3. Find the Baysian Likelihood at the location from collocated data
4. Estimate the posteriori distribution using Baysian theorem
5. Sample form this distribution to determine the value.
![BSGS_Simulation.png](https://bitbucket.org/repo/gABK6j/images/3538075189-BSGS_Simulation.png)


## How to use it ? ##

Several feathure can be turn on/off such as :  variogramm fitting to conditional data, smart neighbooring search for krigeage, constant path and weight for kriging, n-score transform. This feathure are made to make the code running faster.


## About the Developper team ##

*Current developer : RaphaÃ«l Nussbaumer rafnuss@gmail.com*
