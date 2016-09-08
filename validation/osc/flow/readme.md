These simulations are used to compare the different implemenations to solve flow over a cylinder that oscilates on the x axis.  
The simulations are organized and compared on a graph sorted as follows:  
a b c  
d e f  
g h i  
j k l  
The first row is a coarse grid(h = 0.0625) with a low maximum CFL (0.35).  
The bottom three rows have a finer grid(h=0.3125) with an increasing maximum CFL(0.35, 0.7, 1).  
The first column is non-iteratively solved (GN and HN interpolaitons done outside of the linear solve) with an alpha weighting value calculated as discussed in luo et al.  
The second column is iteratively solved with alpha set to 1. This means it is being solved only using the interpolation solutions.  
The third column is the full luo et al. method with correct alpha values and iterative solves.  
