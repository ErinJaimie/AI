# Artificial Intelligence - Assignment 1

## Contributors:
Erin van den Heever, 22547134  

## src Files:
* assignment1.java: contains an implementation of the global best PSO simulation.
* bechmark.java: contains the suite of benchmark functions used for testing the PSO implementation.
* particle.java: contains the implementation of a particle data type.

## Executing the Code:
To run the PSO simulation for a specific benchmark function, open assignment1.java and uncomment the benchmark function that you intend to use.  
$ javac assignment1.java  
$ java assignment 1  

## Data:
* assignment1.java will generate .csv files containing the mean and standard deviation at each iteration for 4 performance measures: swarm diversity, solution quality, percentage of particles to leave the search space and average velocity magnitude.  
* additionally, assignment1.java outputs the global best position found by PSO after 5000 iterations for each independent run of the simulation. The simulation is run without velocity clamping and then with velocity clamping for k values of 0.1, 0.3 and 0.5.  
