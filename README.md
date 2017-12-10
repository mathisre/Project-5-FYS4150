# Project-5-FYS4150

The folder programs contains all the programs used in the project. Some alterations have been made to the classes outside of what was necessary. For example the IO class was changed to use string filenames instead of const char. The main programs are MD_many_temperatures.cpp and MD_one_temperature.cpp. 
The many temperature program runs the molecular dynamics for 10000 time steps for a range of temperatures to be determined from the command line. It then saves the time, temperature and diffusion constant to a text file and makes a movie file for each initial temperature. The movie files contain fewer time steps to save space. 
The single tmeperature program runs the molecular dynamics for a single temperature and it saves the movie file for each timestep. 
The many temperatures program runs with three command line arguments: initial temperature, final temperature and temperature step size. The single temperature program runs with two command line arguments: inintial temperature and number of time steps.

MD_plot.py is used to create the plots used in the project. The program runs by reading directory containing temp_diff_x.txt files, start x, final x and x spacing. It then produces the 4 graphs used in the project. For example it can be run as: python MD_plot.py /Data/ 100 900 10

The Data folder contains two example argon_xx.xyz files. These can be used in Ovito to see how the atoms move in the program.
