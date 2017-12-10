# Project-5-FYS4150

The folder programs contains all the programs used in the project. Some alterations have been made to the classes outside of what was necessary. For example the IO class was changed to use string filenames instead of const char. The main programs are MD_many_temperatures.cpp and MD_one_temperature.cpp. 

MD_plot.py is used to create the plots used in the project. The program runs by reading directory containing temp_diff_x.txt files, start x, final x and x spacing. It then produces the 4 graphs used in the project. For example it can be run as: python MD_plot.py /Data/ 100 900 10
