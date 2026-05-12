##########################################

Bhuvan Ravi Kumar Nozzle project

##########################################

This file includes information on:
1) How to run the project
2) Files included with the project script
3) About the Main script
4) Outputs *(must read for professor/evaluator)
5) Learning outcomes


How to run the project:

Included in the zip folder is a source code folder which contains all the code.
   - The main file to be run is 'Bhuvan_R_main_script.py'
   - Once the python script is loaded on the screen, make sure the src folder is opened in the IDE with the input files and functions. That should be enough to run the main script.
   - If the output or outputs are not generated due to errors, try running the script a couple of times. 
   - Caution: Possibility of errors due to the creation of whitespaces exists if opened in IDE's apart from spyder (as this project was written in spyder)

Files included with the project script:
  1) 'func.py' which is a script with a function to compute the Area ratio where a normal shock stands
  2) 'Ravi Kumar_Bhuvan_input.txt" which is a txt file provided to me for the project.
  3) 'Ravi Kumar_Bhuvan_nozzle_geometry.csv' which contains the nozzle geometry provided to me for the project
  4) 'standard_atmosphere.py' which is a script provided to me that will output the ambient pressure, temperature, and density of the atmosphere as a function of altitude up to 86km accurately


About the main script:

   - The beginning of the script has a section of different functions defined to compute flow property ratios like sub-sonic and super-sonic mach no. from known area mach ratios, T/To ratios, P/Po ratios, density ratios, and thrust coefficient Cf.

   - The next section reads and loads data from the given input files into variables followed by defining physical constants.

   - In the third section, the computation required for Task 1 is performed. First, the Mach no. is calculated by determining the area ratios. Po, To, rho_o, and ho are determined for given chamber pressures. 

   - In the fourth section, computations for Task 2 are performed. The static pressure inside is computed and verified for any presence of a normal shock and based on this, other properties like temperature, mach no., local speed of sound 'a', flow velocity, static density, and enthalpy are computed. 

   - Sizing of nozzle geometry to find the nozzle area is done along with and few task 3 values are determined. This is followed by Task 2 graph plotting.

   - Next, Task 3 tabulation is made followed by Task 4 section which includes computing nozzle flow pressures corresponding to different ambient pressures and plotting them in a single graph.

   - Finally, the Task 5 section computes the Cf and Isp as a function of altitude from sea level up to 86km using the standard atmosphere function and plotting them.

Outputs: 
  1) Task 1 computation and sizing completed
  2) Task 2 plots are generated as expected
  3) Task 3 values are listed
  4) Task 4 plot generated is not the expected result. Despite my great effort of painstakingly evaluating individual code lines, logic, and values and verifying them with hand calculations, there is some discrepancy (arising due to the code or coder) in determining the area location of normal shock and computing the pressure downstream of the normal shock inside the nozzle. I require some evaluation and guidance from the professor to rectify this issue in my project.
  5) Task 5 plot generated with the expected profile 

Learning Outcome:

The greatest learning outcome I have received from working on this project is a firm grasp on the theoretical concept of flow properties inside a de Laval nozzle design.

My second learning outcome has been my confidence in working with python. Learning how to use different libraries and implementing them in a very short time has been challenging as well as rewarding. Overall, from having maybe 1% knowledge of Python to finishing the project in a cumulative effort of ~ 48 hours has been a rollercoaster ride for me with swollen eyes, dizzying headache, and severe lack of sleep, but these learning outcomes are going to benefit me in my career for a long time.

