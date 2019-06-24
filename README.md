1)AnalysisTempPixel:                                                                                                                                                                                                                       
                                                                                                                                                                                                                                           
Analyse data acquired using pixel system to study LO:                                                                                                                                                                                      
-make to compile                                                                                                                                                                                                                           
-to launch the analysis:                                                                                                                                                                                                                    
./AnalysisTemp.exe $Path/To/Data/Directory$ OvValue                                                                                                                                                                                        
ex: ./AnalysisTemp.exe ../TestStabilityNa22PedAllChannels_01_04_2019 7                                                                                                                                                                     
Important:                                                                                                                                                                                                                                 
1) do not put the "/" after the directory                                                                                                                                                                                                 
2) the overvoltage value must be specified in the name of the file                                                                                                                                                                         
                                                                                                                                                                                                                                           
It produce a directory called "Plot" into the directory of the data, inside there is the directory of the analysis on relative LO of the crystal called "EnergyTempCB" containing:                                                         
-The spectra                                                                                                                                                                                                                               
-The plot fitted with the calibration function                                                                                                                                                                                             
-The parameter alfa and beta as function of the temperature                                                                                                                                                                                
-The peak position in temperature                                                                                                                                                                                                          
-The energy resolution in temperature on the 511KeV peak                                                                                                                                                                                   
                                                                                                                                                                                                                                           
It produce also a root file inside the directory RootFileGraphPixel to do additional analysis on all the data acquired.                                                                                                                    
                                                                                                                                                                                                                                           
                                                                                                                                                                                                                                           
2)TimeReolutionPixel:                                                                                                                                                                                                                      
Analyse data acquired using pixel system to study the time resolution:                                                                                                                                                                     
-make to compile                                                                                                                                                                                                                           
-to launch the analysis:                                                                                                                                                                                                                    
./AnalysisTemp.exe $Path/To/Data/Directory$ OvValue                                                                                                                                                                                        
ex: ./AnalysisTemp.exe ../TestStabilityNa22PedAllChannels_01_04_2019 7                                                                                                                                                                     
Important:                                                                                                                                                                                                                                 
1) do not put the "/" after the directory                                                                                                                                                                                                 
2) the overvoltage value must be specified in the name of the file                                                                                                                                                                         
                                                                                                                                                                                                                                           
It produce a directory called "Plot" into the directory of the data, inside there is the directory of the analysis on time resolution of the crystal called "TimeRes" containing:                                                          
-The spectra in coincidence                                                                                                                                                                                                                
-The time resolution with events selected at the PE peaks                                                                                                                                                                                  
-The time resolution as function of temperature (with projection inside the directory Tdiff                                                                                                                                                
-The analysis of time resolution vs energy -> plots of tdiff vs E, the correction and the time resolution vs energy with relative projections in the corresponding directory.                                                              
                                                                                                                                                                                                                                           
It produce also a root file inside the directory RootFileGraphPixel to do additional analysis on all the data acquired.                                                                                                                    
                                                                                                                                                                                                                                           
3)RootFileGraphPixel:                                                                                                                                                                                                                      
Once the 2 program before have performed the analysis of different data acquisition, and the root files have been created (both of the program must have analysed the data and must have created the root file for a single data acquisition ), it performs an analysis of all the data acquired.  
-make to compile 
-to launch the program:
./GlobalPlot.exe “COND” (Cond perform the merging of the data, “1” if the data has been changed from the last run of the program, “0” otherwise).


4)Bar/AnalysisTempNewBar

Analyse data acquired using bar+pixel system to study LO: 
-make new to compile                                                                                                                                                                                                                           
-to launch the analysis:                                                                                                                                                                                                                    
./AnalysisTempNewBar.exe $Path/To/Data/Directory$ OvValue                                                                                                                                                                                        
ex: ./AnalysisTempNewBar.exe ../TestStabilityBar_30_04_2019 -1                                                                                                                                                                     
Important:                                                                                                                                                                                                                                 
1) do not put the "/" after the directory                                                                                                                                                                                                 
2) the overvoltage value must be specified in the name of the file

It produce a directory called "Plot" into the directory of the data, inside there is the directory of the analysis on relative LO of the bar called "EnergyTempCB" containing the same identical analysis done for the pixels.


5)Bar/BarAnalysis

Analyse data acquired using bar+pixel system to study the time resolution:
-make bar or make to compile                                                                                                                                                                                                                           
-to launch the analysis:                                                                                                                                                                                                                    
./BarAnalysis.exe $Path/To/Data/Directory$ OvValue                                                                                                                                                                                        
ex: ./BarAnalysis.exe ../TestStabilityBar_30_04_2019 -1                                                                                                                                                                     
Important:                                                                                                                                                                                                                                 
1) do not put the "/" after the directory                                                                                                                                                                                                 
2) the overvoltage value must be specified in the name of the file

It produce a directory called "Plot" into the directory of the data, inside there is the directory of the analysis on time resolution of the bar called "TR" containing:
-The pedestal plots
-The spectra in coincidence
-The diff distribution at PE peak
-The analysis of time resolution in temperature (with relative directory containing the projections)
-The analysis of time resolution vs Energy
-The analysis of time resolution vs Effective amplitude 
-The analysis of N and C in temperature  (with relative single plot of time res vs E fitted for different temperature)

6) RootFileGraphBar:
Once the analysis of LO has been performed for different sets of data, the global analysis of all the data acquired is performed.
-make to compile 
-to launch the program:
./GlobalPlot.exe “0” (only “0” is a valid option)

   