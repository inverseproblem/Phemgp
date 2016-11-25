# README #

### Phemgp ###

The use of phase equilibrium calculations to compute physical properties of rocks has become commonplace in geophysical modeling. Typically, the phase equilibrium calculations are used to construct two-dimensional tables of rock properties as a function of pressure and temperature. Phemgp is a Fortran program that can be used to assemble a three-dimensional table that accounts for compositional variations from two-dimensional tables. See Zunino et al., 2011. 

A quick guide on how to generate tables to be used with Phemgp is available at http://www.perplex.ethz.ch/phemgp/phemgp.html .
Two tables representing harzburgitic and basaltic end-members are provided here as an example. 


Zunino, A., J. A. D. Connolly, and A. Khan (2011), Precalculated phase equilibrium models for geophysical properties of the crust and mantle as a function of composition, Geochem. Geophys. Geosyst., 12, Q04001, doi:10.1029/2010GC003304


## Installation: ##

Compile Phemgp using a Fortran compiler. For instance, with the GNU Fortran compiler, do:
```
gfortran phemgp_0.5.f90 -o phemgp.x
```

On Linux or MacOS run it with:
```
./phemgp.x
```


## IMPORTANT: ##

In order to use Phemgp with the tables generated from Perple_X 
output (as described in http://www.perplex.ethz.ch/phemgp/phemgp.html) 
the user must add a line containing the basalt fraction in the 
header of the table. 

The line specifying basalt fraction (as a floating point number) has 
to be placed after the row containing the title of the table (e.g., 
"hzmb_200_1.phm"), as shown by the following example (basalt 
fraction=0.2):


Header of the table obtained from Perple_X: 

```
|6.6.6
hzmb_200_1.phm    
           2
T(K)    
   500.002000000000     
   40.8162448980000     
          50
P(bar)  
   1.13999900000000     
   2857.11673473755    
```



modified header:

```
|6.6.6
hzmb_200_1.phm  
0.2d0              !! <<=== here
           2
T(K)    
   500.002000000000     
   40.8162448980000     
          50
P(bar)  
   1.13999900000000     
   2857.11673473755   

```
