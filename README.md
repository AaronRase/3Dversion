# 3D version

### White Noise

Do the following in the input file: 

- In case you want to use white noise, make sure that in the input file you put *Pbool = true*.
- If set to true, also add *Pconst* and set it to the value which represent the square root of the variance of the initial scalar fluctuations.

### Calculating Xi

Do the following in the input file:

- Set *tOutputXiFreq* to the time step at which you want to calculate Xi.
- Define the methods you want to use:
  - *ximethod = 1* calculates Xi using the energy density of the walls.
  - *ximethod = 2* calculates Xi using the area trick.
  - *ximethod = 3* calculates Xi using the line of the sight.
  - Any combination of these numbers (from lowest to highest) will calculate Xi according to the associated methods, e.g. *ximethod = 12*, *ximethod = 123*, *ximethod = 23*, etc.
