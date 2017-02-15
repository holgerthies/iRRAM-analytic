# Evaluation method
All data is output from the ode2 script contained in test/ode2.cc
The script reads the following from the standard input.
- Dimension of the system (1,2 or 3)
- System ID (see below)
- Solver method (0 = Algorithm A or 1 = Algorithm B)
- Step size method (1 = CS-VO, 2 = CO-VS (order 6), 3 = CO-VS (order 50), 4 = VS-VO)
- prune function before running algorithm (0 = do not use prune operator, 1 = use prune operator)
- value t up to which the trajectory should be followed
- output precision in number of decimals (0 = do not specify)
- number of iterations (only if output precision is not specified) 

If the precision is given it outputs for each iteration 
- cpu time for the iteration
- number of steps 
- highest order (number of taylor coefficients of the solution function read) during the iteration
- start precision for the iteration
- mantissa and exponent for the error of the value y(t) of the solution function at the last step
- exponent if error is written in the form 2^n
Otherwise at each step the time t, and all values of y(t) are output with the given precision.

The data in the file data.txt is generated as follows:
the program ode2 is run with parameter --prec_inc=1 and as value for t 1.0 is chosen. 
The other parameters are as given in each entry of the file.

The files outstiff-{a}{b} are generated as follows:
The program ode2 is run with dimension 2, system id 9, solver method a, step size method b and prune=0.
Output precision is set to 22 but lower for step size method 1 (as it does not matter for the step size).

# Test functions
## dimension 1
- 0 : y'(t) = -0.5y^3; y(0) = 1
- 1 : y'(t) = cos(y); y(0) = 0
- 2 : y'(t) = 1/(y+10); y(0) = 0

## dimension 2
- 0 : y1'(t) = 2(y1-y1y2); y2'(t) = y1y2-y2; y(0) = (1,3)
- 1 : y1'(t) = y1cos(y2); y2'(t)=1; y(0) = (1,0)
- 2 : y1'(t) = (y1-y2+4)/(y1+y2+4); y2'(t)=1; y(0) = (0,0)
- 3 : y1'(t) = 1/(y1+y2+10); y2'(t)=1; y(0) = (0,0)
- 4 : y1'(t) = y2; y2'(t) = 10((1-y1^2)y2-y1); y(0) = (1.2, -0.6)
- 5 : y1'(t) = y2; y2'(t) = 100((1-y1^2)y2-y1); y(0) = (1.2, -0.6)
- 6 : y1'(t) = y2; y2'(t) = 1000((1-y1^2)y2-y1); y(0) = (1.2, -0.6)
- 7 : y1'(t) = 1000y2; y2'(t) = 0.1((1-y1^2)1000y2-y1); y(0) = (1.2, -0.0006)
- 8 : y1'(t) = 1000y2; y2'(t) = ((1-y1^2)1000y2-y1); y(0) = (1.2, -0.0006)
- 9 : y1'(t) = 1000y2; y2'(t) = 0.01((1-y1^2)1000y2-y1); y(0) = (1.2, -0.0006)

## dimension 3
- 0 : y1'(t) = y1y2+1; y2'(t) = y1y3+y2; y3'(t) = y1+y2+y3; y(0) = (0,0,0)
