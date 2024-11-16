# MOST
Rcode for mixed-type optimal subsampling technique

This file contains a description of the codes involved in the paper.

qqoss: CPP files needed in the R code.

Note1: By invoking the "qqoss.cpp" file, we can utilize C++ to execute LEV, IBOSS and MOST, enhancing the computational speed.

mysimu: R code for numerical experiments when the full data size N increases

Note1: Using the "mysimu" function, we can simulate three data distribution scenarios (from case1 to case3), and obtain the MSE for the parameter and the MSPE for the response, as the full data N increases, and n=4*10^3.

mysimu2: R code for numerical experiments when the subdata size n increases

Note1: Using the "mysimu2" function, we can simulate three data distribution scenarios (from case1 to case3), and obtain the MSE for the parameter and the MSPE for the response, as the subdata n increases, and N=5*10^4.
