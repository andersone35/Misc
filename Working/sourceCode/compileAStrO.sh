#!/bin/bash
gfortran -static -c sourceCode/AStrO_constantVals.f90
gfortran -static -c sourceCode/AStrO_globalData.f90
gfortran -static -c sourceCode/AStrO_input.f90
gfortran -static -c sourceCode/AStrO_r_overloadFunctions.f90
gfortran -static -c sourceCode/AStrO_c_overloadFunctions.f90
gfortran -static -c sourceCode/AStrO_solvers.f90
gfortran -static -c sourceCode/AStrO_r_designPropertyFunctions.f90
gfortran -static -c sourceCode/AStrO_c_designPropertyFunctions.f90
gfortran -static -c sourceCode/AStrO_r_elementEqns.f90
gfortran -static -c sourceCode/AStrO_c_elementEqns.f90
gfortran -static -c sourceCode/AStrO_output.f90
gfortran -static -c sourceCode/AStrO_bookKeeping.f90
gfortran -static -c sourceCode/AStrO_objective.f90
gfortran -static -c sourceCode/AStrO_commandFunctions.f90
gfortran -static sourceCode/AStrO_runJob.f90 AStrO_constantVals.o AStrO_globalData.o AStrO_input.o AStrO_r_overloadFunctions.o AStrO_c_overloadFunctions.o AStrO_solvers.o AStrO_r_designPropertyFunctions.o AStrO_r_elementEqns.o AStrO_c_designPropertyFunctions.o AStrO_c_elementEqns.o AStrO_output.o AStrO_bookKeeping.o AStrO_objective.o AStrO_commandFunctions.o -o bin/AStrO_runJob.exe

./sourceCode/cleanUp.sh
## End of script
