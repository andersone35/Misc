program AStrO_runJob
    use AStrO_constantVals
	use AStrO_globalData
	use AStrO_input
	use AStrO_output
	use AStrO_r_overloadFunctions
	use AStrO_c_overloadFunctions
	use AStrO_solvers
	use AStrO_r_designPropertyFunctions
	use AStrO_c_designPropertyFunctions
	use AStrO_r_elementEqns
	use AStrO_c_elementEqns
	use AStrO_bookKeeping
	use AStrO_objective
	use AStrO_commandFunctions
    implicit none
	
	integer, allocatable :: allNdSet(:), allElSet(:), writeTSteps(:)
	
	character(len=256) :: fileLine
	character(len=128) :: jobScriptName, inputFileName, outputFileName, fullFileName
	character(len=64) :: writeFields(10)
	character(len=16) :: extension, stepStr
	real*8 :: time
	integer :: readInt(3), errFlag
	integer :: i1, i2, i3, i4, i5, i6, i7, iosVal, bw3, nnd2, readModel, nSInd, eSInd, numFields
	
	call get_command_argument(1,jobScriptName)
	
	open(unit=8, file=jobScriptName, blank='NULL', action='read')
	
	lfUnit = 30
	open(unit=lfUnit, file='jobLogFile.txt', action='write', status='replace')
	
	write(lfUnit,*) 'Beginning Job'
	write(lfUnit,*) 'Script name: ', jobScriptName
	
	readModel = 0
	
	fileLine(1:15) = '               '
	
	read(8,'(A)',iostat=iosVal) fileLine(16:256)
	
	do while(iosVal .eq. 0)
	    !!write(*,*) 'Loop start fileLine : ', fileLine
	    i1 = index(fileLine,'*')
		if(i1 .gt. 0) then
			if(fileLine(i1:i1+14) .eq. '*readModelInput') then
				read(8,'(A)',iostat=iosVal) fileLine(16:256)
				i1 = index(fileLine,'*')
				do while(i1 .eq. 0 .and. iosVal .eq. 0)
					i2 = index(fileLine,':')
					if(i2 .gt. 0) then
						if(fileLine(i2-8:i2) .eq. 'fileName:') then
							read(fileLine(i2+1:i2+128),*) inputFileName
							write(lfUnit,*) 'reading model input, file: ', inputFileName
							call readModelInput(inputFileName)
							readModel = 1
						endif
					endif
					read(8,'(A)',iostat=iosVal) fileLine(16:256)
					i1 = index(fileLine,'*')
				enddo
				outputFileName = 'modelDataEcho.txt'
				call writeModelData(outputFileName)
			elseif(fileLine(i1:i1+18) .eq. '*readDesignVarInput') then
				read(8,'(A)',iostat=iosVal) fileLine(16:256)
				i1 = index(fileLine,'*')
				do while(i1 .eq. 0 .and. iosVal .eq. 0)
					i2 = index(fileLine,':')
					if(i2 .gt. 0) then
						if(fileLine(i2-8:i2) .eq. 'fileName:') then
							read(fileLine(i2+1:i2+128),*) inputFileName
							write(lfUnit,*) 'reading design variable input, file: ', inputFileName
							call readDesignVarInput(inputFileName)
						endif
					endif
					read(8,'(A)',iostat=iosVal) fileLine(16:256)
					i1 = index(fileLine,'*')
				enddo
				outputFileName = 'dVarDataEcho.txt'
				call writeDesignVarData(outputFileName)
			elseif(fileLine(i1:i1+18) .eq. '*readObjectiveInput') then
			    read(8,'(A)',iostat=iosVal) fileLine(16:256)
				i1 = index(fileLine,'*')
				do while(i1 .eq. 0 .and. iosVal .eq. 0)
					i2 = index(fileLine,':')
					if(i2 .gt. 0) then
						if(fileLine(i2-8:i2) .eq. 'fileName:') then
							read(fileLine(i2+1:i2+128),*) inputFileName
							write(lfUnit,*) 'reading objective input, file: ', inputFileName
							call readObjectiveInput(inputFileName)
						endif
					endif
					read(8,'(A)',iostat=iosVal) fileLine(16:256)
					i1 = index(fileLine,'*')
				enddo
				call getdToElComp()
				outputFileName = 'objectiveEcho.txt'
				call writeObjectiveData(outputFileName)
			elseif(fileLine(i1:i1+5) .eq. '*solve') then
			    if(readModel .eq. 0) then
				    write(lfUnit,*) 'Error: *solve command called before model input file read in.'
					read(8,'(A)',iostat=iosVal) fileLine(16:256)
					goto 114
				endif
			    solverBlockDim = numNodes
				nnd2 = 27*numNodes*numNodes
				solverMaxBW = 1
				bw3 = 1
				do while(bw3 .lt. nnd2)
				    solverMaxBW = solverMaxBW + 1
					bw3 = solverMaxBW*solverMaxBW*solverMaxBW
				enddo
				solveThermal = 0
				solveElastic = 1
				nLGeom = 0
				dynamic = 0
				writeSolnHist = 0
				numTSteps = 0
				nMBeta = 0.25d0
				nMGamma = 0.5d0
				delT = 0d0
				simPeriod = 0d0
				rayCoefK = 0d0
				rayCoefM = 0d0
				loadTime = 0d0
				read(8,'(A)',iostat=iosVal) fileLine(16:256)
				i1 = index(fileLine,'*')
				do while(i1 .eq. 0 .and. iosVal .eq. 0)
				    i2 = index(fileLine,':')
					if(i2 .gt. 0) then
					    if(fileLine(i2-7:i2) .eq. 'elastic:') then
						    i3 = index(fileLine,'no')
							if(i3 .gt. 0) then
							    solveElastic = 0
							endif
						elseif(fileLine(i2-7:i2) .eq. 'thermal:') then
						    i3 = index(fileLine,'no')
							if(i3 .gt. 0) then
							    solveThermal = 0
							endif
						elseif(fileLine(i2-13:i2) .eq. 'nonLinearGeom:') then
						    i3 = index(fileLine,'yes')
							if(i3 .gt. 0) then
							    nLGeom = 1
							endif
						elseif(fileLine(i2-14:i2) .eq. 'staticLoadTime:') then
						    read(fileLine(i2+1:i2+64),*,err=137) loadTime
							goto 140
137                         write(lfUnit,*) 'Error: invalid input for staticLoadTime in line:'
                            write(lfUnit,*) fileLine
140                         i1 = i1
                        elseif(fileLine(i2-7:i2) .eq. 'dynamic:') then
						    i3 = index(fileLine,'yes')
							if(i3 .gt. 0) then
							    dynamic = 1
							endif
						elseif(fileLine(i2-8:i2) .eq. 'timeStep:') then
						    read(fileLine(i2+1:i2+64),*,err=150) delT
							goto 152
150                         write(lfUnit,*) 'Error: invalid input for timeStep in line:'
                            write(lfUnit,*) fileLine
152                         i1 = i1
                        elseif(fileLine(i2-11:i2) .eq. 'newmarkBeta:') then
						    read(fileLine(i2+1:i2+64),*,err=156) nMBeta
							goto 158
156                         write(lfUnit,*) 'Error: invalid input for newmarkBeta in line:'
                            write(lfUnit,*) fileLine
158                         i1 = i1
                        elseif(fileLine(i2-12:i2) .eq. 'newmarkGamma:') then
						    read(fileLine(i2+1:i2+64),*,err=162) nMGamma
							goto 164
162                         write(lfUnit,*) 'Error: invalid input for newmarkGamma in line:'
                            write(lfUnit,*) fileLine
164                         i1 = i1
                        elseif(fileLine(i2-9:i2) .eq. 'simPeriod:') then
						    read(fileLine(i2+1:i2+64),*,err=168) simPeriod
							goto 170
168                         write(lfUnit,*) 'Error: invalid input for timeStep in line:'
                            write(lfUnit,*) fileLine
170                         i1 = i1
                        elseif(fileLine(i2-12:i2) .eq. 'saveSolnHist:') then
                            i3 = index(fileLine,'yes')
							if(i3 .gt. 0) then
							    writeSolnHist = 1
							endif		
						endif
					endif
					read(8,'(A)',iostat=iosVal) fileLine(16:256)
					i1 = index(fileLine,'*')
				enddo
				call solve()
114             numNodes = numNodes
            elseif(fileLine(i1:16) .eq. '*writeNodeResults') then
			    read(8,'(A)',iostat=iosVal) fileLine(16:256)
				i1 = index(fileLine,'*')
				nSInd = 0
				do while(i1 .eq. 0 .and. iosVal .eq. 0)
					i2 = index(fileLine,':')
					if(i2 .gt. 0) then
						if(fileLine(i2-8:i2) .eq. 'fileName:') then
							read(fileLine(i2+1:i2+128),*) outputFileName
							read(8,'(A)',iostat=iosVal) fileLine(16:256)
							i1 = index(fileLine,'*')
						elseif(fileLine(i2-7:i2) .eq. 'nodeSet:') then
						    do i3 = 1, numNdSets
							    i4 = index(fileLine,ndSetName(i3))
							    if(i4 .gt. 0) then
								    nSInd = i3
								endif
							enddo
							if(nSInd .eq. 0) then
							    write(lfUnit,*) 'Warning: no node set named ', fileLine(i2+1:i2+64), 'was found'
								write(lfUnit,*) 'Writing results for all nodes'
							endif
							read(8,'(A)',iostat=iosVal) fileLine(16:256)
							i1 = index(fileLine,'*')
						elseif(fileLine(i2-6:i2) .eq. 'fields:') then
						    read(8,'(A)',iostat=iosVal) fileLine(16:256)
							i2 = index(fileLine,':')
							numFields = 0
							do while(i2 .eq. 0 .and. iosVal .eq. 0)
							    i3 = index(fileLine,'-')
								if(i3 .gt. 0) then
								    numFields = numFields + 1
								    read(fileLine(i3+1:i3+64),*) writeFields(numFields)
									if(writeFields(numFields) .eq. 'reactionForce') then
									    elasticLoad(:) = r_0
										if(intVecSize .gt. 0) then
										    intElasticLoad(:) = r_0
										endif
									    call getElasticSolnLoad(0)
									endif
									if(writeFields(numFields) .eq. 'reactionHeatGen') then
									    thermalLoad(:) = r_0
										call getThermalSolnLoad(0)
									endif
								endif
								read(8,'(A)',iostat=iosVal) fileLine(16:256)
								i2 = index(fileLine,':')
							enddo
							i1 = index(fileLine,'*')
						elseif(fileLine(i2-9:i2) .eq. 'timeSteps:') then
						    if(numTSteps .le. 0) then
							    write(lfUnit,*) 'Error: A dynamic analysis must be run prior to writing time step dependent results.'
								write(lfUnit,*) 'Aborting writeNodeResults'
								goto 251
							endif
						    if(.not. allocated(writeTSteps)) then
							    allocate(writeTSteps(numTSteps))
							endif
							writeTSteps(:) = 0
						    i3 = index(fileLine,'all')
							if(i3 .gt. 0) then
							    writeTSteps(:) = 1
							else
							    read(8,'(A)',iostat=iosVal) fileLine(16:256)
								i2 = index(fileLine,':')
								do while(i2 .eq. 0 .and. iosVal .eq. 0)
									i3 = index(fileLine, '-')
									if(i3 .gt. 0) then
										i4 = index(fileLine,'[')
										if(i4 .gt. 0) then
											i5 = index(fileLine,']')
											readInt(3) = 1
											read(fileLine(i4+1:i5-1),*,end=254) readInt(1:3)
254                                         do i6 = readInt(1), readInt(2), readInt(3)
                                                if(i6 .le. numTSteps .and. i6 .gt. 0) then
												    writeTSteps(i6) = 1
												else
												    write(lfUnit,*) 'Error: a time step specified for writeNodeResults is out of'
													write(lfUnit,*) 'the simulated range.  Aborting writeNodeResults.'
													goto 251
												endif
                                            enddo
										else
										    read(fileLine(i3+1:i3+64),*) i5
											if(i5 .le. numTSteps .and. i5 .gt. 0) then
											    writeTSteps(i5) = 1
											else
											    write(lfUnit,*) 'Error: a time step specified for writeNodeResults is out of'
												write(lfUnit,*) 'the simulated range.  Aborting writeNodeResults.'
											    goto 251
											endif
										endif
									endif
									read(8,'(A)',iostat=iosVal) fileLine(16:256)
								    i2 = index(fileLine,':')
								enddo
							endif
						else
						    read(8,'(A)',iostat=iosVal) fileLine(16:256)
							i1 = index(fileLine,'*')
						endif
					endif
				enddo
				if(.not. allocated(writeTSteps)) then
				    if(nSInd .eq. 0) then
					    if(.not. allocated(allNdSet)) then
						    allocate(allNdSet(numNodes))
						endif
						do i3 = 1, numNodes
							allNdSet(i3) = i3
						enddo
						time = numTSteps*delT
						call writeNodeResults(outputFileName,writeFields,numFields,allNdSet,numNodes,time,numTSteps)
					else
					    i3 = ndSetRange(nSInd-1) + 1
						i4 = ndSetRange(nSInd)
						i5 = i4 - i3 + 1
						time = numTSteps*delT
						call writeNodeResults(outputFileName,writeFields,numFields,nodeSets(i3:i4),i5,time,numTSteps)
					endif
				else
				    do i6 = 1, numTSteps
					    if(writeTSteps(i6) .eq. 1) then
						    call readBinarySolution(i6,errFlag)
							if(errFlag .eq. 0) then
							    nodeTemp = prevTemp(:)
								nodeTdot = prevTdot(:)
								nodeDisp = prevDisp(:)
								if(intVecSize .ge. 1) then
									internalDisp(:) = prevIntDisp(:)
								endif
								nodeVel(:) = prevVel(:)
								nodeAcc(:) = prevAcc(:)
								write(stepStr,'(I0)') i6
								i7 = index(outputFileName,'.')
								if(i7 .gt. 0) then
									read(outputFileName(i7:i7+15),*) extension
									read(outputFileName(1:i7-1),*) fullFileName
									outputFileName = fullFileName
								else
									extension = ''
								endif
								fullFileName = outputFileName // '_step' // stepStr // extension
								if(nSInd .eq. 0) then
									if(.not. allocated(allNdSet)) then
										allocate(allNdSet(numNodes))
									endif
									do i3 = 1, numNodes
										allNdSet(i3) = i3
									enddo
									time = numTSteps*delT
									call writeNodeResults(fullFileName,writeFields,numFields,allNdSet,numNodes,time,numTSteps)
								else
									i3 = ndSetRange(nSInd-1) + 1
									i4 = ndSetRange(nSInd)
									i5 = i4 - i3 + 1
									time = numTSteps*delT
									call writeNodeResults(fullFileName,writeFields,numFields,nodeSets(i3:i4),i5,time,numTSteps)
								endif
							else
							    write(lfUnit,*) 'Error: time step ', i6, ' was never recorded in the solution history.'
								write(lfUnit,*) 'Be sure to set the saveSolnHist option to yes under the *solve command'
								write(lfUnit,*) 'if the objective gradient or time-specific intermediate results are desired'
								write(lfUnit,*) 'for a dynamic analysis.'
							endif
						endif
					enddo
					deallocate(writeTSteps)
				endif
				if(allocated(allNdSet)) then
				    deallocate(allNdSet)
				endif
251				i1 = index(fileLine,'*')
			else
				read(8,'(A)',iostat=iosVal) fileLine(16:256)
				i1 = index(fileLine,'*')
			endif
		else
		    read(8,'(A)',iostat=iosVal) fileLine(16:256)
			i1 = index(fileLine,'*')
		endif
	enddo
	
	call deallocateAll()
	
	close(lfUnit)
	close(8)

end program