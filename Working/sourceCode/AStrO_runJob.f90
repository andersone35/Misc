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
    
    character(len=256) :: fileLine, commandString
    character(len=128) :: jobScriptName, inputFileName, outputFileName, fullFileName
    character(len=16) :: fields(10)
    real*8 :: time
    real*8 :: scaleFactor
    integer :: numFields, writeGrad, setDispMode, writeModes
    integer :: readInt(3), errFlag
    integer :: i1, i2, i3, i4, i5, i6, i7, iosVal, bw3, nnd2, readModel, nSInd, eSInd
    
    call get_command_argument(1,jobScriptName)
    
    open(unit=8, file=jobScriptName, blank='NULL', action='read')
    
	abortJob = 0
    lfUnit = 30
    open(unit=lfUnit, file='jobLogFile.txt', action='write', status='replace')
    
    write(lfUnit,*) 'Beginning Job'
    write(lfUnit,*) 'Script name: ', jobScriptName
	write(*,*) 'Beginning Job'
    write(*,*) 'Script name: ', jobScriptName
    
    readModel = 0
    
    fileLine(1:15) = '               '
    
    read(8,'(A)',iostat=iosVal) fileLine(16:256)
    
    do while(iosVal .eq. 0 .and. abortJob .eq. 0)
        !!write(*,*) 'Loop start fileLine : ', fileLine
        i1 = index(fileLine,'*')
        if(i1 .gt. 0) then
            if(fileLine(i1:i1+4) .eq. '*read') then
			    commandString = fileLine
				write(*,*) 'running ', commandString
                call processReadCommand(8,fileLine,iosVal)
                i1 = index(fileLine,'*')
				!! ------
				i2 = index(commandString,'*readModelInput')
				if(i2 .gt. 0) then
				    readModel = 1
					!outputFileName = 'modelInputEcho.txt'
					!call writeModelData(outputFileName)
				endif
				! i2 = index(commandString,'*readDesignVarInput')
				! if(i2 .gt. 0) then
					! outputFileName = 'dVarInputEcho.txt'
					! call writeDesignVarData(outputFileName)
				! endif
				! i2 = index(commandString,'*readObjectiveInput')
				! if(i2 .gt. 0) then
					! outputFileName = 'objectiveEcho.txt'
					! call writeObjectiveData(outputFileName)
				! endif
				!! ------
            elseif(fileLine(i1:i1+5) .eq. '*solve') then
                if(readModel .eq. 0) then
                    write(lfUnit,*) 'Error: *solve command called before model input file read in.'
                    read(8,'(A)',iostat=iosVal) fileLine(16:256)
					abortJob = 1
                    goto 114
                endif
				write(lfUnit,*) 'reading solve command'
				write(*,*) 'reading solve command'
				solverMeth = 'direct'
                solverBlockDim = numNodes
				if(solverBlockDim .lt. 4) then
				    solverBlockdim = 4
				endif
                nnd2 = numNodes*numNodes
                solverMaxBW = 1
                bw3 = 1
                do while(bw3 .lt. nnd2)
                    solverMaxBW = solverMaxBW + 1
                    bw3 = solverMaxBW*solverMaxBW*solverMaxBW
                enddo
				if(solverMaxBW .lt. 6) then
				    solverMaxBW = 6
				endif
                solveThermal = 0
                solveElastic = 1
                nLGeom = 0
				loadTime = 1d0
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
						elseif(fileLine(i2-12:i2) .eq. 'solverMethod:') then
						    read(fileLine(i2+1:i2+64),*) solverMeth
                        elseif(fileLine(i2-14:i2) .eq. 'solverBlockDim:') then
						    read(fileLine(i2+1:i2+64),*,err=184) solverBlockDim
							goto 186
184                         write(lfUnit,*) 'Error: invalid input for solverBlockDim in line:'
                            write(lfUnit,*) fileLine
186							i1 = i1
                        elseif(fileLine(i2-18:i2) .eq. 'solverMaxBandwidth:') then
						    read(fileLine(i2+1:i2+64),*,err=190) solverMaxBW
							goto 192
190                         write(lfUnit,*) 'Error: invalid input for solverBlockDim in line:'
                            write(lfUnit,*) fileLine
192							i1 = i1
                        endif
                    endif
                    read(8,'(A)',iostat=iosVal) fileLine(16:256)
                    i1 = index(fileLine,'*')
                enddo
				if(solverMeth .eq. 'direct') then
				     solverBlockDim = numNodes
					 solverMaxBW = 6*numNodes
				endif
				write(lfUnit,*) 'calling solve'
				close(lfUnit)
				open(unit=lfUnit, file='jobLogFile.txt', action='write', access='append')
                call solve()
114             numNodes = numNodes
            elseif(fileLine(i1:i1+13) .eq. '*modalAnalysis') then
			    if(.not. allocated(currentRank)) then
				    call analysisPrep()
				endif
                modalType = 0 !! 0 = buckling, 1 = vibration
                numEigModes = 10
                read(8,'(A)',iostat=iosVal) fileLine(16:256)
                i1 = index(fileLine,'*')
                do while(i1 .eq. 0 .and. iosVal .eq. 0)
                    i2 = index(fileLine,':')
                    if(i2 .gt. 0) then
                        if(fileLine(i2-8:i2) .eq. 'numModes:') then
                            read(fileLine(i2+1:i2+64),*) numEigModes
                            read(8,'(A)',iostat=iosVal) fileLine(16:256)
                            i1 = index(fileLine,'*')
                        elseif(fileLine(i2-4:i2) .eq. 'type:') then
                            i3 = index(fileLine,'vibration')
                            if(i3 .gt. 0) then
                                modalType = 1
                            endif
                            read(8,'(A)',iostat=iosVal) fileLine(16:256)
                            i1 = index(fileLine,'*')
                        else
                            read(8,'(A)',iostat=iosVal) fileLine(16:256)
                            i1 = index(fileLine,'*')
                        endif
                    else
                        read(8,'(A)',iostat=iosVal) fileLine(16:256)
                        i1 = index(fileLine,'*')
                    endif
                enddo
                allocate(eigenModes(elMatDim,numEigModes))
                allocate(eigenVals(numEigModes))
                allocate(eigenFactors(numEigModes))
                allocate(diagMassMat(elMatDim))
				eigenModes(:,:) = r_0
                eigenVals(:) = r_0
                eigenFactors(:) = r_0
                diagMassMat(:) = r_0
                i3 = dynamic
                dynamic = 0
                call getElasticSolnLoad(1)
                dynamic = i3
                if(modalType .eq. 1) then
                    call getDiagMassMat()
                    do i3 = 1, elMatDim
                        diagMassMat(i3) = r_1/diagMassMat(i3)
                    enddo
                else
                    diagMassMat(:) = r_1
                endif
				if(solverMeth .eq. 'iterative') then
					call getSparseEigenModes(eigenModes,eigenVals,numEigModes,elMatDim,elasticMat,elMatCols,elMatRange,elMatSize, &
						diagMassMat,elMPCMat,elMPCMatCols,elMPCMatRange,elMPCDim,elMPCSize)
				else
				    call getEigModesLDL(eigenModes,eigenVals,numEigModes,elMatDim,elMatLT,elMatLTRange,elMatLTSize,diagMassMat)
				endif
                if(modalType .eq. 0) then
                    call getBucklingLoadFactors(eigenFactors,numEigModes)
                else
                    do i3 = 1, numEigModes
                        eigenFactors(i3) = r_p5*sqrt(eigenVals(i3))/r_pi
                    enddo
                endif
            elseif(fileLine(i1:i1+13) .eq. '*setDispToMode') then
                scaleFactor = 0.01d0
                setDispMode = 1
                read(8,'(A)',iostat=iosVal) fileLine(16:256)
                i1 = index(fileLine,'*')
                do while(i1 .eq. 0 .and. iosVal .eq. 0)
                    i2 = index(fileLine,':')
                    if(i2 .gt. 0) then
                        if(fileLine(i2-4:i2) .eq. 'mode:') then
                            read(fileLine(i2+1:i2+64),*) setDispMode
                            read(8,'(A)',iostat=iosVal) fileLine(16:256)
                            i1 = index(fileLine,'*')
                        elseif(fileLine(i2-11:i2) .eq. 'scaleFactor:') then
                            read(fileLine(i2+1:i2+64),*) scaleFactor
                            read(8,'(A)',iostat=iosVal) fileLine(16:256)
                            i1 = index(fileLine,'*')
                        else
                            read(8,'(A)',iostat=iosVal) fileLine(16:256)
                            i1 = index(fileLine,'*')
                        endif
                    else
                        read(8,'(A)',iostat=iosVal) fileLine(16:256)
                        i1 = index(fileLine,'*')
                    endif
                enddo
                call setSolnToMode(nodeDisp,elMatDim,setDispMode,scaleFactor)
            elseif(fileLine(i1:i1+13) .eq. '*calcObjective') then
                call getCompleteObjective()
                read(8,'(A)',iostat=iosVal) fileLine(16:256)
                i1 = index(fileLine,'*')
            elseif(fileLine(i1:i1+15) .eq. '*calcObjGradient') then
                call getTotaldLdD()
                read(8,'(A)',iostat=iosVal) fileLine(16:256)
                i1 = index(fileLine,'*')
            elseif(fileLine(i1:i1+16) .eq. '*writeNodeResults' .or. fileLine(i1:i1+19) .eq. '*writeElementResults') then
			    write(*,*) 'writing node results'
                call processWriteNodeEl(8,fileLine,iosVal)
                i1 = index(fileLine,'*')
            elseif(fileLine(i1:i1+17) .eq. '*writeModalResults') then
                writeModes = 1
                read(8,'(A)',iostat=iosVal) fileLine(16:256)
                i1 = index(fileLine,'*')
                do while(i1 .eq. 0 .and. iosVal .eq. 0)
                    i2 = index(fileLine,':')
                    if(i2 .gt. 0) then
                        if(fileLine(i2-8:i2) .eq. 'fileName:') then
                            read(fileLine(i2+1:i2+128),*) outputFileName
							outputFileName = adjustl(outputFileName)
							outputFileName = trim(outputFileName)
                            read(8,'(A)',iostat=iosVal) fileLine(16:256)
                            i1 = index(fileLine,'*')
                        elseif(fileLine(i2-10:i2) .eq. 'writeModes:') then
                            i2 = index(fileLine,'no')
                            if(i2 .gt. 1) then
                                writeModes = 0
                            endif
                            read(8,'(A)',iostat=iosVal) fileLine(16:256)
                            i1 = index(fileLine,'*')
                        else
                            read(8,'(A)',iostat=iosVal) fileLine(16:256)
                            i1 = index(fileLine,'*')
                        endif
                    else
                        read(8,'(A)',iostat=iosVal) fileLine(16:256)
                        i1 = index(fileLine,'*')
                    endif
                enddo
                call writeModalResults(outputFileName,writeModes)
            elseif(fileLine(i1:i1+20) .eq. '*writeNodeCoordinates' .or. fileLine(i1:i1+22) .eq. '*writeElementProperties') then
                call processWriteNdElProps(8,fileLine,iosVal)
                i1 = index(fileLine,'*')
            elseif(fileLine(i1:i1+14) .eq. '*writeObjective') then
                writeGrad = 1
                read(8,'(A)',iostat=iosVal) fileLine(16:256)
                i1 = index(fileLine,'*')
                do while(i1 .eq. 0 .and. iosVal .eq. 0)
                    i2 = index(fileLine,':')
                    if(i2 .gt. 0) then
                        if(fileLine(i2-8:i2) .eq. 'fileName:') then
                            read(fileLine(i2+1:i2+128),*) outputFileName
							outputFileName = adjustl(outputFileName)
							outputFileName = trim(outputFileName)
                            read(8,'(A)',iostat=iosVal) fileLine(16:256)
                            i1 = index(fileLine,'*')
                        elseif(fileLine(i2-7:i2) .eq. 'include:') then
                            read(8,'(A)',iostat=iosVal) fileLine(16:256)
                            i2 = index(fileLine,':')
                            numFields = 0
                            do while(i2 .eq. 0 .and. iosVal .eq. 0)
                                i3 = index(fileLine,'-')
                                if(i3 .gt. 0) then
                                    numFields = numFields + 1
                                    read(fileLine(i3+1:i3+15),*) fields(numFields)
                                endif
                                read(8,'(A)',iostat=iosVal) fileLine(16:256)
                                i2 = index(fileLine,':')
                            enddo
                            i1 = index(fileLine,'*')
                        elseif(fileLine(i2-13:i2) .eq. 'writeGradient:') then
                            i3 = index(fileLine,'no')
                            if(i3 .gt. 0) then
                                writeGrad = 0
                            endif
                            read(8,'(A)',iostat=iosVal) fileLine(16:256)
                            i1 = index(fileLine,'*')
                        else
                            read(8,'(A)',iostat=iosVal) fileLine(16:256)
                            i1 = index(fileLine,'*')
                        endif
                    else
                        read(8,'(A)',iostat=iosVal) fileLine(16:256)
                        i1 = index(fileLine,'*')
                    endif
                enddo
                call writeObjective(outputFileName,fields,numFields,writeGrad)
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