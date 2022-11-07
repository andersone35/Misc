module AStrO_commandFunctions
    use AStrO_globalData
	use AStrO_constantVals
	use AStrO_input
	use AStrO_output
	use AStrO_r_elementEqns
	use AStro_c_elementEqns
	use AStrO_r_designPropertyFunctions
	use AStrO_c_designPropertyFunctions
	use AStrO_solvers
	use AStrO_bookKeeping
	use AStrO_objective
	
	contains
	
	subroutine readAnyInput(fileName,commandTag)
	    implicit none
		
		character(len=128), intent(in) :: fileName
		character(len=32), intent(in) :: commandTag
		
		if(commandTag(1:15) .eq. '*readModelInput') then
		    call readModelInput(fileName)
			call readConstraints(fileName)
			call readLoads(fileName)
			call readInitialState(fileName)
		elseif(commandTag(1:10) .eq. '*readLoads') then
		    call readLoads(fileName)
		elseif(commandTag(1:16) .eq. '*readConstraints') then
		    call readConstraints(fileName)
		elseif(commandTag(1:17) .eq. '*readInitialState') then
		    call readInitialState(fileName)
		elseif(commandTag(1:19) .eq. '*readDesignVarInput') then
		    call readDesignVarInput(fileName)
		elseif(commandTag(1:20) .eq. '*readDesignVarValues') then
		    call readDesignVarValues(fileName)
		elseif(commandTag(1:19) .eq. '*readObjectiveInput') then
		    call readObjectiveInput(fileName)
		elseif(commandTag(1:19) .eq. '*readNodeResults') then
		    call readNodeResults(fileName)
		endif
		
	end subroutine readAnyInput
	
	subroutine processReadCommand(jobUnit,fileLine,iosVal)
	    implicit none
		
		integer, intent(in) :: jobUnit
		character(len=256), intent(out) :: fileLine
		integer, intent(out) :: iosVal
	
	    character(len=32) :: commandTag
		character(len=128) :: inputFileName
	    integer :: i1, i2
	
	    i1 = index(fileLine,'*')
		commandTag = fileLine(i1:i1+31)
	    read(jobUnit,'(A)',iostat=iosVal) fileLine(16:256)
		i1 = index(fileLine,'*')
		do while(i1 .eq. 0 .and. iosVal .eq. 0)
			i2 = index(fileLine,':')
			if(i2 .gt. 0) then
				if(fileLine(i2-8:i2) .eq. 'fileName:') then
					read(fileLine(i2+1:i2+128),*) inputFileName
					write(lfUnit,*) 'calling ', commandTag, 'file: ', inputFileName
					call readAnyInput(inputFileName,commandTag)
				endif
			endif
			read(jobUnit,'(A)',iostat=iosVal) fileLine(16:256)
			i1 = index(fileLine,'*')
		enddo
	
	end subroutine processReadCommand
	
	subroutine getThermalSolnLoad(buildMat)
	    implicit none
		
		integer, intent(in) :: buildMat
		
		integer :: numNds, dofPerNd, numIntDof, numIntPts
		integer :: dofTable(2,33)
		real*8 :: intPts(3,8), ipWt(8)
		real*8 :: locNds(3,10), globNds(3,10), orient(3,3), cMat(6,6), den, tExp(6), tCond(3,3), sHeat
		real*8 :: ABD(9,9), sMass(6,6), stExp(6), stCond(3,3), ssHeat
		real*8 :: bStiff(6,6), bMass(6,6), btExp(6), btCond(3,3), bsHeat
		real*8 :: temp(10), Tdot(10), pTemp(10), pTdot(10)
		real*8 :: disp(6,11), vel(6,11), acc(6,11), pDisp(6,11), pVel(6,11), pAcc(6,11)
		
		real*8 :: Rt(10), dRdT(10,10)
		integer :: i1, i2, i3, i4, i5, i6, Kgi, Kgj, inserted, eType
		
		if(buildMat .eq. 1) then
		    thermMat(:) = r_0
		endif
		do i1 = 1, numEls
		    call r_getElementData(numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt, &
				locNds, globNds, orient, cMat, den, tExp, tCond, sHeat, & 
				ABD, sMass, stExp, stCond, ssHeat, &
				bStiff, bMass, btExp, btCond, bsHeat, &
				temp,Tdot,disp,vel,acc,pTemp,pTdot,pDisp,pVel,pAcc,i1)
			eType = elementType(i1)
			call r_getElRt(Rt,dRdT,buildMat,numNds,numIntPts,intPts,ipWt, &
				locNds,orient,tCond,stCond,btCond,den,sHeat,ssHeat,bsHeat, &
				temp,pTemp,pTDot,eType)
			do i2 = 1, numNds
			    i3 = currentRank(elementList(i2,i1))
				thermalLoad(i3) = thermalLoad(i3) - Rt(i2)
			enddo
			if(buildMat .eq. 1) then
				do i2 = 1, numNds
					Kgi = currentRank(elementList(i2,i1))
					do i3 = 1, numNds
						Kgj = currentRank(elementList(i3,i1))
						i4 = thermMatRange(Kgi-1)+1
						inserted = 0
						do while(inserted .eq. 0 .and. i4 .lt. thermMatRange(Kgi))
							i6 = thermMatCols(i4)
							if(i6 .eq. 0) then
								thermMat(i4) = dRdT(i2,i3)
								thermMatCols(i4) = Kgj
								inserted = 1
							elseif(i6 .eq. Kgj) then
								thermMat(i4) = thermMat(i4) + dRdT(i2,i3)
								inserted = 1
							endif
							i4 = i4 + 1
						enddo
					enddo
				enddo			
			endif
		enddo
		
	end subroutine getThermalSolnLoad
	
	subroutine getThermalAppliedLoad(time)
	    implicit none
		
		real*8, intent(in) :: time
		
		integer :: dofTable(2,33), ndDof
		real*8 :: ndFrc(33), bdyFrc(6)
		real*8 :: trac(6), nDir(3), tTol	
		
		real*8 :: ndLd(7)
		character(len=16) :: lt
		integer :: i1, i2, i3, i4, i5, i6, i7, i8, i9, i10
		
		do i1 = 1, numLds
		    if(time .ge. loadsActTime(1,i1) .and. time .le. loadsActTime(2,i1)) then
			    lt = loadType(i1)
			    if(lt .eq. 'nodal') then
					do i2 = loadsRange(i1-1) + 1, loadsRange(i1)
						read(loadNodes(i2),*,err=55) i3
						i5 = currentRank(i3)
						thermalLoad(i5) = thermalLoad(i5) + inputLoads(7,i2)
						goto 63
55                      do i4 = 1, numNdSets
							if(ndSetName(i4) .eq. loadNodes(i2)) then
								do i5 = ndSetRange(i4-1)+1, ndSetRange(i4)
									i8 = currentRank(nodeSets(i5))
									thermalLoad(i8) = thermalLoad(i8) + inputLoads(7,i2)
								enddo
							endif
						enddo
63                      i4 = 1
				    enddo
				elseif(lt .eq. 'bodyHeatGen') then
					do i2 = loadsRange(i1-1) + 1, loadsRange(i1)
						read(loadNodes(i2),*,err=117) i3
						bdyFrc(1) = inputLoads(1,i2)
						bdyFrc(2:6) = r_0
						call r_getElBodyForce(ndFrc,dofTable,ndDof,bdyFrc,0,i3)
						do i4 = 1, ndDof
						    i5 = dofTable(1,i4)
							if(i5 .eq. 1) then
								i6 = currentRank(elementList(dofTable(2,i4),i3))
								thermalLoad(i6) = thermalLoad(i6) + ndFrc(i4)
							endif
						enddo
						goto 134
117                     do i4 = 1, numElSets
							if(elSetName(i4) .eq. loadNodes(i2)) then
								do i5 = elSetRange(i4-1)+1, elSetRange(i4)
									i6 = elementSets(i5)
								    bdyFrc(1) = inputLoads(1,i2)
									bdyFrc(2:6) = r_0
									call r_getElBodyForce(ndFrc,dofTable,ndDof,bdyFrc,0,i6)
									do i7 = 1, ndDof
										i8 = dofTable(1,i7)
										if(i8 .eq. 1) then
											i9 = currentRank(elementList(dofTable(2,i7),i6))
											thermalLoad(i9) = thermalLoad(i9) + ndFrc(i7)
										endif
									enddo
								enddo
							endif
						enddo
134                     i4 = 1
					enddo
                elseif(lt .eq. 'surfaceFlux') then
					do i2 = loadsRange(i1-1) + 1, loadsRange(i1)
						read(loadNodes(i2),*,err=147) i3
						trac(1) = inputLoads(5,i2)
						trac(2:6) = r_0
						do i4 = 1, 6
						    if(elementSurfaces(i4,i3) .eq. 1) then
							    nDir(:) = inputLoads(1:3,i2)
								tTol = inputLoads(4,i2)
							    call r_getElSurfaceTraction(ndFrc,dofTable,ndDof,trac,r_0,i3,i4,nDir,tTol)
							endif
						enddo
						do i4 = 1, ndDof
						    i5 = dofTable(1,i4)
							if(i5 .eq. 1) then
								i6 = currentRank(elementList(dofTable(2,i4),i3))
								thermalLoad(i6) = thermalLoad(i6) - ndFrc(i4)
							endif
						enddo
						goto 170
147                     do i4 = 1, numElSets
							if(elSetName(i4) .eq. loadNodes(i2)) then
								do i5 = elSetRange(i4-1)+1, elSetRange(i4)
									i6 = elementSets(i5)
									trac(1) = inputLoads(5,i2)
						            trac(2:6) = r_0
									do i7 = 1, 6
										if(elementSurfaces(i7,i6) .eq. 1) then
											nDir(:) = inputLoads(1:3,i2)
											tTol = inputLoads(4,i2)
											call r_getElSurfaceTraction(ndFrc,dofTable,ndDof,trac,r_0,i6,i7,nDir,tTol)
										endif
									enddo
									do i7 = 1, ndDof
										i8 = dofTable(1,i7)
										if(i8 .eq. 1) then
											i9 = currentRank(elementList(dofTable(2,i7),i6))
											thermalLoad(i9) = thermalLoad(i9) - ndFrc(i7)
										endif
									enddo
								enddo
							endif
						enddo
170                     i4 = 1
					enddo				
				endif
			endif
		enddo
		
		do i2 = 1, numNodes
		    call r_getNodeDLoad(time,ndLd,i2)
			i3 = currentRank(i2)
			thermalLoad(i3) = thermalLoad(i3) + ndLd(7)
		enddo
	
	end subroutine getThermalAppliedLoad
	
	subroutine getThermalConstraintLoad()
	    implicit none
		
	    integer :: i1, i2, i3
		real*8, allocatable :: tempV(:)
		
		allocate(tempV(thermMPCDim))
		
		tempV(:) = thermMPCRHS(:)
		do i1 = 1, thermMPCDim
		    do i2 = thermMPCMatRange(i1-1)+1,thermMPCMatRange(i1)
			    i3 = thermMPCMatCols(i2)
				if(i3 .ne. 0) then
				    tempV(i1) = tempV(i1) - thermMPCMat(i2)*nodeTemp(i3)
				endif
			enddo
		enddo
		
		do i1 = 1, thermMPCDim
		    do i2 = thermMPCMatRange(i1-1)+1,thermMPCMatRange(i1)
			    i3 = thermMPCMatCols(i2)
				if(i3 .ne. 0) then
				    thermalLoad(i3) = thermalLoad(i3) + thermMPCMat(i2)*tempV(i1)
				endif
			enddo
		enddo
		
		deallocate(tempV)
		
	end subroutine getThermalConstraintLoad
	
	subroutine scaleThermalMPC()
	    implicit none
		
		real*8 :: maxMat, rowMag
		integer :: i1, i2, i3, i4
		
		maxMat = r_0
		do i1 = 1, thermMatSize
		    rowMag = abs(thermMat(i1))
		    if(rowMag .gt. maxMat) then
			    maxMat = rowMag
			endif
		enddo
		
		do i1 = 1, thermMPCDim
		    i3 = thermMPCMatRange(i1-1)+1
			i4 = thermMPCMatRange(i1)
		    rowMag = r_0
			do i2 = i3, i4
			    rowMag = rowMag + thermMPCMat(i2)*thermMPCMat(i2)
			enddo
			rowMag = sqrt(rowMag)
			thermMPCMat(i3:i4) = (r_1/rowMag)*thermMPCMat(i3:i4)
			thermMPCRHS(i1) = (r_1/rowMag)*thermMPCRHS(i1)
		enddo
		
		thermMPCMat = 100d0*maxMat*thermMPCMat
		thermMPCRHS = 100d0*maxMat*thermMPCRHS
	
	end subroutine scaleThermalMPC
	
	subroutine getElasticSolnLoad(buildMat)
	    implicit none
		
		integer, intent(in) :: buildMat
		
		real*8 :: Ru(33), dRdU(33,33), intPts(3,8), ipWt(8)
		real*8 :: locNds(3,10), globNds(3,10), orient(3,3), cMat(6,6), tExp(6), tCond(3,3), sHeat, ABD(9,9), stExp(6)
		real*8 :: stCond(3,3), ssHeat, btCond(3,3), bsHeat
		real*8 :: bStiff(6,6), btExp(6), temp(10), disp(6,11), vel(6,11), acc(6,11), pDisp(6,11), pVel(6,11), pAcc(6,11)
		real*8 :: Tdot(10), pTemp(10), pTdot(10)
		real*8 :: statInOri(3,48), statdrIdrG(3,16), statRot(3), den, sMass(6,6), bMass(6,6)
		real*8 :: nnInv
		real*8 :: resVec(33), bVec(33), zVec(33)
		integer :: eType, numNds, dofPerNd, numIntDof, numIntPts, dofTable(2,33)
		integer :: i1, i2, i3, i4, i5, i6, kGi, kGj, inserted
		
		if(buildMat .eq. 1) then
		    elasticMat(:) = r_0
		endif
		do i1 = 1, numEls
		    call r_getElementData(numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt, &
				locNds, globNds, orient, cMat, den, tExp, tCond, sHeat, & 
				ABD, sMass, stExp, stCond, ssHeat, &
				bStiff, bMass, btExp, btCond, bsHeat, &
				temp,Tdot,disp,vel,acc,pTemp,pTdot,pDisp,pVel,pAcc,i1)
			!! -------------------
			! if(i1 .eq. 1) then
			    ! write(lfUnit,*) 'first element data'
				! write(lfUnit,*) 'numNds: ', numNds, 'dofPerNd: ', dofPerNd, 'numIntDof: ', numIntDof, 'numIntPts: ', numIntPts
				! write(lfUnit,*) 'intPts: '
				! write(lfUnit,*) intPts
				! write(lfUnit,*) 'ipWt'
				! write(lfUnit,*) ipWt
			    ! write(lfUnit,*) 'locNds: '
				! write(lfUnit,*) locNds
				! write(lfUnit,*) 'globNds: '
				! write(lfUnit,*) globNds
				! write(lfUnit,*) 'orient: '
				! write(lfUnit,*) orient
				! write(lfUnit,*) 'ABD: '
				! write(lfUnit,*) ABD
				! write(lfUnit,*) 'stExp: '
				! write(lfUnit,*) stExp
				! write(lfUnit,*) 'temp: '
				! write(lfUnit,*) temp
				! write(lfUnit,*) 'disp: '
				! write(lfUnit,*) disp
				! close(lfUnit)
				! open(unit=lfUnit, file='jobLogFile.txt', action='write', access='append')
			! endif
			!! ----------------------------
			eType = elementType(i1)
			statRot(:) = r_0
			do i2 = 1, numNds
			    statRot(:) = statRot(:) + disp(4:6,i2)
			enddo
			nnInv = numNds
			statRot(:) = (r_1/nnInv)*statRot(:)
			call r_getInstOrient(statInOri,statdrIdrG,orient,statRot)
			statInOri(:,4:48) = r_0
			statdrIdrG(:,4:16) = r_0
		    call r_getElRu(Ru,dRdU,buildMat,eType,numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt, &
				locNds,globNds,orient,cMat,tExp,ABD,stExp,bStiff,btExp,temp,disp,vel,acc,pDisp,pVel,pAcc, &
				statInOri,statdrIdrG,statRot,den,sMass,bMass)
			i3 = numNds*dofPerNd
			do i2 = 1, i3
			    i4 = dofTable(1,i2)
				i5 = elementList(dofTable(2,i2),i1)
				i6 = nDofIndex(i4,i5)
				elasticLoad(i6) = elasticLoad(i6) - Ru(i2)
			enddo
			if(numIntDof .ne. 0) then
			    i4 = intVecRange(i1-1) + 1
				i5 = intVecRange(i1)
				intElasticLoad(i4:i5) = intElasticLoad(i4:i5) - Ru(i3+1:i3+numIntDof)
			endif
			if(buildMat .eq. 1) then
			    ! -------------------
			    ! if(i1 .eq. 1) then
					! write(lfUnit,*) 'first element elastic matrix:'
					! i3 = numNds*dofPerNd + numIntDof
					! do i2 = 1, i3
						! write(lfUnit,*) dRdU(1:i3,i2)
					! enddo
					! close(lfUnit)
					! open(unit=lfUnit, file='jobLogFile.txt', action='write', access='append')
				! endif
				! ---------------------
				if(numIntDof .gt. 0) then
					i2 = numNds*dofPerNd + 1
					i3 = numNds*dofPerNd + numIntDof
					call r_condenseMat(dRdU,33,resVec,bVec,zVec,i2,i3)
					i4 = intMatRange(i1-1) + 1
					do i5 = 1, i3
						do i6 = i2, i3
							intElasticMat(i4) = dRdU(i6,i5)
							i4 = i4 + 1
						enddo
					enddo
				endif
				!! ---------------
				! if(i1 .eq. 1) then
					! dRdU(1,1) = 10000d0*dRdU(1,1)
					! dRdU(4,4) = 10000d0*dRdU(4,4)
					! dRdU(5,5) = 10000d0*dRdU(5,5)
					! dRdU(8,8) = 10000d0*dRdU(8,8)
					! dRdU(9,9) = 10000d0*dRdU(9,9)
					! dRdU(12,12) = 10000d0*dRdU(12,12)
					! dRdU(13,13) = 10000d0*dRdU(13,13)
					! dRdU(16,16) = 10000d0*dRdU(16,16)
					! dRdU(17,17) = 10000d0*dRdU(17,17)
					! dRdU(20,20) = 10000d0*dRdU(20,20)
					! dRdU(21,21) = 10000d0*dRdU(21,21)
					! dRdU(24,24) = 10000d0*dRdU(24,24)
					! bVec(:) = r_0
					! bVec(2) = r_1
					! bVec(3) = r_1
					! resVec(:) = r_0
					! write(lfUnit,*) 'first el mat:'
					! do i2 = 1, 24
					    ! write(lfUnit,*) dRdU(1:24,i2)
					! enddo
					! call rFactorMat(dRdU(1:24,1:24),24,24,0)
					! write(lfUnit,*) 'factored first el:'
					! do i2 = 1, 24
					    ! write(lfUnit,*) dRdU(1:24,i2)
					! enddo
					! call solveRFactor(resVec(1:24),dRdU(1:24,1:24),bVec(1:24),24,24,0)
					! write(lfUnit,*) 'single el result:'
					! do i2 = 1, numNds
					    ! write(lfUnit,*) resVec(i2), resVec(i2+4), resVec(i2+8), resVec(i2+12), resVec(i2+16), resVec(i2+20)
					! enddo
				! endif
				!! -------------------
				i2 = numNds*dofPerNd
				do i3 = 1, i2
					i5 = dofTable(1,i3)
					i6 = elementList(dofTable(2,i3),i1)
					kGi = nDofIndex(i5,i6)
					do i4 = 1, i2
						i5 = dofTable(1,i4)
						i6 = elementList(dofTable(2,i4),i1)
						kGj = nDofIndex(i5,i6)
						i5 = elMatRange(kGi-1) + 1
						inserted = 0
						do while(inserted .eq. 0 .and. i5 .le. elMatRange(kGi))
							i6 = elMatCols(i5)
							if(i6 .eq. 0) then
								elasticMat(i5) = dRdU(i3,i4)
								elMatCols(i5) = kGj
								inserted = 1
							elseif(i6 .eq. kGj) then
								elasticMat(i5) = elasticMat(i5) + dRdU(i3,i4)
								inserted = 1
							endif
							i5 = i5 + 1
						enddo
					enddo
				enddo
			endif
		enddo
		
	end subroutine getElasticSolnLoad
	
	subroutine getElasticAppliedLoad(time)
	    implicit none
		
		real*8, intent(in) :: time
		
		integer :: isGrav, dofTable(2,33), ndDof
		real*8 :: ndLd(7), nnInv, ndFrc(33), bdyFrc(6)
		real*8 :: trac(6), pressure, nDir(3), tTol
		real*8 :: center(3), axis(3), angVel
		character(len=16) :: lt
		integer :: i1, i2, i3, i4, i5, i6, i7, i8, i9, i10
		
		write(lfUnit,*) 'building elastic applied load'
		!close(lfUnit)
		!open(unit=lfUnit, file='jobLogFile.txt', action='write', access='append')
		
		do i1 = 1, numLds
		    if(time .ge. loadsActTime(1,i1) .and. time .le. loadsActTime(2,i1)) then
			    lt = loadType(i1)
			    if(lt .eq. 'nodal') then
				    write(lfUnit,*) 'building nodal load'
					write(lfUnit,*) 'numNdSets: ', numNdSets
					write(lfUnit,*) 'nDofIndex: '
					write(lfUnit,*) nDofIndex
				    close(lfUnit)
			        open(unit=lfUnit, file='jobLogFile.txt', action='write', access='append')
					do i2 = loadsRange(i1-1) + 1, loadsRange(i1)
						read(loadNodes(i2),*,err=67) i3
						do i4 = 1, 6
							i5 = nDofIndex(i4,i3)
							if(i5 .ne. 0) then
								elasticLoad(i5) = elasticLoad(i5) + inputLoads(i4,i2)
							endif
						enddo
						goto 75
67                      do i4 = 1, numNdSets
							if(ndSetName(i4) .eq. loadNodes(i2)) then
							    write(lfUnit,*) 'applying to node set ', ndSetName(i4)
						        close(lfUnit)
						        open(unit=lfUnit, file='jobLogFile.txt', action='write', access='append')
								do i5 = ndSetRange(i4-1)+1, ndSetRange(i4)
									i6 = nodeSets(i5)
									do i7 = 1, 6
										i8 = nDofIndex(i7,i6)
										if(i8 .ne. 0) then
											elasticLoad(i8) = elasticLoad(i8) + inputLoads(i7,i2)
										endif
									enddo
								enddo
							endif
						enddo
75                      i4 = 1
					enddo
				elseif(lt .eq. 'bodyForce' .or. lt .eq. 'gravitational' .or. lt .eq. 'centrifugal') then
				    if(lt .eq. 'bodyForce') then
					    isGrav = 0
				    else
					    isGrav = 1
					endif
					do i2 = loadsRange(i1-1) + 1, loadsRange(i1)
						read(loadNodes(i2),*,err=308) i3
						if(lt .eq. 'centrifugal') then
						    center(:) = inputLoads(1:3,i2)
							axis(:) = inputLoads(4:6,i2)
							angVel = inputLoads(7,i2)
						    call r_getElCentAcc(bdyFrc(1:3),center,axis,angVel,i3)
							bdyFrc(4:6) = r_0
						else
						    bdyFrc(:) = inputLoads(1:6,i2)
						endif
						call r_getElBodyForce(ndFrc,dofTable,ndDof,bdyFrc,isGrav,i3)
						do i4 = 1, ndDof
						    i5 = dofTable(1,i4)
							i6 = elementList(dofTable(2,i4),i3)
							i7 = nDofIndex(i5,i6)
							elasticLoad(i7) = elasticLoad(i7) + ndFrc(i4)
						enddo
						goto 328
308                     do i4 = 1, numElSets
							if(elSetName(i4) .eq. loadNodes(i2)) then
								do i5 = elSetRange(i4-1)+1, elSetRange(i4)
									i6 = elementSets(i5)
									if(lt .eq. 'centrifugal') then
									    center(:) = inputLoads(1:3,i2)
							            axis(:) = inputLoads(4:6,i2)
							            angVel = inputLoads(7,i2)
										call r_getElCentAcc(bdyFrc(1:3),center,axis,angVel,i6)
										bdyFrc(4:6) = r_0
									else
										bdyFrc(:) = inputLoads(1:6,i2)
									endif
									call r_getElBodyForce(ndFrc,dofTable,ndDof,bdyFrc,isGrav,i6)
									do i7 = 1, ndDof
										i8 = dofTable(1,i7)
										i9 = elementList(dofTable(2,i7),i6)
										i10 = nDofIndex(i8,i9)
										elasticLoad(i10) = elasticLoad(i10) + ndFrc(i7)
									enddo
								enddo
							endif
						enddo
328                     i4 = 1
					enddo
                elseif(lt .eq. 'surfacePressure' .or. lt .eq. 'surfaceTraction') then
					do i2 = loadsRange(i1-1) + 1, loadsRange(i1)
					    if(lt .eq. 'surfacePressure') then
						    pressure = inputLoads(5,i2)
							trac(:) = r_0
						else
						    pressure = r_0
							trac(1:3) = inputLoads(5:7,i2)
							trac(4:6) = r_0
						endif
						read(loadNodes(i2),*,err=367) i3
						do i4 = 1, 6
						    if(elementSurfaces(i4,i3) .eq. 1) then
							    nDir(:) = inputLoads(1:3,i2)
								tTol = inputLoads(4,i2)
							    call r_getElSurfaceTraction(ndFrc,dofTable,ndDof,trac,pressure,i3,i4,nDir,tTol)
							endif
						enddo
						do i4 = 1, ndDof
						    i5 = dofTable(1,i4)
							i6 = elementList(dofTable(2,i4),i3)
							i7 = nDofIndex(i5,i6)
							elasticLoad(i7) = elasticLoad(i7) + ndFrc(i4)
						enddo
						goto 387
367                     do i4 = 1, numElSets
							if(elSetName(i4) .eq. loadNodes(i2)) then
								do i5 = elSetRange(i4-1)+1, elSetRange(i4)
									i6 = elementSets(i5)
									do i7 = 1, 6
										if(elementSurfaces(i7,i6) .eq. 1) then
											nDir(:) = inputLoads(1:3,i2)
											tTol = inputLoads(4,i2)
											call r_getElSurfaceTraction(ndFrc,dofTable,ndDof,trac,pressure,i6,i7,nDir,tTol)
										endif
									enddo
									do i7 = 1, ndDof
										i8 = dofTable(1,i7)
										i9 = elementList(dofTable(2,i7),i6)
										i10 = nDofIndex(i8,i9)
										elasticLoad(i10) = elasticLoad(i10) + ndFrc(i7)
									enddo
								enddo
							endif
						enddo
387                     i4 = 1
					enddo			
				endif
			endif
		enddo
		
		write(lfUnit,*) 'building elastic dVar load'
		!close(lfUnit)
		!open(unit=lfUnit, file='jobLogFile.txt', action='write', access='append')
		
		do i2 = 1, numNodes
		    call r_getNodeDLoad(time,ndLd,i2)
			do i3 = 1, 6
			    i4 = nDofIndex(i3,i2)
				if(i4 .gt. 0) then
			        elasticLoad(i4) = elasticLoad(i4) + ndLd(i3)
				endif
			enddo
		enddo
		
		write(lfUnit,*) 'finished building elastic applied load'
		close(lfUnit)
	    open(unit=lfUnit, file='jobLogFile.txt', action='write', access='append')
		
	end subroutine getElasticAppliedLoad
	
	subroutine getElasticConstraintLoad()
	    implicit none
		
		integer :: i1, i2, i3
		real*8, allocatable :: tempV(:)
		
		allocate(tempV(elMPCDim))
		
		tempV(:) = elMPCRHS(:)
		do i1 = 1, elMPCDim
		    do i2 = elMPCMatRange(i1-1)+1,elMPCMatRange(i1)
			    i3 = elMPCMatCols(i2)
				if(i3 .ne. 0) then
				    tempV(i1) = tempV(i1) - elMPCMat(i2)*nodeDisp(i3)
				endif
			enddo
		enddo
		
		do i1 = 1, elMPCDim
		    do i2 = elMPCMatRange(i1-1)+1,elMPCMatRange(i1)
			    i3 = elMPCMatCols(i2)
				if(i3 .ne. 0) then
				    elasticLoad(i3) = elasticLoad(i3) + elMPCMat(i2)*tempV(i1)
				endif
			enddo
		enddo
		
		deallocate(tempV)
		
	end subroutine getElasticConstraintLoad
	
	subroutine scaleElasticMPC()
	    implicit none
		
		real*8 :: maxMat, rowMag
		integer :: i1, i2, i3, i4
		
		maxMat = r_0
		do i1 = 1, elMatSize
		    rowMag = abs(elasticMat(i1))
		    if(rowMag .gt. maxMat) then
			    maxMat = rowMag
			endif
		enddo
		
		do i1 = 1, elMPCDim
		    i3 = elMPCMatRange(i1-1)+1
			i4 = elMPCMatRange(i1)
		    rowMag = r_0
			do i2 = i3, i4
			    rowMag = rowMag + elMPCMat(i2)*elMPCMat(i2)
			enddo
			rowMag = sqrt(rowMag)
			elMPCMat(i3:i4) = (r_1/rowMag)*elMPCMat(i3:i4)
			elMPCRHS(i1) = (r_1/rowMag)*elMPCRHS(i1)
		enddo
		
		elMPCMat = 100d0*maxMat*elMPCMat
		elMPCRHS = 100d0*maxMat*elMPCRHS
		
	end subroutine scaleElasticMPC
	
	subroutine getDiagMassMat()
	    implicit none
		
		integer :: numNds, dofPerNd, numIntDof, numIntPts
		integer :: dofTable(2,33)
		real*8 :: intPts(3,8), ipWt(8)
		real*8 :: locNds(3,10), globNds(3,10), orient(3,3), cMat(6,6), den, tExp(6), tCond(3,3), sHeat
		real*8 :: ABD(9,9), sMass(6,6), sTELd(6), stCond(3,3), ssHeat
		real*8 :: bStiff(6,6), bMass(6,6), bTELd(6), btCond(3,3), bsHeat
		real*8 :: temp(10), Tdot(10), pTemp(10), pTdot(10)
		real*8 :: disp(6,11), vel(6,11), acc(6,11), pDisp(6,11), pVel(6,11), pAcc(6,11)		
		
		integer :: eType
		real*8 :: Rum(33), dRdA(33,33)
		real*8 :: statRot(3), statInOri(3,48), statdrIdrG(3,16)
		real*8 :: nnInv
		integer :: i1, i2, i3, i4, i5, i6

        if(.not. allocated(diagMassMat)) then
		    allocate(diagMassMat(elMatDim))
		endif

        diagMassMat(:) = r_0		
		do i1 = 1, numEls
		    call r_getElementData(numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt, &
				locNds, globNds, orient, cMat, den, tExp, tCond, sHeat, & 
				ABD, sMass, sTELd, stCond, ssHeat, &
				bStiff, bMass, bTELd, btCond, bsHeat, &
				temp,Tdot,disp,vel,acc,pTemp,pTdot,pDisp,pVel,pAcc,i1)
			eType = elementType(i1)
			statRot(:) = r_0
			do i2 = 1, numNds
			    statRot(:) = statRot(:) + disp(4:6,i2)
			enddo
			nnInv = numNds
			statRot(:) = (r_1/nnInv)*statRot(:)
			call r_getInstOrient(statInOri,statdrIdrG,orient,statRot)
			statInOri(:,4:48) = r_0
			statdrIdrG(:,4:16) = r_0
			acc(:,:) = r_0
			i3 = numNds*dofPerNd
			do i2 = 1, i3
			    i4 = dofTable(1,i2)
				i5 = dofTable(2,i2)
				acc(i4,i5) = r_1
			enddo
		    call r_getElRum(Rum,dRdA,0,eType,numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt, &
	            locNds, globNds, orient, statInOri, statdrIdrG, statRot, den, sMass, bMass, disp, acc)
			do i2 = 1, i3
			    i4 = dofTable(1,i2)
				i5 = dofTable(2,i2)
				i6 = nDofIndex(i4,i5)
				diagMassMat(i6) = diagMassMat(i6) + Rum(i2)
			enddo
		enddo
		
	end subroutine getDiagMassMat
	
	subroutine updateExternalRHS(extVec,extDim,intVec,intDim)
	    implicit none
		
		integer, intent(in) :: extDim, intDim
		real*8, intent(in) :: intVec(intDim)
		real*8, intent(out) :: extVec(extDim) 
		
		real*8 :: elIntMat(33,33), temp1(33), temp2(33), zVec(33)
		real*8 :: intPts(3,8), ipWt(8)
		integer :: numNds, dofPerNd, numIntDof, numIntPts, dofTable(2,33)
		integer :: i1, i2, i3, i4, i5, i6, nodalDof, totalDof, eType
	    
		do i1 = 1, numEls
		    call r_getElementProfile(numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt,i1)
			if(numIntDof .ne. 0) then
				nodalDof = numNds*dofPerNd
				totalDof = nodalDof + numIntDof
				i4 = intMatRange(i1-1) + 1
				do i2 = 1, totalDof
				    do i3 = nodalDof+1, numIntDof
					    elIntMat(i3,i2) = intElasticMat(i4)
						i4 = i4 + 1
					enddo
				enddo
				i4 = intVecRange(i1-1) + 1
				i5 = intVecRange(i1)
				temp2(nodalDof+1:totalDof) = intVec(i4:i5)
				call r_SolveLUxb(temp1, elIntMat, temp2, zVec, 33, nodalDof+1, totalDof)
				do i2 = 1, nodalDof
				    i4 = elementList(dofTable(2,i2),i1)
					i5 = dofTable(1,i2)
					i6 = nDofIndex(i5,i4)
				    do i3 = nodalDof+1, numIntDof
					    extVec(i6) = extVec(i6) - elIntMat(i3,i2)*temp1(i3)
					enddo
				enddo
			endif
		enddo
	
	end subroutine updateExternalRHS
	
	subroutine updateInternalSoln(extSoln,extDim,intSoln,intRHS,intDim)
	    implicit none
		
		integer, intent(in) :: extDim, intDim
		real*8, intent(in) :: extSoln(extDim), intRHS(intDim)
		real*8, intent(out) :: intSoln(intDim)
		
		real*8 :: elIntMat(33,33), temp1(33), temp2(33), zVec(33)
		real*8 :: intPts(3,8), ipWt(8)
		integer :: numNds, dofPerNd, numIntDof, numIntPts, dofTable(2,33)
		integer :: i1, i2, i3, i4, i5, i6, nodalDof, totalDof, eType
	    
		do i1 = 1, numEls
		    call r_getElementProfile(numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt,i1)
			if(numIntDof .ne. 0) then
				nodalDof = numNds*dofPerNd
				totalDof = nodalDof + numIntDof
				i4 = intMatRange(i1-1) + 1
				do i2 = 1, totalDof
				    do i3 = nodalDof+1, numIntDof
					    elIntMat(i3,i2) = intElasticMat(i4)
						i4 = i4 + 1
					enddo
				enddo
				i4 = intVecRange(i1-1) + 1
				i5 = intVecRange(i1)
				temp2(nodalDof+1:totalDof) = intRHS(i4:i5)
				do i2 = 1, nodalDof
				    i4 = elementList(dofTable(2,i2),i1)
					i5 = dofTable(1,i2)
					i6 = nDofIndex(i5,i4)
				    do i3 = nodalDof+1, numIntDof
					    temp2(i3) = temp2(i3) - elIntMat(i3,i2)*extSoln(i6)
					enddo
				enddo
				call r_SolveLUxb(temp1, elIntMat, temp2, zVec, 33, nodalDof+1, totalDof)
				i4 = intVecRange(i1-1) + 1
				i5 = intVecRange(i1)
				intSoln(i4:i5) = intSoln(i4:i5) + temp1(nodalDof+1:totalDof)
			endif
		enddo
		
	end subroutine updateInternalSoln
	
	subroutine loadInitialState()
	    implicit none
		
		integer :: i1, i2, i3
		
		prevTemp(:) = r_0
		prevTdot(:) = r_0
		prevDisp(:) = r_0
		prevVel(:) = r_0
		prevAcc(:) = r_0
		prevIntDisp(:) = r_0
		
		do i1 = 1, numNodes
		    i3 = currentRank(i1)
			prevTemp(i3) = initialTemp(i1)
			prevTdot(i3) = initialTdot(i1)
		    do i2 = 1, 6
			    i3 = nDofIndex(i2,i1)
				if(i3 .ne. 0) then
				    prevDisp(i3) = initialDisp(i2,i1)
					prevVel(i3) = initialVel(i2,i1)
					prevAcc(i3) = initialAcc(i2,i1)
				endif
			enddo
		enddo
		
	end subroutine loadInitialState
	
	subroutine getModelDimension(modDim)
	    implicit none
		
		real*8, intent(out) :: modDim
		
		real*8 :: xMin, xMax, yMin, yMax, zMin, zMax
		real*8 :: xLen, yLen, zLen, md2
		real*8 :: ndCrd(3)
		integer :: i1, i2, i3, i4
		
		xMin = 1e+100_8
		yMin = xMin
		zMin = xMin
		xMax = -xMin
		yMax = xMax
		zMax = xMax
		
		do i1 = 1, numNodes
		    ndCrd(:) = nodeList(:,i1)
			do i2 = ndToDRange(i1-1)+1, ndToDRange(i1)
				i3 = ndToD(i2)
				if(dCategory(i3) .eq. 'nodeCoord') then
					i4 = dComponent(i3)
					ndCrd(i4) = ndCrd(i4) + ndToCoef(i2)*r_dVec(i3)
				endif
			enddo
			if(ndCrd(1) .gt. xMax) then
			    xMax = ndCrd(1)
			endif
			if(ndCrd(1) .lt. xMin) then
			    xMin = ndCrd(1)
			endif
			if(ndCrd(2) .gt. yMax) then
			    yMax = ndCrd(2)
			endif
			if(ndCrd(2) .lt. yMin) then
			    yMin = ndCrd(2)
			endif
			if(ndCrd(3) .gt. zMax) then
			    zMax = ndCrd(3)
			endif
			if(ndCrd(3) .lt. zMin) then
			    zMin = ndCrd(3)
			endif
		enddo
		
		xLen = xMax - xMin
		yLen = yMax - yMin
		zLen = zMax - zMin
		
		md2 = xLen*xLen + yLen*yLen + zLen*zLen
		modDim = sqrt(md2)
		
	end subroutine getModelDimension
	
	subroutine setSolnToMode(solnVec,sVdim,modeNum,scaleFact)
	    implicit none
		
		integer, intent(in) :: modeNum, sVdim
		real*8, intent(in) :: scaleFact
		real*8, intent(out) :: solnVec(sVdim)
		
		real*8 :: modDim, eVMax, multFact, absVi
		integer :: i1, i2
		
		call getModelDimension(modDim)
		
		eVMax = r_0
		do i1 = 1, sVdim
		    absVi = abs(eigenModes(i1,modeNum))
		    if(absVi .gt. eVMax) then
			    eVMax = absVi
			endif
		enddo
		
		multFact = scaleFact*modDim/eVMax
		
		solnVec(:) = multFact*eigenModes(:,modeNum)
		
	end subroutine setSolnToMode
	
	subroutine getBucklingLoadFactors(bFact,numFact)
	    implicit none
		
		integer, intent(in) :: numFact
		real*8, intent(out) :: bFact(numFact)
		
		real*8, allocatable :: Ruk(:)
		real*8 :: numer, denom
		integer :: i1, i2, dynCopy
		
		allocate(Ruk(elMatDim))
		
		dynCopy = dynamic
		dynamic = 0
		elasticLoad(:) = r_0
		call getElasticSolnLoad(0)
		Ruk(:) = -elasticLoad(:)
		
		elasticLoad(:) = r_0
		call getElasticAppliedLoad(loadTime)
		
		do i1 = 1, numFact
		    numer = r_0
			denom = r_0
			do i2 = 1, elMatDim
			    numer = numer + Ruk(i2)*eigenModes(i2,i1)
				denom = denom + elasticLoad(i2)*eigenModes(i2,i1)
			enddo
			if(abs(denom) .gt. r_0) then
			    bFact(i1) = numer/denom
			else
			    bFact(i1) = 1e+100_8
			endif
		enddo
		
		dynamic = dynCopy
		deallocate(Ruk)
		
	end subroutine getBucklingLoadFactors
	
	subroutine stepSolve(time)
	    implicit none
		
		real*8, intent(in) :: time
		
		real*8 :: delUNorm, mag
		integer :: i1, i2, maxNLIt, numVecs
		
		if(solveThermal .eq. 1) then
		    if(thermMatDim .lt. 200) then
				numVecs = thermMatDim/2
			else
				numVecs = 100
			endif
			thermalLoad(:) = r_0
			call getThermalSolnLoad(0)
			call getThermalAppliedLoad(time)
			call getThermalConstraintLoad()
			call gMRes(delDisp(1:numNodes),thermMat,thermMatCols,thermMatRange,thermMatSize,thermMatDim, &
			    thermMPCMat,thermMPCMatCols,thermMPCMatRange,thermMPCSize,thermMPCDim,thermMatLT,thermMatLTRange, &
				thermMatLTSize,thermalLoad,numVecs,thermMatDim)
			nodeTemp(:) = nodeTemp(:) + delDisp(1:numNodes)		
		endif
	
	    if(solveElastic .eq. 1) then
			if(nLGeom .eq. 0) then
				maxNLIt = 1
			else
			    maxNLIt = 20
			endif
			if(elMatDim .lt. 200) then
				numVecs = elMatDim/2
			else
				numVecs = 100
			endif
			i1 = 0
			delUNorm = r_1
			do while(i1 .lt. maxNLIt .and. delUNorm .gt. 1e-12)
			    elasticLoad(:) = r_0
				intElasticLoad(:) = r_0
			    call getElasticSolnLoad(nLGeom)
				write(lfUnit,*) 'finished building elastic solution load'
				!close(lfUnit)
			    !open(unit=lfUnit, file='jobLogFile.txt', action='write', access='append')
				call getElasticAppliedLoad(time)
				write(lfUnit,*) 'finished building elastic applied load'
				!close(lfUnit)
			    !open(unit=lfUnit, file='jobLogFile.txt', action='write', access='append')
				call getElasticConstraintLoad()
				write(lfUnit,*) 'finished building elastic constraint load'
				write(lfUnit,*) 'elasticLoad:'
				do i2 = 1, elMatDim
				    write(lfUnit,*) elasticLoad(i2)
				enddo
				!close(lfUnit)
			    !open(unit=lfUnit, file='jobLogFile.txt', action='write', access='append')
				if(nLGeom .ne. 0) then
					call convertToLTri(elMatLT,elMatLTRange,elMatLTSize,elasticMat,elMatCols,elMatRange, &
						 elMatSize,elMatDim,elMPCMat,elMPCMatCols,elMPCMatRange,elMPCSize,elMPCDim)
					call getSparseLDLFact(elMatLT,elMatLTSize,elMatLTRange,elMatDim)
				endif
				write(lfUnit,*) 'finished factoring'
				!close(lfUnit)
			    !open(unit=lfUnit, file='jobLogFile.txt', action='write', access='append')
				if(intVecSize .gt. 0) then
				    call updateExternalRHS(elasticLoad,elMatDim,intElasticLoad,intVecSize)
					write(lfUnit,*) 'finished updating external RHS'
				    close(lfUnit)
			        open(unit=lfUnit, file='jobLogFile.txt', action='write', access='append')
				endif
				write(lfUnit,*) 'calling gMRes, elastic solution'
			    !close(lfUnit)
			    !open(unit=lfUnit, file='jobLogFile.txt', action='write', access='append')
				!call gMRes(delDisp,elasticMat,elMatCols,elMatRange,elMatSize,elMatDim,elMPCMat,elMPCMatCols, &
				!    elMPCMatRange,elMPCSize,elMPCDim,elMatLT,elMatLTRange,elMatLTSize,elasticLoad,numVecs,elMatDim)
				call solveRFactor(delDisp,aFull,elasticLoad,elMatDim,elMatDim,0)
				if(intVecSize .gt. 0) then
				    call updateInternalSoln(delDisp,elMatDim,internalDisp,intElasticLoad,intVecSize)
				endif
				nodeDisp(:) = nodeDisp(:) + delDisp(:)
				delUNorm = r_0
				do i2 = 1, elMatDim
				    mag = abs(delDisp(i2))
				    if(mag .gt. delUNorm) then
					    delUNorm = mag
					endif
				enddo
				i1 = i1 + 1
			enddo
		endif
	
	end subroutine stepSolve
	
	subroutine solve()
	    implicit none
		
		character(len=128) :: outfileName
		real*8 :: time, c1, c2
		integer :: i1, i2
		
		if(.not. allocated(currentRank)) then
		    call analysisPrep()
		endif
		
		write(lfUnit,*) 'finished analysisPrep'
		close(lfUnit)
		open(unit=lfUnit, file='jobLogFile.txt', action='write', access='append')
		
		call loadInitialState()
		
		if(solveThermal .eq. 1) then
		    call getThermalSolnLoad(1)
			call scaleThermalMPC()
			call convertToLTri(thermMatLT,thermMatLTRange,thermMatLTSize,thermMat,thermMatCols,thermMatRange, &
			    thermMatSize,thermMatDim,thermMPCMat,thermMPCMatCols,thermMPCMatRange,thermMPCSize,thermMPCDim)
			call getSparseLDLFact(thermMatLT,thermMatLTSize,thermMatLTRange,thermMatDim)
		endif
		
		if(solveElastic .eq. 1) then
		    call getElasticSolnLoad(1)
			write(lfUnit,*) 'finished building elastic matrix'
			close(lfUnit)
			open(unit=lfUnit, file='jobLogFile.txt', action='write', access='append')
			call scaleElasticMPC()
			if(nLGeom .eq. 0) then
			    !! --------------------------
				outfileName = 'sparseStiffness.txt'
				call writeSparseMatrix(outfileName,elasticMat,elMatCols,elMatRange,elMatSize,elMatDim)
				outfileName = 'elasticMPC.txt'
				call writeSparseMatrix(outfileName,elMPCMat,elMPCMatCols,elMPCMatRange,elMPCSize,elMPCDim)
				!! ---------------------------------
				
				! call convertToLTri(elMatLT,elMatLTRange,elMatLTSize,elasticMat,elMatCols,elMatRange, &
					 ! elMatSize,elMatDim,elMPCMat,elMPCMatCols,elMPCMatRange,elMPCSize,elMPCDim)
			    allocate(aFull(elMatDim,elMatDim))
			    call convertToFullPop(aFull,elasticMat,elMatCols,elMatRange,elMatSize,elMatDim,elMPCMat, &
				    elMPCMatCols,elMPCMatRange,elMPCSize,elMPCDim)
				
				!! ------------------------
				!outfileName = 'elasticLTStiffess.txt'
				!call writeLowerTriMatrix(outfileName,elMatLT,elMatLTRange,elMatLTSize,elMatDim)
				!! -------------------------
				
				!call getSparseLDLFact(elMatLT,elMatLTSize,elMatLTRange,elMatDim)
				
				call rFactorMat(aFull,elMatDim,elMatDim,0)
				
				write(lfUnit,*) 'Factored diagonal terms: '
				do i1 = 1, elMatDim
				    write(lfUnit,*) aFull(i1,i1)
				enddo
				
				!! ------------------------
				!outfileName = 'elasticLTFactored.txt'
				!call writeLowerTriMatrix(outfileName,elMatLT,elMatLTRange,elMatLTSize,elMatDim)
				!! -------------------------
			endif
		endif
		
		if(dynamic .eq. 0) then
		    call stepSolve(loadTime)
		else
		    c1 = r_p5/nMBeta
			c2 = r_2/(delT*delT)
		    i1 = 0
			time = r_0
			if(writeSolnHist .eq. 1) then
				call writeBinarySolution(i1)
			endif
		    do while(time .lt. simPeriod)
			    time = time + delT
				i1 = i1 + 1
				call stepSolve(time)
				do i2 = 1, numNodes
				    nodeTdot(i2) = (r_1/nMGamma)*((r_1/delT)*(nodeTemp(i2) - prevTemp(i2)) - (r_1-nMGamma)*prevTdot(i2))
				enddo
				do i2 = 1, elMatDim
					nodeAcc(i2) = c1*(c2*(nodeDisp(i2)-prevDisp(i2)-delT*prevVel(i2)) - (r_1 - r_2*nMBeta)*prevAcc(i2))
					nodeVel(i2) = prevVel(i2) + delT*((r_1-nMGamma)*prevAcc(i2) + nMGamma*nodeAcc(i2))
				enddo
				prevTemp(:) = nodeTemp(:)
				prevTdot(:) = nodeTdot(:)
				prevDisp(:) = nodeDisp(:)
				if(intVecSize .gt. 0) then
				    prevIntDisp(:) = internalDisp(:)
				endif
				prevVel(:) = nodeVel(:)
				prevAcc(:) = nodeAcc(:)
				if(writeSolnHist .eq. 1) then
					call writeBinarySolution(i1)
				endif
			enddo
			numTSteps = i1
		endif
		
		deallocate(aFull)
		
	end subroutine solve
	
	subroutine getCompleteObjective()
	    implicit none
		
		real*8 :: time
		integer :: i1, i2, i3, errFlag
		
		objVal = r_0
		
		if(dynamic .eq. 0) then
		    call getObjective(loadTime)
		else
		    i1 = numTSteps
		    call readBinarySolution(numTSteps,errFlag)
			if(errFlag .eq. 1) then
			    goto 896
			endif
			do i1 = numTSteps-1, 0, -1
			    nodeTemp(:) = prevTemp(:)
				nodeTdot(:) = prevTdot(:)
			    nodeDisp(:) = prevDisp(:)
				if(intVecSize .gt. 0) then
				    internalDisp(:) = prevIntDisp(:)
				endif
				nodeVel(:) = prevVel(:)
				nodeAcc(:) = prevAcc(:)
				call readBinarySolution(i1,errFlag)
				if(errFlag .eq. 1) then
				   goto 896
				endif
				time = delT*(i1+1)
				call getObjective(time)
			enddo
			return
896			write(lfUnit,*) 'Error: could not read the solution history file for time step ', i1
			write(lfUnit,*) 'Be sure to run a dynamic analysis with *solve with the saveSolnHistory option set to yes'
			write(lfUnit,*) 'if objective and sensitivities are desired for dynamic problems.'
		endif
		
	end subroutine getCompleteObjective
	
	subroutine augmentdLdu()
	    implicit none
		
		real*8 :: multFact
		integer :: i1
		
		multFact = 1e-20_8
		
		nodeTemp(:) = multFact*nodeTemp(:)
		nodeTdot(:) = multFact*nodeTdot(:)
		nodeDisp(:) = multFact*nodeDisp(:)
		nodeVel(:) = multFact*nodeVel(:)
		nodeAcc(:) = multFact*nodeAcc(:)
		prevTemp(:) = multFact*prevTemp(:)
		prevTdot(:) = multFact*prevTdot(:)
		prevDisp(:) = multFact*prevDisp(:)
		prevVel(:) = multFact*prevVel(:)
		prevAcc(:) = multFact*prevAcc(:)
		if(intVecSize .gt. 0) then
		    internalDisp(:) = multFact*internalDisp(:)
		    prevIntDisp(:) = multFact*prevIntDisp(:)
		endif
		
		if(solveThermal .eq. 1) then
		    !! Block 1
		    swapVec(1:numNodes) = prevTemp(:)
			prevTemp(:) = tempAdj(:)
			thermalLoad(:) = r_0
			call getThermalSolnLoad(0)
			dLdt(:) = dLdt(:) + thermalLoad(:) + tdotAdj(:)
			prevTemp(:) = swapVec(1:numNodes)
			
			!! Block 2
			swapVec(1:numNodes) = prevTdot(:)
			prevTdot(:) = tempAdj(:)
			thermalLoad(:) = r_0
			call getThermalSolnLoad(0)
			dLdtdot(:) = dLdtdot(:) + thermalLoad(:) + delT*(r_1-nMGamma)*tdotAdj(:)
			prevTdot(:) = swapVec(1:numNodes)
		endif
		
		if(solveElastic .eq. 1) then
		    !! Block 3
			swapVec(:) = prevDisp(:)
			prevDisp(:) = dispAdj(:)
			if(intVecSize .gt. 0) then
			    intswapVec(:) = prevIntDisp(:)
			    prevIntDisp(:) = intDispAdj(:)
			endif
			elasticLoad(:) = r_0
			intElasticLoad(:) = r_0
			call getElasticSolnLoad(0)
			dLdu(:) = dLdu(:) + elasticLoad(:) + accAdj(:)
			prevDisp(:) = swapVec(:)
			if(intVecSize .gt. 0) then
			    intdLdu(:) = intdLdu(:) + intElasticLoad(:)
			    prevIntDisp(:) = intswapVec(:)
			endif
			
			!! Block 4
			swapVec(:) = prevAcc(:)
			prevAcc(:) = dispAdj(:)
			elasticLoad(:) = r_0
			intElasticLoad(:) = r_0
			call getElasticSolnLoad(0)
			
			do i1 = 1, elMatDim
			    dLda(i1) = dLda(i1) + elasticLoad(i1) + r_p5*delT*delT*(r_1-r_2*nMBeta)*accAdj(i1) + delT*(r_1-nMGamma)*velAdj(i1)
			enddo
			prevAcc(:) = swapVec(:)
			
			!! Block 5
			swapVec(:) = prevVel(:)
			prevVel(:) = dispAdj(:)
			elasticLoad(:) = r_0
			intElasticLoad(:) = r_0
			call getElasticSolnLoad(0)
			dLdv(:) = dLdv(:) + elasticLoad(:) + delT*accAdj(:) + velAdj(:)
			prevVel(:) = swapVec(:)
		endif
		
		multFact = r_1/multFact
		
		nodeTemp(:) = multFact*nodeTemp(:)
		nodeTdot(:) = multFact*nodeTdot(:)
		nodeDisp(:) = multFact*nodeDisp(:)
		nodeVel(:) = multFact*nodeVel(:)
		nodeAcc(:) = multFact*nodeAcc(:)
		prevTemp(:) = multFact*prevTemp(:)
		prevTdot(:) = multFact*prevTdot(:)
		prevDisp(:) = multFact*prevDisp(:)
		prevVel(:) = multFact*prevVel(:)
		prevAcc(:) = multFact*prevAcc(:)
		internalDisp(:) = multFact*internalDisp(:)
		prevIntDisp(:) = multFact*prevIntDisp(:)
		
	end subroutine augmentdLdu
	
	subroutine solveDynamicAdjoint()
	    implicit none
		
		integer :: elNum
		integer :: numNds, dofPerNd, numIntDof, numIntPts
		integer :: dofTable(2,33)
		real*8 :: intPts(3,8), ipWt(8)
		real*8 :: locNds(3,10), globNds(3,10), orient(3,3), cMat(6,6), den, tExp(6), tCond(3,3), sHeat
		real*8 :: ABD(9,9), sMass(6,6), stExp(6), stCond(3,3), ssHeat
		real*8 :: bStiff(6,6), bMass(6,6), btExp(6), btCond(3,3), bsHeat
		real*8 :: temp(10), Tdot(10), pTemp(10), pTdot(10)
		real*8 :: disp(6,11), vel(6,11), acc(6,11), pDisp(6,11), pVel(6,11), pAcc(6,11)
		
		real*8 :: Ruk(33), dRdU(33,33), dRdT(33,10), stEn, dSEdT(10)
		
		integer :: i1, i2, i3, i4, i5, i6, i7, eType, numVecs
		
		if(solveElastic .eq. 1) then
			velAdj(:) = dLdv(:)
			accAdj(:) = (-r_1/(delT*delT*nMBeta))*(dLda(:) + delT*nMGamma*velAdj(:))
			
			dLdu(:) = dLdu(:) - accAdj(:)
			call getElasticSolnLoad(1)
			call scaleElasticMPC()
			call convertToLTri(elMatLT,elMatLTRange,elMatLTSize,elasticMat,elMatCols,elMatRange, &
				 elMatSize,elMatDim,elMPCMat,elMPCMatCols,elMPCMatRange,elMPCSize,elMPCDim)
			call getSparseLDLFact(elMatLT,elMatLTSize,elMatLTRange,elMatDim)
			if(intVecSize .gt. 0) then
			    call updateExternalRHS(dLdu,elMatDim,intdLdu,intVecSize)
			endif
			if(elMatDim .lt. 200) then
				numVecs = elMatDim/2
			else
				numVecs = 100
			endif
			call gMRes(dispAdj,elasticMat,elMatCols,elMatRange,elMatSize,elMatDim,elMPCMat,elMPCMatCols, &
				    elMPCMatRange,elMPCSize,elMPCDim,elMatLT,elMatLTRange,elMatLTSize,dLdu,numVecs,elMatDim)
			if(intVecSize .gt. 0) then
			    call updateInternalSoln(dispAdj,elMatDim,intDispAdj,intdLdu,intVecSize)
			endif
        endif

        if(solveThermal .eq. 1) then
		    tdotAdj(:) = (-r_1/(delT*nMGamma))*dLdtdot(:)
			
			if(solveElastic .eq. 1) then
			    do i1 = 1, numEls
				    call r_getElementData(numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt, &
						locNds, globNds, orient, cMat, den, tExp, tCond, sHeat, & 
						ABD, sMass, stExp, stCond, ssHeat, &
						bStiff, bMass, btExp, btCond, bsHeat, &
						temp,Tdot,disp,vel,acc,pTemp,pTdot,pDisp,pVel,pAcc,i1)
					eType = elementType(i1)
					call r_getElRuk(Ruk,dRdU,dRdT,stEn,dSEdT,1,eType,numNds,dofPerNd,numIntDof,numIntPts,dofTable, &
	                    intPts,ipWt,locNds,globNds,orient,cMat,tExp,ABD,stExp,bStiff,btExp,temp,disp)
					i2 = numNds*dofPerNd
					do i3 = 1, i2
					    i4 = dofTable(1,i3)
						i5 = elementList(dofTable(2,i3),i1)
						i6 = nDofIndex(i4,i5)
						do i4 = 1, numNds
						    i5 = currentRank(elementList(i4,i1))
						    dLdt(i5) = dLdt(i5) - dRdT(i3,i4)*dispAdj(i6)
						enddo
					enddo
					if(numIntDof .gt. 0) then
					    i2 = numNds*dofPerNd + 1
						i3 = numNds*dofPerNd + numIntDof
						i4 = intVecRange(i1-1) + 1
						do i5 = i2, i3
						    do i6 = 1, numNds
							    i7 = currentRank(elementList(i6,i1))
							    dLdt(i7) = dLdt(i7) - dRdT(i5,i6)*intDispAdj(i4)
							enddo
							i4 = i4 + 1
						enddo
					endif
				enddo
			endif
			
			dLdt(:) = dLdt(:) - tdotAdj(:)
			
			call getThermalSolnLoad(1)
			call scaleThermalMPC()
			call convertToLTri(thermMatLT,thermMatLTRange,thermMatLTSize,thermMat,thermMatCols,thermMatRange, &
			    thermMatSize,thermMatDim,thermMPCMat,thermMPCMatCols,thermMPCMatRange,thermMPCSize,thermMPCDim)
			call getSparseLDLFact(thermMatLT,thermMatLTSize,thermMatLTRange,thermMatDim)
			if(thermMatDim .lt. 200) then
				numVecs = thermMatDim/2
			else
				numVecs = 100
			endif
			call gMRes(tempAdj(:),thermMat,thermMatCols,thermMatRange,thermMatSize,thermMatDim, &
			    thermMPCMat,thermMPCMatCols,thermMPCMatRange,thermMPCSize,thermMPCDim,thermMatLT,thermMatLTRange, &
				thermMatLTSize,dLdt,numVecs,thermMatDim)
        endif		
		
	end subroutine solveDynamicAdjoint
	
	subroutine getdRthermaldD(time,dVarNum)
	    implicit none
		
		integer, intent(in) :: dVarNum
		real*8, intent(in) :: time
		
		integer :: eType
		integer :: numNds, dofPerNd, numIntDof, numIntPts
		integer :: dofTable(2,33)
		complex*16 :: intPts(3,8), ipWt(8)
		complex*16 :: locNds(3,10), globNds(3,10), orient(3,3), cMat(6,6), den, tExp(6), tCond(3,3), sHeat
		complex*16 :: ABD(9,9), sMass(6,6), sTELd(6), stCond(3,3), ssHeat
		complex*16 :: bStiff(6,6), bMass(6,6), bTELd(6), btCond(3,3), bsHeat
		complex*16 :: temp(10), Tdot(10), pTemp(10), pTdot(10)
		complex*16 :: disp(6,11), vel(6,11), acc(6,11), pDisp(6,11), pVel(6,11), pAcc(6,11)

		integer :: ndDof
		complex*16 :: ndFrc(33), bdyFrc(6)
		complex*16 :: trac(6), nDir(3), tTol
        character(len=16) :: lt		
		
		complex*16 :: Rt(10), dRdt(10,10), ndLd(7)
		integer :: i1, i2, i3, i4, i5, i6, i7, i8, i9
		
		dRtdD(:) = r_0
		
	    elDepends(:) = 0
		do i1 = dToElCRange(dVarNum-1)+1, dToElCRange(dVarNum)
		    i2 = dToElComp(i1)
			if(i2 .ne. 0) then
			    elDepends(i2) = 1
			endif
		enddo
		
		c_dVec(:) = c_1*r_dVec(:)
		c_dVec(dVarNum) = c_dVec(dVarNum) + compStep
		
		do i1 = dToElCRange(dVarNum-1) + 1, dToElCRange(dVarNum)
		    i2 = dToElComp(i1)
			call c_getElementData(numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt, &
				locNds, globNds, orient, cMat, den, tExp, tCond, sHeat, & 
				ABD, sMass, sTELd, stCond, ssHeat, &
				bStiff, bMass, bTELd, btCond, bsHeat, &
	            temp,Tdot,disp,vel,acc,pTemp,pTdot,pDisp,pVel,pAcc,i2)
			eType = elementType(i2)
			call c_getElRt(Rt,dRdt,0,numNds,numIntPts,intPts,ipWt, &
				locNds,orient,tCond,stCond,btCond,den,sHeat,ssHeat,bsHeat, &
				temp,pTemp,pTDot,eType)
			do i3 = 1, numNds
			    i4 = currentRank(elementList(i3,i2))
				dRtdD(i4) = dRtdD(i4) + compStepInv*imag(Rt(i3))
			enddo
		enddo
		
		do i1 = 1, numLds
		    if(time .ge. loadsActTime(1,i1) .and. time .le. loadsActTime(2,i1)) then
			    lt = loadType(i1)
				if(lt .eq. 'bodyHeatGen') then
					do i2 = loadsRange(i1-1) + 1, loadsRange(i1)
						read(loadNodes(i2),*,err=1082) i3
						if(elDepends(i3) .eq. 1) then
							bdyFrc(1) = inputLoads(1,i2)
							bdyFrc(2:6) = c_0
							call c_getElBodyForce(ndFrc,dofTable,ndDof,bdyFrc,0,i3)
							do i4 = 1, ndDof
								i5 = dofTable(1,i4)
								if(i5 .eq. 1) then
									i6 = currentRank(elementList(dofTable(2,i4),i3))
									dRtdD(i6) = dRtdD(i6) - compStepInv*imag(ndFrc(i4))
								endif
							enddo
						endif
						goto 1099
1082                    do i4 = 1, numElSets
							if(elSetName(i4) .eq. loadNodes(i2)) then
								do i5 = elSetRange(i4-1)+1, elSetRange(i4)
									i6 = elementSets(i5)
									if(elDepends(i6) .eq. 1) then
										bdyFrc(1) = inputLoads(1,i2)
										bdyFrc(2:6) = c_0
										call c_getElBodyForce(ndFrc,dofTable,ndDof,bdyFrc,0,i6)
										do i7 = 1, ndDof
											i8 = dofTable(1,i7)
											if(i8 .eq. 1) then
												i9 = currentRank(elementList(dofTable(2,i7),i6))
												dRtdD(i9) = dRtdD(i9) - compStepInv*imag(ndFrc(i7))
											endif
										enddo
									endif
								enddo
							endif
						enddo
1099                    i4 = 1
					enddo
                elseif(lt .eq. 'surfaceFlux') then
					do i2 = loadsRange(i1-1) + 1, loadsRange(i1)
						read(loadNodes(i2),*,err=1121) i3
						if(elDepends(i3) .eq. 1) then
							trac(1) = inputLoads(5,i2)
							trac(2:6) = c_0
							do i4 = 1, 6
								if(elementSurfaces(i4,i3) .eq. 1) then
									nDir(:) = inputLoads(1:3,i2)
									tTol = inputLoads(4,i2)
									call c_getElSurfaceTraction(ndFrc,dofTable,ndDof,trac,c_0,i3,i4,nDir,tTol)
								endif
							enddo
							do i4 = 1, ndDof
								i5 = dofTable(1,i4)
								if(i5 .eq. 1) then
									i6 = currentRank(elementList(dofTable(2,i4),i3))
									dRtdD(i6) = dRtdD(i6) + compStepInv*imag(ndFrc(i4))
								endif
							enddo
						endif
						goto 1144
1121                    do i4 = 1, numElSets
							if(elSetName(i4) .eq. loadNodes(i2)) then
								do i5 = elSetRange(i4-1)+1, elSetRange(i4)
									i6 = elementSets(i5)
									if(elDepends(i6) .eq. 1) then
										trac(1) = inputLoads(5,i2)
										trac(2:6) = c_0
										do i7 = 1, 6
											if(elementSurfaces(i7,i6) .eq. 1) then
												nDir(:) = inputLoads(1:3,i2)
												tTol = inputLoads(4,i2)
												call c_getElSurfaceTraction(ndFrc,dofTable,ndDof,trac,c_0,i6,i7,nDir,tTol)
											endif
										enddo
										do i7 = 1, ndDof
											i8 = dofTable(1,i7)
											if(i8 .eq. 1) then
												i9 = currentRank(elementList(dofTable(2,i7),i6))
												dRtdD(i9) = dRtdD(i9) + compStepInv*imag(ndFrc(i7))
											endif
										enddo
									endif
								enddo
							endif
						enddo
1144                    i4 = 1
					enddo				
				endif
			endif
		enddo	
		
		do i1 = dToNdRange(dVarNum-1) + 1, dToNdRange(dVarNum)
		    i2 = dToNd(i1)
		    call c_getNodeDLoad(c_1*time,ndLd,i2)
			i3 = currentRank(i2)
			dRtdD(i3) = dRtdD(i3) - compStepInv*imag(ndLd(7))
		enddo
		
	end subroutine getdRthermaldD
	
	subroutine getdRelasticdD(time,dVarNum)
	    implicit none
		
		integer, intent(in) :: dVarNum
		real*8, intent(in) :: time
		
		integer :: eType
		integer :: numNds, dofPerNd, numIntDof, numIntPts
		integer :: dofTable(2,33)
		complex*16 :: intPts(3,8), ipWt(8)
		complex*16 :: locNds(3,10), globNds(3,10), orient(3,3), cMat(6,6), den, tExp(6), tCond(3,3), sHeat
		complex*16 :: ABD(9,9), sMass(6,6), stExp(6), stCond(3,3), ssHeat
		complex*16 :: bStiff(6,6), bMass(6,6), btExp(6), btCond(3,3), bsHeat
		complex*16 :: temp(10), Tdot(10), pTemp(10), pTdot(10)
		complex*16 :: disp(6,11), vel(6,11), acc(6,11), pDisp(6,11), pVel(6,11), pAcc(6,11)

  		integer :: isGrav, ndDof
		complex*16 :: ndFrc(33), bdyFrc(6)
		complex*16 :: trac(6), pressure, nDir(3), tTol
		complex*16 :: center(3), axis(3), angVel
		character(len=16) :: lt
		
		complex*16 :: Ru(33), dRdU(33,33), ndLd(7)
		complex*16 :: statRot(3), nnInv, statInOri(3,48), statdrIdrG(3,16)
		integer :: i1, i2, i3, i4, i5, i6, i7, i8, i9, i10
		
	    dRudD(:) = r_0
		if(intVecSize .gt. 0) then
	        intdRudD(:) = r_0
		endif
		
		elDepends(:) = 0
		do i1 = dToElCRange(dVarNum-1)+1, dToElCRange(dVarNum)
		    i2 = dToElComp(i1)
			if(i2 .ne. 0) then
			    elDepends(i2) = 1
			endif
		enddo
		
		c_dVec(:) = c_1*r_dVec(:)
		c_dVec(dVarNum) = c_dVec(dVarNum) + compStep
		
		do i1 = dToElCRange(dVarNum-1) + 1, dToElCRange(dVarNum)
		    i2 = dToElComp(i1)
			call c_getElementData(numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt, &
				locNds, globNds, orient, cMat, den, tExp, tCond, sHeat, & 
				ABD, sMass, stExp, stCond, ssHeat, &
				bStiff, bMass, btExp, btCond, bsHeat, &
	            temp,Tdot,disp,vel,acc,pTemp,pTdot,pDisp,pVel,pAcc,i2)
			eType = elementType(i2)
			statRot(:) = c_0
			do i2 = 1, numNds
			    statRot(:) = statRot(:) + disp(4:6,i2)
			enddo
			nnInv = c_1*numNds
			statRot(:) = (c_1/nnInv)*statRot(:)
			call c_getInstOrient(statInOri,statdrIdrG,orient,statRot)
			statInOri(:,4:48) = c_0
			statdrIdrG(:,4:16) = c_0
			call c_getElRu(Ru,dRdU,0,eType,numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt, &
				locNds,globNds,orient,cMat,tExp,ABD,stExp,bStiff,btExp,temp,disp,vel,acc,pDisp,pVel,pAcc, &
				statInOri,statdrIdrG,statRot,den,sMass,bMass)
			i3 = numNds*dofPerNd
			do i4 = 1, i3
			    i5 = dofTable(1,i4)
				i6 = elementList(dofTable(2,i4),i2)
			    i7 = nDofIndex(i5,i6)
				dRudD(i4) = dRudD(i4) + compStepInv*imag(Ru(i4))
			enddo
			if(numIntDof .gt. 0) then
			    i3 = numNds*dofPerNd + 1
				i4 = numNds*dofPerNd + numIntDof
				i5 = intVecRange(i2-1) + 1
				do i6 = i3, i4
				    intdRudD(i5) = intdRudD(i5) + compStepInv*imag(Ru(i6))
					i5 = i5 + 1
				enddo
			endif
		enddo
		
		do i1 = 1, numLds
		    if(time .ge. loadsActTime(1,i1) .and. time .le. loadsActTime(2,i1)) then
			    lt = loadType(i1)
				if(lt .eq. 'bodyForce' .or. lt .eq. 'gravitational' .or. lt .eq. 'centrifugal') then
				    if(lt .eq. 'bodyForce') then
					    isGrav = 0
				    else
					    isGrav = 1
					endif
					do i2 = loadsRange(i1-1) + 1, loadsRange(i1)
						read(loadNodes(i2),*,err=1087) i3
						if(elDepends(i3) .eq. 1) then
							if(lt .eq. 'centrifugal') then
								center(:) = inputLoads(1:3,i2)
								axis(:) = inputLoads(4:6,i2)
								angVel = inputLoads(7,i2)
								call c_getElCentAcc(bdyFrc(1:3),center,axis,angVel,i3)
								bdyFrc(4:6) = c_0
							else
								bdyFrc(:) = inputLoads(1:6,i2)
							endif
							call c_getElBodyForce(ndFrc,dofTable,ndDof,bdyFrc,isGrav,i3)
							do i4 = 1, ndDof
								i5 = dofTable(1,i4)
								i6 = elementList(dofTable(2,i4),i3)
								i7 = nDofIndex(i5,i6)
								dRudD(i7) = dRudD(i7) - compStepInv*imag(ndFrc(i4))
							enddo
						endif
						goto 1110
1087                    do i4 = 1, numElSets
							if(elSetName(i4) .eq. loadNodes(i2)) then
								do i5 = elSetRange(i4-1)+1, elSetRange(i4)
									i6 = elementSets(i5)
									if(elDepends(i6) .eq. 1) then
										if(lt .eq. 'centrifugal') then
											center(:) = inputLoads(1:3,i2)
											axis(:) = inputLoads(4:6,i2)
											angVel = inputLoads(7,i2)
											call c_getElCentAcc(bdyFrc(1:3),center,axis,angVel,i6)
											bdyFrc(4:6) = c_0
										else
											bdyFrc(:) = inputLoads(1:6,i2)
										endif
										call c_getElBodyForce(ndFrc,dofTable,ndDof,bdyFrc,isGrav,i6)
										do i7 = 1, ndDof
											i8 = dofTable(1,i7)
											i9 = elementList(dofTable(2,i7),i6)
											i10 = nDofIndex(i8,i9)
											dRudD(i10) = dRudD(i10) - compStepInv*imag(ndFrc(i7))
										enddo
									endif
								enddo
							endif
						enddo
1110                    i4 = 1
					enddo
                elseif(lt .eq. 'surfacePressure' .or. lt .eq. 'surfaceTraction') then
					do i2 = loadsRange(i1-1) + 1, loadsRange(i1)
					    if(lt .eq. 'surfacePressure') then
						    pressure = inputLoads(5,i2)
							trac(:) = c_0
						else
						    pressure = c_0
							trac(1:3) = inputLoads(5:7,i2)
							trac(4:6) = c_0
						endif
						read(loadNodes(i2),*,err=1137) i3
						if(elDepends(i3) .eq. 1) then
							do i4 = 1, 6
								if(elementSurfaces(i4,i3) .eq. 1) then
									nDir(:) = inputLoads(1:3,i2)
									tTol = inputLoads(4,i2)
									call c_getElSurfaceTraction(ndFrc,dofTable,ndDof,trac,pressure,i3,i4,nDir,tTol)
								endif
							enddo
							do i4 = 1, ndDof
								i5 = dofTable(1,i4)
								i6 = elementList(dofTable(2,i4),i3)
								i7 = nDofIndex(i5,i6)
								dRudD(i7) = dRudD(i7) - compStepInv*imag(ndFrc(i4))
							enddo
						endif
						goto 1157
1137                    do i4 = 1, numElSets
							if(elSetName(i4) .eq. loadNodes(i2)) then
								do i5 = elSetRange(i4-1)+1, elSetRange(i4)
									i6 = elementSets(i5)
									if(elDepends(i6) .eq. 1) then
										do i7 = 1, 6
											if(elementSurfaces(i7,i6) .eq. 1) then
												nDir(:) = inputLoads(1:3,i2)
												tTol = inputLoads(4,i2)
												call c_getElSurfaceTraction(ndFrc,dofTable,ndDof,trac,pressure,i6,i7,nDir,tTol)
											endif
										enddo
										do i7 = 1, ndDof
											i8 = dofTable(1,i7)
											i9 = elementList(dofTable(2,i7),i6)
											i10 = nDofIndex(i8,i9)
											dRudD(i10) = dRudD(i10) - compStepInv*imag(ndFrc(i7))
										enddo
									endif
								enddo
							endif
						enddo
1157                    i4 = 1
					enddo			
				endif
			endif
		enddo
		
		do i1 = dToNdRange(dVarNum-1) + 1, dToNdRange(dVarNum)
		    i2 = dToNd(i1)
		    call c_getNodeDLoad(c_1*time,ndLd,i2)
			do i3 = 1, 6
			    i4 = nDofIndex(i3,i2)
				if(i4 .gt. 0) then
			        dRudD(i4) = dRudD(i4) - compStepInv*imag(ndLd(i3))
				endif
			enddo
		enddo
		
	end subroutine getdRelasticdD
	
	subroutine getTotaldLdD()
	    implicit none
		
		real*8 :: time
		integer :: i1, i2, i3, i4, i5, i6, errFlag
		
		if(dynamic .eq. 0) then
		    dLdt(:) = r_0
			dLdtdot(:) = r_0
			dLdu(:) = r_0
			if(intVecSize .gt. 0) then
			    intdLdu(:) = r_0
			endif
			dLdv(:) = r_0
			dLda(:) = r_0
			call getdObjdU(loadTime)
			call solveDynamicAdjoint()
			dLdD(:) = r_0
			call getdObjdD(loadTime)
			do i1 = 1, numDVar
			    if(solveThermal .eq. 1) then
				    call getdRthermaldD(loadTime,i1)
					do i2 = 1, numNodes
					    dLdD(i1) = dLdD(i1) - tempAdj(i2)*dRtdD(i2)
					enddo
				endif
				if(solveElastic .eq. 1) then
				    call getdRelasticdD(loadTime,i1)
					do i2 = 1, elMatDim
					    dLdD(i1) = dLdD(i1) - dispAdj(i2)*dRudD(i2)
					enddo
					if(intVecSize .gt. 0) then
					    do i2 = 1, intVecSize
						    dLdD(i1) = dLdD(i1) - intDispAdj(i2)*intdRudD(i2)
						enddo
					endif
				endif
			enddo
		else
		    dLdD(:) = r_0
		    call readBinarySolution(numTSteps,errFlag)
			time = delT*numTSteps
			do i1 = numTSteps, 1, -1
			    dLdt(:) = r_0
				dLdtdot(:) = r_0
				dLdu(:) = r_0
				if(intVecSize .gt. 0) then
					intdLdu(:) = r_0
				endif
				dLdv(:) = r_0
				dLda(:) = r_0
				if(i1 .lt. numTSteps) then
				    call augmentdLdu()
				endif
				nodeTemp(:) = prevTemp(:)
				nodeTdot(:) = prevTdot(:)
				nodeDisp(:) = prevDisp(:)
				nodeVel(:) = prevVel(:)
				nodeAcc(:) = prevAcc(:)
				if(intVecSize .gt. 0) then
					internalDisp(:) = prevIntDisp(:)
				endif
				call readBinarySolution(i1-1,errFlag)
				call getdObjdU(time)
				call solveDynamicAdjoint()
				call getdObjdD(time)
				do i2 = 1, numDVar
					if(solveThermal .eq. 1) then
						call getdRthermaldD(time,i2)
						do i3 = 1, numNodes
							dLdD(i2) = dLdD(i2) - tempAdj(i3)*dRtdD(i3)
						enddo
					endif
					if(solveElastic .eq. 1) then
						call getdRelasticdD(time,i2)
						do i3 = 1, elMatDim
							dLdD(i2) = dLdD(i2) - dispAdj(i3)*dRudD(i3)
						enddo
						if(intVecSize .gt. 0) then
							do i3 = 1, intVecSize
								dLdD(i2) = dLdD(i2) - intDispAdj(i3)*intdRudD(i3)
							enddo
						endif
					endif
				enddo
				time = time - delT
			enddo
		endif
		
	end subroutine getTotaldLdD
	
	subroutine processWriteNodeEl(jobUnit,fileLine,iosVal)
	    implicit none
		
		integer, intent(in) :: jobUnit
		character(len=256), intent(out) :: fileLine
		integer, intent(out) :: iosVal
		
		integer, allocatable :: allNdSet(:), allElSet(:), writeTSteps(:)
		
		character(len=32) :: commandTag
		character(len=128) :: outputFileName, fullFileName
		
		character(len=64) :: writeFields(10)
	    character(len=16) :: extension, stepStr
		real*8 :: time
		
		integer :: nSInd, eSInd, numFields, readInt(3)
		integer :: i1, i2, i3, i4, i5, i6, i7, errFlag
		
		i1 = index(fileLine,'*')
		commandTag = fileLine(i1:i1+31)
		read(jobUnit,'(A)',iostat=iosVal) fileLine(16:256)
		i1 = index(fileLine,'*')
		numFields = 0
		nSInd = 0
		eSInd = 0
		do while(i1 .eq. 0 .and. iosVal .eq. 0)
			i2 = index(fileLine,':')
			if(i2 .gt. 0) then
				if(fileLine(i2-8:i2) .eq. 'fileName:') then
					read(fileLine(i2+1:i2+128),*) outputFileName
					read(jobUnit,'(A)',iostat=iosVal) fileLine(16:256)
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
					read(jobUnit,'(A)',iostat=iosVal) fileLine(16:256)
					i1 = index(fileLine,'*')
				elseif(fileLine(i2-10:i2) .eq. 'elementSet:') then
				    do i3 = 1, numElSets
						i4 = index(fileLine,elSetName(i3))
						if(i4 .gt. 0) then
							eSInd = i3
						endif
					enddo
					if(eSInd .eq. 0) then
						write(lfUnit,*) 'Warning: no element set named ', fileLine(i2+1:i2+64), 'was found'
						write(lfUnit,*) 'Writing results for all elements'
					endif
					read(jobUnit,'(A)',iostat=iosVal) fileLine(16:256)
					i1 = index(fileLine,'*')
				elseif(fileLine(i2-6:i2) .eq. 'fields:') then
					read(jobUnit,'(A)',iostat=iosVal) fileLine(16:256)
					i2 = index(fileLine,':')
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
						read(jobUnit,'(A)',iostat=iosVal) fileLine(16:256)
						i2 = index(fileLine,':')
					enddo
					i1 = index(fileLine,'*')
				elseif(fileLine(i2-9:i2) .eq. 'timeSteps:') then
					if(numTSteps .le. 0) then
						write(lfUnit,*) 'Error: A dynamic analysis must be run prior to writing time step dependent results.'
						write(lfUnit,*) 'Aborting ', commandTag
						goto 251
					endif
					if(.not. allocated(writeTSteps)) then
						allocate(writeTSteps(numTSteps))
					endif
					writeTSteps(:) = 0
					i3 = index(fileLine,'all')
					if(i3 .gt. 0) then
						writeTSteps(:) = 1
						read(jobUnit,'(A)',iostat=iosVal) fileLine(16:256)
						i1 = index(fileLine,'*')
					else
						read(jobUnit,'(A)',iostat=iosVal) fileLine(16:256)
						i2 = index(fileLine,':')
						do while(i2 .eq. 0 .and. iosVal .eq. 0)
							i3 = index(fileLine, '-')
							if(i3 .gt. 0) then
								i4 = index(fileLine,'[')
								if(i4 .gt. 0) then
									i5 = index(fileLine,']')
									readInt(3) = 1
									read(fileLine(i4+1:i5-1),*,end=254) readInt(1:3)
254                                 do i6 = readInt(1), readInt(2), readInt(3)
										if(i6 .le. numTSteps .and. i6 .gt. 0) then
											writeTSteps(i6) = 1
										else
											write(lfUnit,*) 'Error: a time step specified for writeNodeResults is out of'
											write(lfUnit,*) 'the simulated range.  Aborting ', commandTag
											goto 251
										endif
									enddo
								else
									read(fileLine(i3+1:i3+64),*) i5
									if(i5 .le. numTSteps .and. i5 .gt. 0) then
										writeTSteps(i5) = 1
									else
										write(lfUnit,*) 'Error: a time step specified for writeNodeResults is out of'
										write(lfUnit,*) 'the simulated range.  Aborting ', commandTag
										goto 251
									endif
								endif
							endif
							read(jobUnit,'(A)',iostat=iosVal) fileLine(16:256)
							i1 = index(fileLine,'*')
							i2 = index(fileLine,':')
						enddo
					endif
				else
					read(jobUnit,'(A)',iostat=iosVal) fileLine(16:256)
					i1 = index(fileLine,'*')
				endif
			endif
		enddo
		if(.not. allocated(writeTSteps)) then
		    if(commandTag(1:17) .eq. '*writeNodeResults') then
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
				if(eSInd .eq. 0) then
					if(.not. allocated(allElSet)) then
						allocate(allElSet(numEls))
					endif
					do i3 = 1, numEls
						allElSet(i3) = i3
					enddo
					time = numTSteps*delT
					call writeElementResults(outputFileName,writeFields,numFields,allElSet,numEls,time,numTSteps)
				else
					i3 = elSetRange(eSInd-1) + 1
					i4 = elSetRange(eSInd)
					i5 = i4 - i3 + 1
					time = numTSteps*delT
					call writeElementResults(outputFileName,writeFields,numFields,elementSets(i3:i4),i5,time,numTSteps)
				endif
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
						if(commandTag(1:17) .eq. '*writeNodeResults') then
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
							if(eSInd .eq. 0) then
								if(.not. allocated(allElSet)) then
									allocate(allElSet(numEls))
								endif
								do i3 = 1, numEls
									allElSet(i3) = i3
								enddo
								time = numTSteps*delT
								call writeElementResults(fullFileName,writeFields,numFields,allElSet,numEls,time,numTSteps)
							else
								i3 = elSetRange(eSInd-1) + 1
								i4 = elSetRange(eSInd)
								i5 = i4 - i3 + 1
								time = numTSteps*delT
								call writeElementResults(fullFileName,writeFields,numFields,elementSets(i3:i4),i5,time,numTSteps)
							endif
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
		if(allocated(allElSet)) then
		    deallocate(allElSet)
		endif
251		i1 = index(fileLine,'*')
		
	end subroutine processWriteNodeEl
	
	subroutine processWriteNdElProps(jobUnit,fileLine,iosVal)
	    implicit none
		
		integer, intent(in) :: jobUnit
		character(len=256), intent(out) :: fileLine
		integer, intent(out) :: iosVal
		
		integer, allocatable :: allNdSet(:), allElSet(:)
		character(len=16) :: propList(32)
		character(len=128) :: fileName
		character(len=64) :: setName
		integer :: i1, i2, i3, i4, i5, nSInd, eSInd, numProps, wtNd
		
		i1 = index(fileLine,'*')
		if(fileLine(i1:i1+20) .eq. '*writeNodeCoordinates') then
		    wtNd = 1
		else
		    wtNd = 0
		endif
		
		numProps = 0
		nSInd = 0
		eSInd = 0
		
		read(jobUnit,'(A)',iostat=iosVal) fileLine(16:256)
		i1 = index(fileLine,'*')
		do while(i1 .eq. 0 .and. iosVal .eq. 0)
		    i2 = index(fileLine,':')
			if(i2 .gt. 0) then
			    if(fileLine(i2-8:i2) .eq. 'fileName:') then
				    read(fileLine(i2+1:i2+128),*) fileName
					read(jobUnit,'(A)',iostat=iosVal) fileLine(16:256)
		            i1 = index(fileLine,'*')
				elseif(fileLine(i2-7:i2) .eq. 'nodeSet:') then
				    read(fileLine(i2+1:i2+64),*) setName
					do i3 = 1, numNdSets
					    if(setName .eq. ndSetName(i3)) then
						    nSInd = i3
						endif
					enddo
					if(nSInd .eq. 0) then
					    write(lfUnit,*) 'Warning: there is no node set with the name ', setName
						write(lfUnit,*) 'Writing coordinates of all nodes.'
					endif
					read(jobUnit,'(A)',iostat=iosVal) fileLine(16:256)
		            i1 = index(fileLine,'*')
				elseif(fileLine(i2-10:i2) .eq. 'elementSet:') then
				    read(fileLine(i2+1:i2+64),*) setName
					do i3 = 1, numElSets
					    if(setName .eq. elSetName(i3)) then
						    eSInd = i3
						endif
					enddo
					if(eSInd .eq. 0) then
					    write(lfUnit,*) 'Warning: there is no element set with the name ', setName
						write(lfUnit,*) 'Writing properties of all elements.'
					endif
					read(jobUnit,'(A)',iostat=iosVal) fileLine(16:256)
		            i1 = index(fileLine,'*')
				elseif(fileLine(i2-10:i2) .eq. 'properties:') then
				    read(jobUnit,'(A)',iostat=iosVal) fileLine(16:256)
		            i1 = index(fileLine,'*')
					i2 = index(fileLine,':')
					do while(i2 .eq. 0 .and. iosVal .eq. 0)
					    i3 = index(fileLine,'-')
						if(i3 .gt. 0) then
						    numProps = numProps + 1
							read(fileLine(i3+1:i3+16),*) propList(numProps)
						endif
						read(jobUnit,'(A)',iostat=iosVal) fileLine(16:256)
						i2 = index(fileLine,':')
					enddo
					i1 = index(fileLine,'*')
				else
				    read(jobUnit,'(A)',iostat=iosVal) fileLine(16:256)
		            i1 = index(fileLine,'*')
				endif
			else
			    read(jobUnit,'(A)',iostat=iosVal) fileLine(16:256)
		        i1 = index(fileLine,'*')
			endif
		enddo
		if(wtNd .eq. 1) then
		    if(nSInd .eq. 0) then
			    allocate(allNdSet(numNodes))
				do i3 = 1, numNodes
				    allNdSet(i3) = i3
				enddo
				call writeNodeCoord(fileName,allNdSet,numNodes)
				deallocate(allNdSet)
			else
			    i3 = ndSetRange(nSInd-1) + 1
				i4 = ndSetRange(nSInd)
				i5 = i4 - i3 + 1
				call writeNodeCoord(fileName,nodeSets(i3:i4),i5)
			endif
		else
		    if(eSInd .eq. 0) then
			    allocate(allElSet(numEls))
				do i3 = 1, numEls
				    allElSet(i3) = i3
				enddo
				call writeElProperties(fileName,allElSet,numEls,propList,numProps)
				deallocate(allElSet)
			else
			    i3 = elSetRange(eSInd-1) + 1
				i4 = elSetRange(eSInd)
				i5 = i4 - i3 + 1
				call writeElProperties(fileName,elementSets(i3:i4),i5,propList,numProps)
			endif
		endif
		
	end subroutine processWriteNdElProps
	
end module AStrO_commandFunctions