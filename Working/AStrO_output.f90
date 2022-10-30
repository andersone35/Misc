module AStrO_output
    use AStrO_globalData
	use AStrO_constantVals
	use AStrO_r_elementEqns
	
	contains
	
!! General output
	
	subroutine writeNodeResults(fileName,fields,numFields,ndSet,nSLen,time,tStep)
	    implicit none
		
		integer, intent(in) :: numFields, nSlen, tStep
		real*8, intent(in) :: time
		character(len=128), intent(in) :: fileName
		character(len=16), intent(in) :: fields(numFields)
		integer, intent(in) :: ndSet(nSLen)
		
		character(len=16) :: intStr
		character(len=128) :: fltStr
		integer :: i1, i2, i3, i4, i5, i6
		
		open(unit=15, file=fileName, status='replace', action='write')
		
		write(15,*) 'time: ', time
		write(15,*) 'timeStep: ', tStep
		
		do i1 = 1, numFields
		    if(fields(i1) .eq. 'temperature') then
			    write(15,*) 'temperature:'
			    do i2 = 1, nSLen
				    i3 = ndSet(i2)
					i4 = currentRank(i3)
					write(intStr,'(I0)') i3
					write(fltStr,'(G18.12)') nodeTemp(i4)
					write(15,*) '    - [', intStr, ', ', fltStr, ']'
				enddo
			elseif(fields(i1) .eq. 'tdot') then
			    write(15,*) 'tdot:'
			    do i2 = 1, nSLen
				    i3 = ndSet(i2)
					i4 = currentRank(i3)
					write(intStr,'(I0)') i3
					write(fltStr,'(G18.12)') nodeTdot(i4)
					write(15,*) '    - [', intStr, ', ', fltStr, ']'
				enddo
			elseif(fields(i1) .eq. 'reactionHeatGen') then
			    write(15,*) 'reactionHeatGen:'
			    do i2 = 1, nSLen
				    i3 = ndSet(i2)
					i4 = currentRank(i3)
					write(intStr,'(I0)') i3
					write(fltStr,'(G18.12)') -thermalLoad(i4)
					write(15,*) '    - [', intStr, ', ', fltStr, ']'
				enddo
			elseif(fields(i1) .eq. 'displacement') then
			    write(15,*) 'displacement:'
			    do i2 = 1, nSLen
				    i3 = ndSet(i2)
					write(intStr,'(I0)') i3
					i6 = 1
					do i4 = 1, 6
					    i5 = nDofIndex(i4,i3)
						if(i5 .gt. 0) then
						    write(fltStr(i6:i6+17),'(G18.12)') nodeDisp(i5)
							i6 = i6 + 18
							fltStr(i6:i6) = ','
							i6 = i6 + 1
						else
						    write(fltStr(i6:i6+17),'(G18.12)') r_0
							i6 = i6 + 18
							fltStr(i6:i6) = ','
							i6 = i6 + 1
						endif
					enddo
					write(15,*) '    - [', intStr, ', ', fltStr(1:i6-2), ']'
				enddo
			elseif(fields(i1) .eq. 'velocity') then
			    write(15,*) 'velocity:'
			    do i2 = 1, nSLen
				    i3 = ndSet(i2)
					write(intStr,'(I0)') i3
					i6 = 1
					do i4 = 1, 6
					    i5 = nDofIndex(i4,i3)
						if(i5 .gt. 0) then
						    write(fltStr(i6:i6+17),'(G18.12)') nodeVel(i5)
							i6 = i6 + 18
							fltStr(i6:i6) = ','
							i6 = i6 + 1
						else
						    write(fltStr(i6:i6+17),'(G18.12)') r_0
							i6 = i6 + 18
							fltStr(i6:i6) = ','
							i6 = i6 + 1
						endif
					enddo
					write(15,*) '    - [', intStr, ', ', fltStr(1:i6-2), ']'
				enddo
			elseif(fields(i1) .eq. 'acceleration') then
			    write(15,*) 'acceleration:'
			    do i2 = 1, nSLen
				    i3 = ndSet(i2)
					write(intStr,'(I0)') i3
					i6 = 1
					do i4 = 1, 6
					    i5 = nDofIndex(i4,i3)
						if(i5 .gt. 0) then
						    write(fltStr(i6:i6+17),'(G18.12)') nodeAcc(i5)
							i6 = i6 + 18
							fltStr(i6:i6) = ','
							i6 = i6 + 1
						else
						    write(fltStr(i6:i6+17),'(G18.12)') r_0
							i6 = i6 + 18
							fltStr(i6:i6) = ','
							i6 = i6 + 1
						endif
					enddo
					write(15,*) '    - [', intStr, ', ', fltStr(1:i6-2), ']'
				enddo
			elseif(fields(i1) .eq. 'reactionForce') then
			    write(15,*) 'reactionForce:'
			    do i2 = 1, nSLen
				    i3 = ndSet(i2)
					write(intStr,'(I0)') i3
					i6 = 1
					do i4 = 1, 6
					    i5 = nDofIndex(i4,i3)
						if(i5 .gt. 0) then
						    write(fltStr(i6:i6+17),'(G18.12)') -elasticLoad(i5)
							i6 = i6 + 18
							fltStr(i6:i6) = ','
							i6 = i6 + 1
						else
						    write(fltStr(i6:i6+17),'(G18.12)') r_0
							i6 = i6 + 18
							fltStr(i6:i6) = ','
							i6 = i6 + 1
						endif
					enddo
					write(15,*) '    - [', intStr, ', ', fltStr(1:i6-2), ']'
				enddo
			else
			    write(lfUnit,*) 'Warning: Unrecognized field name ', fields(i1), ' listed in writeNodeResults'
				write(lfUnit,*) ' '
			endif
		enddo
		
		close(15)
		
	end subroutine writeNodeResults
	
	subroutine writeElementResults(fileName,fields,numFields,elSet,eSLen,time,tStep)
	    implicit none
		
		integer, intent(in) :: numFields, eSLen, tStep
		real*8, intent(in) :: time
		character(len=128), intent(in) :: fileName
		character(len=16), intent(in) :: fields(numFields)
		integer, intent(in) :: elSet(eSLen)
		
		integer :: numNds, dofPerNd, numIntDof, numIntPts
		integer :: dofTable(2,33)
		real*8 :: intPts(3,8), ipWt(8)
		
		real*8 :: flux(3), dfluxdT(3,10), Tx(3), dTxdT(3,10)
		
		real*8 :: frcMom(9), shDef(9), ptT
		
		character(len=16) :: elNumStr, secNumStr, ipStr, layNumStr
        character(len=32) :: intStr
		character(len=256) :: fltStr
		real*8 :: stress(6), strain(6), sEDen, strength, fIndex, maxFInd
		integer :: i1, i2, i3, i4, i5, i6, i7, eType, secNum, matID, layerInd, numLayers
		
		open(unit=138, file=fileName, status='replace', action='write')
		
		do i1 = 1, numFields
		    if(fields(i1) .eq. 'stress') then
			    write(138,*) 'stress:'
				write(138,*) '## element, section, integration point, layer, S11, S22, S33, S12, S13, S23, max index '
				do i2 = 1, eSLen
				    i3 = elSet(i2)
					eType = elementType(i3)
					secNum = elementSection(i3)
					numLayers = secLayupRange(secNum) - secLayupRange(secNum-1)
					call r_getElementProfile(numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt,i3)
					if(eType .eq. 41 .or. eType .eq. 3) then
					    do i4 = 1, numIntPts
						    do i5 = 1, numLayers
							    call r_getElStress(stress,strain,i3,i5,intPts(:,i4))
								write(elNumStr,'(I0)') i3
								write(secNumStr,'(I0)') secNum
								write(ipStr,'(I0)') i4
								write(layNumStr,'(I0)') i5
								intStr = elNumStr // ', ' // secNumStr // ', ' // ipStr // ', ' // layNumStr
                                i6 = 1
                                do i7 = 1, 6
								    write(fltStr(i6:i6+17),'(G18.12)') stress(i7)
									i6 = i6 + 18
									fltStr(i6:i6) = ','
									i6 = i6 + 1
                                enddo
								layerInd = secLayupRange(secNum-1) + i5
								matID = layupMatId(layerInd)
								maxFInd = r_0
								do i7 = 1, 3
								    strength = materialMaxStress(i7,matID)
								    if(strength .ne. r_0) then
									    fIndex = stress(i7)/strength
										if(fIndex .gt. maxFInd) then
										    maxFInd = fIndex
										endif
									endif
									strength = materialMaxStress(i7+3,matID)
								    if(strength .ne. r_0) then
									    fIndex = stress(i7)/strength
										if(fIndex .gt. maxFInd) then
										    maxFInd = fIndex
										endif
									endif
									strength = materialMaxStress(i7+6,matID)
								    if(strength .ne. r_0) then
									    fIndex = abs(stress(i7+3))/strength
										if(fIndex .gt. maxFInd) then
										    maxFInd = fIndex
										endif
									endif
								enddo
								if(maxFInd .eq. r_0) then
								    fltStr(i6:i6+3) = ' NA,'
									i6 = i6 + 4
								else
								    write(fltStr(i6:i6+17),'(G18.12)') maxFInd
									i6 = i6 + 18
									fltStr(i6:i6) = ','
									i6 = i6 + 1
								endif
								write(138,*) '    - [', intStr, ', ', fltStr(1:i6-2), ']'
							enddo
						enddo
					elseif(eType .eq. 2) then
					else
					    do i4 = 1, numIntPts
							call r_getElStress(stress,strain,i3,0,intPts(:,i4))
							write(elNumStr,'(I0)') i3
							write(secNumStr,'(I0)') secNum
							write(ipStr,'(I0)') i4
							layNumStr = ' NA'
							intStr = elNumStr // ', ' // secNumStr// ', ' // ipStr // ', ' // layNumStr
							i6 = 1
							do i7 = 1, 6
								write(fltStr(i6:i6+17),'(G18.12)') stress(i7)
								i6 = i6 + 18
								fltStr(i6:i6) = ','
								i6 = i6 + 1
							enddo
							matID = sectionMatId(secNum)
							maxFInd = r_0
							do i7 = 1, 3
								strength = materialMaxStress(i7,matID)
								if(strength .ne. r_0) then
									fIndex = stress(i7)/strength
									if(fIndex .gt. maxFInd) then
										maxFInd = fIndex
									endif
								endif
								strength = materialMaxStress(i7+3,matID)
								if(strength .ne. r_0) then
									fIndex = stress(i7)/strength
									if(fIndex .gt. maxFInd) then
										maxFInd = fIndex
									endif
								endif
								strength = materialMaxStress(i7+6,matID)
								if(strength .ne. r_0) then
									fIndex = abs(stress(i7+3))/strength
									if(fIndex .gt. maxFInd) then
										maxFInd = fIndex
									endif
								endif
							enddo
							if(maxFInd .eq. r_0) then
								fltStr(i6:i6+3) = ' NA,'
								i6 = i6 + 4
							else
								write(fltStr(i6:i6+17),'(G18.12)') maxFInd
								i6 = i6 + 18
								fltStr(i6:i6) = ','
								i6 = i6 + 1
							endif
							write(138,*) '    - [', intStr, ', ', fltStr(1:i6-2), ']'
						enddo
					endif
				enddo
			elseif(fields(i1) .eq. 'strain') then
			    write(138,*) 'strain:'
				write(138,*) '## element, section, integration point, layer, E11, E22, E33, E12, E13, E23, max index '
				do i2 = 1, eSLen
				    i3 = elSet(i2)
					eType = elementType(i3)
					secNum = elementSection(i3)
					numLayers = secLayupRange(secNum) - secLayupRange(secNum-1)
					call r_getElementProfile(numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt,i3)
					if(eType .eq. 41 .or. eType .eq. 3) then
					    do i4 = 1, numIntPts
						    do i5 = 1, numLayers
							    call r_getElStress(stress,strain,i3,i5,intPts(:,i4))
								write(elNumStr,'(I0)') i3
								write(secNumStr,'(I0)') secNum
								write(ipStr,'(I0)') i4
								write(layNumStr,'(I0)') i5
								intStr = elNumStr // ', ' // secNumStr // ', ' // ipStr // ', ' // layNumStr
                                i6 = 1
                                do i7 = 1, 6
								    write(fltStr(i6:i6+17),'(G18.12)') strain(i7)
									i6 = i6 + 18
									fltStr(i6:i6) = ','
									i6 = i6 + 1
                                enddo
								layerInd = secLayupRange(secNum-1) + i5
								matID = layupMatId(layerInd)
								maxFInd = r_0
								do i7 = 1, 3
								    strength = materialMaxStrain(i7,matID)
								    if(strength .ne. r_0) then
									    fIndex = strain(i7)/strength
										if(fIndex .gt. maxFInd) then
										    maxFInd = fIndex
										endif
									endif
									strength = materialMaxStrain(i7+3,matID)
								    if(strength .ne. r_0) then
									    fIndex = strain(i7)/strength
										if(fIndex .gt. maxFInd) then
										    maxFInd = fIndex
										endif
									endif
									strength = materialMaxStrain(i7+6,matID)
								    if(strength .ne. r_0) then
									    fIndex = abs(strain(i7+3))/strength
										if(fIndex .gt. maxFInd) then
										    maxFInd = fIndex
										endif
									endif
								enddo
								if(maxFInd .eq. r_0) then
								    fltStr(i6:i6+3) = ' NA,'
									i6 = i6 + 4
								else
								    write(fltStr(i6:i6+17),'(G18.12)') maxFInd
									i6 = i6 + 18
									fltStr(i6:i6) = ','
									i6 = i6 + 1
								endif
								write(138,*) '    - [', intStr, ', ', fltStr(1:i6-2), ']'
							enddo
						enddo
					elseif(eType .eq. 2) then
					else
					    do i4 = 1, numIntPts
							call r_getElStress(stress,strain,i3,0,intPts(:,i4))
							write(elNumStr,'(I0)') i3
							write(secNumStr,'(I0)') secNum
							write(ipStr,'(I0)') i4
							layNumStr = ' NA'
							intStr = elNumStr // ', ' // secNumStr// ', ' // ipStr // ', ' // layNumStr
							i6 = 1
							do i7 = 1, 6
								write(fltStr(i6:i6+17),'(G18.12)') strain(i7)
								i6 = i6 + 18
								fltStr(i6:i6) = ','
								i6 = i6 + 1
							enddo
							matID = sectionMatId(secNum)
							maxFInd = r_0
							do i7 = 1, 3
								strength = materialMaxStrain(i7,matID)
								if(strength .ne. r_0) then
									fIndex = strain(i7)/strength
									if(fIndex .gt. maxFInd) then
										maxFInd = fIndex
									endif
								endif
								strength = materialMaxStrain(i7+3,matID)
								if(strength .ne. r_0) then
									fIndex = strain(i7)/strength
									if(fIndex .gt. maxFInd) then
										maxFInd = fIndex
									endif
								endif
								strength = materialMaxStrain(i7+6,matID)
								if(strength .ne. r_0) then
									fIndex = abs(strain(i7+3))/strength
									if(fIndex .gt. maxFInd) then
										maxFInd = fIndex
									endif
								endif
							enddo
							if(maxFInd .eq. r_0) then
								fltStr(i6:i6+3) = ' NA,'
								i6 = i6 + 4
							else
								write(fltStr(i6:i6+17),'(G18.12)') maxFInd
								i6 = i6 + 18
								fltStr(i6:i6) = ','
								i6 = i6 + 1
							endif
							write(138,*) '    - [', intStr, ', ', fltStr(1:i6-2), ']'
						enddo
					endif
				enddo
            elseif(fields(i1) .eq. 'strainEnergyDen') then
			    write(138,*) 'strainEnergyDen:'
				write(138,*) '## element, section, integration point, layer, strain energy density, max Index '
				do i2 = 1, eSLen
				    i3 = elSet(i2)
					eType = elementType(i3)
					secNum = elementSection(i3)
					numLayers = secLayupRange(secNum) - secLayupRange(secNum-1)
					call r_getElementProfile(numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt,i3)
					if(eType .eq. 41 .or. eType .eq. 3) then
					    do i4 = 1, numIntPts
						    do i5 = 1, numLayers
							    call r_getElStress(stress,strain,i3,i5,intPts(:,i4))
								write(elNumStr,'(I0)') i3
								write(secNumStr,'(I0)') secNum
								write(ipStr,'(I0)') i4
								write(layNumStr,'(I0)') i5
								intStr = elNumStr // ', ' // secNumStr // ', ' // ipStr // ', ' // layNumStr
                                sEDen = r_0
                                do i7 = 1, 6
								    sEDen = sEDen + r_p5*stress(i7)*strain(i7)
                                enddo
								write(fltStr(1:18),'(G18.12)') sEDen
								fltStr(19:19) = ','
								layerInd = secLayupRange(secNum-1) + i5
								matID = layupMatId(layerInd)
								strength = materialMaxStrnEngy(matID)
								if(strength .ne. r_0) then
								    maxFInd = sEDen/strength
									write(fltStr(20:37),'(G18.12)') maxFInd
								else
								    fltStr(20:37) = ' NA               '
								endif
								write(138,*) '    - [', intStr, ', ', fltStr(1:37), ']'
							enddo
						enddo
					elseif(eType .eq. 2) then
					else
					    do i4 = 1, numIntPts
							call r_getElStress(stress,strain,i3,0,intPts(:,i4))
							write(elNumStr,'(I0)') i3
							write(secNumStr,'(I0)') secNum
							write(ipStr,'(I0)') i4
							layNumStr = ' NA'
							intStr = elNumStr // ', ' // secNumStr // ', ' // ipStr // ', ' // layNumStr
							sEDen = r_0
							do i7 = 1, 6
								sEDen = sEDen + r_p5*stress(i7)*strain(i7)
							enddo
							write(fltStr(1:18),'(G18.12)') sEDen
							fltStr(19:19) = ','
							matID = sectionMatId(secNum)
							strength = materialMaxStrnEngy(matID)
							if(strength .ne. r_0) then
								maxFInd = sEDen/strength
								write(fltStr(20:37),'(G18.12)') maxFInd
							else
								fltStr(20:37) = ' NA               '
							endif							
						    write(138,*) '    - [', intStr, ', ', fltStr(1:37), ']'
						enddo
					endif
				enddo
			elseif(fields(i1) .eq. 'heatFlux') then
			    write(138,*) 'heatFlux:'
				write(138,*) '## element, section, integration point, q1, q2, q3 '
				do i2 = 1, eSLen
				    i3 = elSet(i2)
					eType = elementType(i3)
					secNum = elementSection(i3)
					call r_getElementProfile(numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt,i3)
					do i4 = 1, numIntPts
						call r_getElFlux(flux,dfluxdT,Tx,dTxdT,i3,intPts(:,i4))
						write(elNumStr,'(I0)') i3
						write(secNumStr,'(I0)') secNum
						write(ipStr,'(I0)') i4
						intStr = elNumStr // ', ' // secNumStr // ', ' // ipStr
						i6 = 1
						do i7 = 1, 3
							write(fltStr(i6:i6+17),'(G18.12)') flux(i7)
							i6 = i6 + 18
							fltStr(i6:i6) = ','
							i6 = i6 + 1
						enddo
						write(138,*) '    - [', intStr, ', ', fltStr(1:i6-2), ']'
					enddo
				enddo
			elseif(fields(i1) .eq. 'tempGradient') then
			    write(138,*) 'tempGradient:'
				write(138,*) '## element, section, integration point, Tg1, Tg2, Tg3 '
				do i2 = 1, eSLen
				    i3 = elSet(i2)
					eType = elementType(i3)
					secNum = elementSection(i3)
					call r_getElementProfile(numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt,i3)
					do i4 = 1, numIntPts
						call r_getElFlux(flux,dfluxdT,Tx,dTxdT,i3,intPts(:,i4))
						write(elNumStr,'(I0)') i3
						write(secNumStr,'(I0)') secNum
						write(ipStr,'(I0)') i4
						intStr = elNumStr // ', ' // secNumStr // ', ' // ipStr
						i6 = 1
						do i7 = 1, 3
							write(fltStr(i6:i6+17),'(G18.12)') Tx(i7)
							i6 = i6 + 18
							fltStr(i6:i6) = ','
							i6 = i6 + 1
						enddo
						write(138,*) '    - [', intStr, ', ', fltStr(1:i6-2), ']'
					enddo
				enddo
			elseif(fields(i1) .eq. 'shellDeformation') then
			    write(138,*) 'shellDeformation:'
				write(138,*) '## element, section, integration point, E11, E22, E12, K11, K22, K12 '
				do i2 = 1, eSLen
				    i3 = elSet(i2)
					eType = elementType(i3)
					if(eType .eq. 41 .or. eType .eq. 3) then
						secNum = elementSection(i3)
						call r_getElementProfile(numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt,i3)
						do i4 = 1, numIntPts
							call r_getShellFrcMom(frcMom,shDef,ptT,i3,intPts(:,i4))
							write(elNumStr,'(I0)') i3
							write(secNumStr,'(I0)') secNum
							write(ipStr,'(I0)') i4
							intStr = elNumStr // ', ' // secNumStr // ', ' // ipStr
							i6 = 1
							do i7 = 1, 6
								write(fltStr(i6:i6+17),'(G18.12)') shDef(i7)
								i6 = i6 + 18
								fltStr(i6:i6) = ','
								i6 = i6 + 1
							enddo
							write(138,*) '    - [', intStr, ', ', fltStr(1:i6-2), ']'
						enddo
					endif
				enddo
			elseif(fields(i1) .eq. 'shellForceMoment') then
			    write(138,*) 'shellForceMoment:'
				write(138,*) '## element, section, integration point, F11, F22, F12, M11, M22, M12 '
				do i2 = 1, eSLen
				    i3 = elSet(i2)
					eType = elementType(i3)
					if(eType .eq. 41 .or. eType .eq. 3) then
						secNum = elementSection(i3)
						call r_getElementProfile(numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt,i3)
						do i4 = 1, numIntPts
							call r_getShellFrcMom(frcMom,shDef,ptT,i3,intPts(:,i4))
							write(elNumStr,'(I0)') i3
							write(secNumStr,'(I0)') secNum
							write(ipStr,'(I0)') i4
							intStr = elNumStr // ', ' // secNumStr // ', ' // ipStr
							i6 = 1
							do i7 = 1, 6
								write(fltStr(i6:i6+17),'(G18.12)') frcMom(i7)
								i6 = i6 + 18
								fltStr(i6:i6) = ','
								i6 = i6 + 1
							enddo
							write(138,*) '    - [', intStr, ', ', fltStr(1:i6-2), ']'
						enddo
					endif
				enddo
			elseif(fields(i1) .eq. 'beamDeformation') then
			    write(138,*) 'beamDeformation:'
				write(138,*) '## element, section, integration point, E11, E21, E31, K11, K21, K31 '
				do i2 = 1, eSLen
				    i3 = elSet(i2)
					eType = elementType(i3)
					if(eType .eq. 2) then
						secNum = elementSection(i3)
						call r_getElementProfile(numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt,i3)
						do i4 = 1, numIntPts
							call r_getBeamFrcMom(frcMom,shDef,i3,intPts(:,i4))
							write(elNumStr,'(I0)') i3
							write(secNumStr,'(I0)') secNum
							write(ipStr,'(I0)') i4
							intStr = elNumStr // ', ' // secNumStr // ', ' // ipStr
							i6 = 1
							do i7 = 1, 6
								write(fltStr(i6:i6+17),'(G18.12)') shDef(i7)
								i6 = i6 + 18
								fltStr(i6:i6) = ','
								i6 = i6 + 1
							enddo
							write(138,*) '    - [', intStr, ', ', fltStr(1:i6-2), ']'
						enddo
					endif
				enddo
			elseif(fields(i1) .eq. 'beamForceMoment') then
			    write(138,*) 'beamForceMoment:'
				write(138,*) '## element, section, integration point, F11, F21, F31, M11, M21, M31 '
				do i2 = 1, eSLen
				    i3 = elSet(i2)
					eType = elementType(i3)
					if(eType .eq. 2) then
						secNum = elementSection(i3)
						call r_getElementProfile(numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt,i3)
						do i4 = 1, numIntPts
							call r_getBeamFrcMom(frcMom,shDef,i3,intPts(:,i4))
							write(elNumStr,'(I0)') i3
							write(secNumStr,'(I0)') secNum
							write(ipStr,'(I0)') i4
							intStr = elNumStr // ', ' // secNumStr // ', ' // ipStr
							i6 = 1
							do i7 = 1, 6
								write(fltStr(i6:i6+17),'(G18.12)') frcMom(i7)
								i6 = i6 + 18
								fltStr(i6:i6) = ','
								i6 = i6 + 1
							enddo
							write(138,*) '    - [', intStr, ', ', fltStr(1:i6-2), ']'
						enddo
					endif
				enddo
			else
			    write(lfUnit,*) 'Warning: Unrecognized field name ', fields(i1), ' listed in writeNodeResults'
				write(lfUnit,*) ' '
			endif
		enddo
		
		close(138)
		
	end subroutine writeElementResults
	
	subroutine writeDVarValues(fileName,fields,numFields)
	    implicit none
		
		character(len=128), intent(in) :: fileName
		integer, intent(in) :: numFields
		character(len=16), intent(in) :: fields(numFields)
		
		character(len=128) :: writeLine
		integer :: i1, i2, i3
		
		open(unit=610, file=fileName, status='replace', action='write')
		
		write(610,*) 'designVariableValues:'
		
		do i1 = 1, numDVar
		    writeLine(1:13) = '    - value: '
			write(writeLine(14:31),'(G18.12)') r_dVec(i1)
			write(610,*) writeLine(1:31)
			do i2 = 1, numFields
			    if(fields(i2) .eq. 'category') then
					write(610,*) '      category: ', dCategory(i1)
				elseif(fields(i2) .eq. 'subCategory') then
					write(610,*) '      subCategory: ', dCategory(i1)
				elseif(fields(i2) .eq. 'component') then
					write(610,*) '      component: ', dComponent(i1)
				elseif(fields(i2) .eq. 'layer') then
					write(610,*) '      layer: ', dLayer(i1)
				elseif(fields(i2) .eq. 'activeTime') then
				    writeLine(1:19) = '      activeTime: ['
					write(writeLine(20:37),'(G18.12)') dActTime(1,i1)
					writeLine(38:39) = ', '
					write(writeLine(40:57),'(G18.12)') dActTime(2,i1)
					writeLine(58:58) = ']'
					write(610,*) writeLine(1:58)
				endif
			enddo
		enddo
		
		close(610)
		
	end subroutine writeDVarValues
	
	subroutine writeObjective(fileName,fields,numFields,writeGrad)
	    implicit none
		
		character(len=128), intent(in) :: fileName
		integer, intent(in) :: numFields
		character(len=16), intent(in) :: fields(numFields)
		integer, intent(in) :: writeGrad
		
		character(len=128) :: writeLine
		character(len=32) :: intStr, fltStr
		
		real*8 :: totalObj
		integer :: i1, i2, i3
		
		open(unit=653, file=fileName, status='replace', action='write')
		
		write(653,*) 'objective:'
		
		totalObj = r_0
		do i1 = 1, numObjTerms
		    totalObj = totalObj + objVal(i1)
		enddo
		
        writeLine = '    totalValue: '
        write(writeLine(19:36),'(G18.12)') totalObj
		write(653,*) writeLine

        write(653,*) '    termValues:'
        do i1 = 1, numObjTerms
		    writeLine = '        - value: '
            write(writeLine(18:35),'(G18.12)') objVal(i1)
		    write(653,*) writeLine(1:35)
			do i2 = 1, numFields
			    if(fields(i2) .eq. 'category') then
				    write(653,*) '          category: ', objCategory(i1)
				elseif(fields(i2) .eq. 'operator') then
				    write(653,*) '          operator: ', objOperator(i1)
				elseif(fields(i2) .eq. 'component') then
				    write(653,*) '          component: ', objComponent(i1)
				elseif(fields(i2) .eq. 'layer') then
				    write(653,*) '          layer: ', objLayer(i1)
				elseif(fields(i2) .eq. 'elementSet') then
				    write(653,*) '          elementSet: ', objElSet(i1)
				elseif(fields(i2) .eq. 'nodeSet') then
				    write(653,*) '          nodeSet: ', objNdSet(i1)
				elseif(fields(i2) .eq. 'coefficient') then
				    writeLine = '          coefficient: '
					write(writeLine(23:40),'(G18.12)') objCoef(i1)
				    write(653,*) writeLine(1:40)
				elseif(fields(i2) .eq. 'exponent') then
				    writeLine = '          exponent: '
					write(writeLine(21:38),'(G18.12)') objExp(i1)
				    write(653,*) writeLine(1:38)
				elseif(fields(i2) .eq. 'activeTime') then
				    writeLine(1:19) = '      activeTime: ['
					write(writeLine(20:37),'(G18.12)') dActTime(1,i1)
					writeLine(38:39) = ', '
					write(writeLine(40:57),'(G18.12)') dActTime(2,i1)
					writeLine(58:58) = ']'
					write(653,*) writeLine(1:58)
				endif
			enddo
		enddo
		
		if(writeGrad .eq. 1) then
			write(653,*) 'objectiveGradient:'
			
			do i1 = 1, numDVar
				write(intStr,'(I0)') i1
				write(fltStr,'(G18.12)') dLdD(i1)
				writeLine = '    - [' // intStr // ', ' // fltStr // ']'
				write(653,*) writeLine
			enddo
		endif
		
		close(653)
		
	end subroutine writeObjective
	
!! Binary
	subroutine writeBinarySolution(stepNum)
	    implicit none
		
		integer, intent(in) :: stepNum
		
		character(len=16) :: stepString
		character(len=32) :: fullFileName
		
		write(stepString,'(I0)') stepNum
        fullFileName = 'SolutionHistory' // stepString // '.out'
	    open(unit=17, file=fullFileName, form='unformatted', status='replace', action='write')
		
		write(17) prevTemp(:)
		write(17) prevTdot(:)
		write(17) prevDisp(:)
		if(intVecSize .ge. 1) then
		    write(17) prevIntDisp(:)
		endif
		write(17) prevVel(:)
		write(17) prevAcc(:)
		
		close(17) 
		
	end subroutine writeBinarySolution
	
!! Debugging
	
	subroutine writeModelData(fileName)
	    implicit none
		
		character(len=128), intent(in) :: fileName
		integer :: i1
		
		open(unit=37, file=fileName, status='replace', action='write')
		
		write(37,*) 'nodeList'
		do i1 = 1, numNodes
		    write(37,*) nodeList(:,i1)
		enddo
		write(37,*) 'elementList'
		do i1 = 1, numEls
		    write(37,*) elementList(:,i1)
		enddo
		write(37,*) 'elementSection'
		write(37,*) elementSection(:)
		write(37,*) 'elementType'
		write(37,*) elementType(:)
		write(37,*) 'elementSurfaces'
		write(37,*) elementSurfaces(:,:)
		write(37,*) 'elementSets'
		write(37,*) elementSets(:)
		write(37,*) 'elSetName'
		write(37,*) elSetName(:)
		write(37,*) 'elSetRange'
		write(37,*) elSetRange(:)
		write(37,*) 'nodeSets'
		write(37,*) nodeSets(:)
		write(37,*) 'ndSetName'
		write(37,*) ndSetName(:)
		write(37,*) 'ndSetRange'
		write(37,*) ndSetRange(:)
		write(37,*) 'sectionType'
		write(37,*) sectionType(:)
		write(37,*) 'sectionMatId'
		write(37,*) sectionMatId(:)
		write(37,*) 'sectionMatName'
		write(37,*) sectionMatName(:)
		write(37,*) 'sectionOrient'
		write(37,*) sectionOrient(:,:)
		write(37,*) 'layupMatName'
		write(37,*) layupMatName(:)
		write(37,*) 'layupMatId'
		write(37,*) layupMatId(:)
		write(37,*) 'layupThickness'
		write(37,*) layupThickness(:)
		write(37,*) 'layupAngle'
		write(37,*) layupAngle(:)
		write(37,*) 'secLayupRange'
		write(37,*) secLayupRange(:)
		write(37,*) 'sectionZOffset'
		write(37,*) sectionZOffset(:)
		write(37,*) 'beamProperties'
		write(37,*) beamProperties(:,:)
		write(37,*) 'beamStiffness'
		write(37,*) beamStiffness(:,:)
		write(37,*) 'beamExpLoadCoeff'
		write(37,*) beamExpLoadCoef(:,:)
		write(37,*) 'beamMass'
		write(37,*) beamMass(:,:)
		write(37,*) 'beamThermCond'
		write(37,*) beamThermCond(:,:)
		write(37,*) 'beamSpecHeat'
		write(37,*) beamSpecHeat(:)
		write(37,*) 'materialName'
		write(37,*) materialName(:)
		write(37,*) 'materialDensity'
		write(37,*) materialDensity(:)
		write(37,*) 'materialElastic'
		write(37,*) materialElastic(:,:)
		write(37,*) 'materialStiffMat'
		write(37,*) materialStiffMat(:,:)
		write(37,*) 'materialThermCond'
		write(37,*) materialThermCond(:,:)
		write(37,*) 'materialThermExp'
		write(37,*) materialThermExp(:,:)
		write(37,*) 'materialSpecHeat'
		write(37,*) materialSpecHeat(:)
		write(37,*) 'materialMaxStress'
		write(37,*) materialMaxStress(:,:)
		write(37,*) 'materialMaxStrain'
		write(37,*) materialMaxStrain(:,:)
		write(37,*) 'materialMaxStrnEngy'
		write(37,*) materialMaxStrnEngy(:)
		
		if(numMPC .gt. 0) then
		    write(37,*) 'mpcEqn'
		    write(37,*) mpcEqn(:)
			write(37,*) 'mpcNode'
			write(37,*) mpcNode(:)
			write(37,*) 'mpcDof'
			write(37,*) mpcDof(:)
			write(37,*) 'mpcCoef'
			write(37,*) mpcCoef(:)
			write(37,*) 'mpcRHS'
			write(37,*) mpcRHS(:)
		endif
		
		if(numLds .gt. 0) then
		    write(37,*) 'loadNodes'
		    write(37,*) loadNodes(:)
			write(37,*) 'inputLoads'
	        write(37,*) inputLoads(:,:)
			write(37,*) 'loadsRange'
	        write(37,*) loadsRange(:)
			write(37,*) 'loadsActTime'
	        write(37,*) loadsActTime(:,:)
			write(37,*) 'loadType'
	        write(37,*) loadType(:)
		endif
		
		write(37,*) 'initialDisp'
		write(37,*) initialDisp(:,:)
		write(37,*) 'initialVel'
	    write(37,*) initialVel(:,:)
		write(37,*) 'initalAcc'
	    write(37,*) initialAcc(:,:)
		write(37,*) 'initialTemp'
	    write(37,*) initialTemp(:)
		write(37,*) 'initialTdot'
	    write(37,*) initialTdot(:)
		
		close(37)
	end subroutine writeModelData
	
	subroutine writeDesignVarData(fileName)
	    implicit none
		
		character(len=128), intent(in) :: fileName
		
		open(unit=122, file=fileName, status='replace', action='write')
		
		write(122,*) 'r_dVec'
		write(122,*) r_dVec(:)
		write(122,*) 'dCategory'
		write(122,*) dCategory(:)
		write(122,*) 'dSubCat'
		write(122,*) dSubCat(:)
		write(122,*) 'dComponent'
		write(122,*) dComponent(:)
		write(122,*) 'dLayer'
		write(122,*) dLayer(:)
		write(122,*) 'dActTime'
		write(122,*) dActTime(:,:)
		write(122,*) 'elToD'
		write(122,*) elToD(:)
		write(122,*) 'ndToD'
		write(122,*) ndToD(:)
		write(122,*) 'dToNd'
		write(122,*) dToNd(:)
		write(122,*) 'elToDComp'
		write(122,*) elToDComp(:)
		write(122,*) 'dToElComp'
		write(122,*) dToElComp(:)
		write(122,*) 'elToDRange'
		write(122,*) elToDRange(:)
		write(122,*) 'ndToDRange'
		write(122,*) ndToDRange(:)
		write(122,*) 'dToNdRange'
		write(122,*) dToNdRange(:)
		write(122,*) 'elToDCRange'
		write(122,*) elToDCRange(:)
		write(122,*) 'dToElCRange'
		write(122,*) dToElCRange(:)
		write(122,*) 'elToCoef'
		write(122,*) elToCoef(:)
		write(122,*) 'ndToCoef'
		write(122,*) ndToCoef(:)
		write(122,*) 'dToNdCoef'
		write(122,*) dToNdCoef(:)		
		
		close(122)
	end subroutine writeDesignVarData
	
	subroutine writeObjectiveData(fileName)
	    implicit none
		
		character(len=128), intent(in) :: fileName
		
		open(unit=220, file=fileName, status='replace', action='write')
		
		write(220,*) 'objCategory'
		write(220,*) objCategory(:)
		write(220,*) 'objOperator'
	    write(220,*) objOperator(:)
		write(220,*) 'objActTime'
		write(220,*) objActTime(:,:)
		write(220,*) 'objComponent'
		write(220,*) objComponent(:)
		write(220,*) 'objLayer'
		write(220,*) objLayer(:)
		write(220,*) 'objCoef'
		write(220,*) objCoef(:)
		write(220,*) 'objExp'
		write(220,*) objExp(:)
		write(220,*) 'objElSet'
		write(220,*) objElSet(:)
		write(220,*) 'objNdSet'
		write(220,*) objNdSet(:)
		write(220,*) 'objTgtTag'
		write(220,*) objTgtTag(:)
		write(220,*) 'objTgtVal'
		write(220,*) objTgtVal(:)
		write(220,*) 'objTgtRange'
		write(220,*) objTgtRange(:)
		
		close(220)
		
	end subroutine writeObjectiveData
	
	
end module AStrO_output