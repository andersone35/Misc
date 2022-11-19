program testSectionEls
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
	
	real*8 :: orient(3,3), instOri(3,3), rot(3)
	complex*16 :: c_orient(3,3), c_instOri(3,3), c_rot(3)
	
	real*8 :: allInstOri(3,240)
	complex*16 :: c_allInstOri(3,240)
	
	real*8 :: instU(6,11), globU(6,11), globNds(3,10)
	complex*16 :: c_instU(6,11), c_globU(6,11), c_globNds(3,10)
	
	real*8 :: instUCopy(6,11)

    integer :: numNds, dofTable(2,33)
	
	integer :: i1, i2, i3, i4, i5, i6, i7
	
	open(unit=1, file='sectionTestOut.txt', status='replace', action='write')
	
	orient(:,1) = (/1d0,0d0,0d0/)
	orient(:,2) = (/0d0,1d0,0d0/)
	orient(:,3) = (/0d0,0d0,1d0/)
	c_orient(:,:) = c_1*orient(:,:)
	
	rot(:) = (/0.7d0,0.5d0,0d0/)
	c_rot(:) = c_1*rot(:)
	
	call r_rotateAlpha(instOri,orient,rot,0,0)
	
	write(1,*) 'instOri:'
	do i1 = 1, 3
	    write(1,*) instOri(i1,:)
	enddo
	
	do i1 = 1, 3
	    call r_rotateAlpha(instOri,orient,rot,i1,0)
		write(1,*) 'd(instOri)/d', i1
		do i2 = 1, 3
		    write(1,*) instOri(i2,:)
		enddo
		c_rot(i1) = c_rot(i1) + compStep
		write(1,*) 'complex check'
		call c_rotateAlpha(c_instOri,c_orient,c_rot,0,0)
		c_rot(i1) = c_rot(i1) - compStep
		do i2 = 1, 3
		    write(1,*) compStepInv*imag(c_instOri(i2,1)), compStepInv*imag(c_instOri(i2,2)), compStepInv*imag(c_instOri(i2,3))
		enddo
	enddo
	
	do i1 = 1, 3
	    do i3 = 1, 3
			call r_rotateAlpha(instOri,orient,rot,i1,i3)
			write(1,*) 'd(instOri)/d', i1, i3
			do i2 = 1, 3
				write(1,*) instOri(i2,:)
			enddo
			c_rot(i1) = c_rot(i1) + compStep
			write(1,*) 'complex check'
			call c_rotateAlpha(c_instOri,c_orient,c_rot,i3,0)
			c_rot(i1) = c_rot(i1) - compStep
			do i2 = 1, 3
				write(1,*) compStepInv*imag(c_instOri(i2,1)), compStepInv*imag(c_instOri(i2,2)), compStepInv*imag(c_instOri(i2,3))
			enddo
		enddo
	enddo

    !! Check instantaneous dof
	numNds = 4

    i3 = 1	
	do i1 = 1, 6
	    do i2 = 1, numNds
		    dofTable(1,i3) = i1
			dofTable(2,i3) = i2
			i3 = i3 + 1
		enddo
	enddo
	
	dofTable(:,25) = (/3,7/)
	dofTable(:,26) = (/3,8/)
	dofTable(:,27) = (/3,9/)
	dofTable(:,28) = (/3,10/)
	
	! dofTable(:,25) = (/1,5/)
	! dofTable(:,26) = (/1,6/)
	! dofTable(:,27) = (/2,5/)
	! dofTable(:,28) = (/2,6/)
	! dofTable(:,29) = (/3,7/)
	! dofTable(:,30) = (/3,8/)
	! dofTable(:,31) = (/3,9/)
	! dofTable(:,32) = (/3,10/)
	
	globNds(:,:) = r_0
	globNds(:,1) = (/0d0,-0.5d0,0d0/)
	globNds(:,2) = (/1d0,-0.5d0,0d0/)
	globNds(:,3) = (/1d0,0.5d0,0d0/)
	globNds(:,4) = (/0d0,0.5d0,0d0/)
	
	c_globNds(:,:) = c_1*globNds(:,:)
	
	globU(:,:) = r_0
	globU(:,1) = (/0d0,0d0,0d0,0d0,-0.99d0,0d0/)
	globU(:,2) = (/-0.45d0,0d0,0.85d0,0d0,-1.01d0,0d0/)
	globU(:,3) = (/-0.45d0,0d0,0.85d0,0d0,-1.01d0,0d0/)
	globU(:,4) = (/0d0,0d0,0d0,0d0,-0.99d0,0d0/)
	globU(3,7) = -0.01d0
	globU(3,9) = -0.01d0
	
	c_globU(:,:) = c_1*globU(:,:)
	
	rot(:) = r_0
	do i1 = 1, 4
	    rot(:) = rot(:) + (r_1/numNds)*globU(4:6,i1)
	enddo
	
	call r_getInstOrient(allInstOri,orient,globU,numNds)
	! write(1,*) 'allInstOri:'
	! do i1 = 0, 3
	    ! do i2 = 0, 3
		    ! i3 = 12*i1 + 3*i2
			! write(1,*) 'dv: ', i1, i2
			! do i4 = 1, 3
			    ! write(1,*) allInstOri(i4,i3+1:i3+3)
			! enddo
		! enddo
	! enddo
	! write(1,*) ' '
	
	call r_getInstDof(instU,globU,globNds,allInstOri,orient,numNds,dofTable,0,0)
	
	write(1,*) 'instU:'
	do i1 = 1, 6
	    write(1,*) instU(i1,:)
	enddo
	write(1,*) ' '
	
	c_allInstOri = c_1*allInstOri
	c_orient = c_1*orient
	c_rot = c_1*rot
	do i1 = 1, 28
	    call r_getInstDof(instU,globU,globNds,allInstOri,orient,numNds,dofTable,i1,0)
		write(1,*) 'd(instU)d', i1
		do i2 = 1, 6
		    write(1,*) instU(i2,:)
		enddo
		i4 = dofTable(1,i1)
		i5 = dofTable(2,i1)
        c_globU(i4,i5) = c_globU(i4,i5) + compStep
		if(i4 .gt. 3) then
			call c_getInstOrient(c_allInstOri,c_orient,c_globU,numNds)
		endif
        call c_getInstDof(c_instU,c_globU,c_globNds,c_allInstOri,c_orient,numNds,dofTable,0,0)
        c_globU(i4,i5) = c_globU(i4,i5) - compStep
		c_allInstOri = c_1*allInstOri
		do i2 = 1,6
		    do i3 = 1, 11
			    instUCopy(i2,i3) = compStepInv*imag(c_instU(i2,i3))
				if(abs(instUCopy(i2,i3) - instU(i2,i3)) .gt. 1e-10) then
				    write(1,*) 'discrepancy: ', i2, i3
				endif
			enddo
		enddo
		write(1,*) 'complex check:'
		do i2 = 1, 6
		    write(1,*) instUCopy(i2,:)
		enddo
		write(1,*) ' '
		do i7 = i1, 28
			call r_getInstDof(instU,globU,globNds,allInstOri,orient,numNds,dofTable,i1,i7)
			write(1,*) 'd(instU)d', i1, i7
			do i2 = 1, 6
				write(1,*) instU(i2,:)
			enddo
			i4 = dofTable(1,i1)
			i5 = dofTable(2,i1)
			c_globU(i4,i5) = c_globU(i4,i5) + compStep
			if(i4 .gt. 3) then
				call c_getInstOrient(c_allInstOri,c_orient,c_globU,numNds)
			endif
			call c_getInstDof(c_instU,c_globU,c_globNds,c_allInstOri,c_orient,numNds,dofTable,i7,0)
			c_globU(i4,i5) = c_globU(i4,i5) - compStep
			c_allInstOri = c_1*allInstOri
			do i2 = 1,6
				do i3 = 1, 11
					instUCopy(i2,i3) = compStepInv*imag(c_instU(i2,i3))
					if(abs(instUCopy(i2,i3) - instU(i2,i3)) .gt. 1e-6) then
						write(1,*) 'discrepancy: ', i2, i3
					endif
				enddo
			enddo
			write(1,*) 'complex check:'
			do i2 = 1, 6
				write(1,*) instUCopy(i2,:)
			enddo
			write(1,*) ' '		
		enddo
	enddo
	
	close(1)
	
end program