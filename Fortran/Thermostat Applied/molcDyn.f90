program main

    use class_module

    implicit none
    integer :: i,j,k,monteCarloSteps = 20000
    real*8 :: pe, TE, ii,jj,kk
    !Initializing everyhing
    call posInit()
    call initForce()
    call initVelocities()
	print *, "20000 Steps Simulation"
    open(20, file = 'KE.dat')
    open(30, file = 'PE.dat')
    open(40, file = 'TE.dat')
    !open(50, file = 'diffT.dat')
    thermoCounter = 500
    do j = 1, monteCarloSteps
      call verletAlgorithm()
        TE = kineticEnergy + potential
        write(20, *) j,kineticEnergy/dfloat(particles)
        write(30, *) j,potential/dfloat(particles)
        write(40, *) j, TE/dfloat(particles)
        updateCounter = updateCounter + 1
    end do

    close(20)
    close(30)
    close(40)
	print *, "Simulation Completed"
    !close(50)
end program main
