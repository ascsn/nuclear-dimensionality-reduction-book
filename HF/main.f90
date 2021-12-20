program main
    use init
    use solver
    use wavefunctions
    use fields
    integer :: iq,i,npr,ir
    real(wp) :: j,rneu,rpro,rtot,rch,npro,nneu
    open(unit=5,file='in',status='old',form='formatted')
    open(unit=6,file='out',form='formatted')
    open(unit=13,file='plt',form='formatted')
    open(unit=14,file='densities',form='formatted')
    open(unit=15,file='fields',form='formatted')

    call init_params
    call init_grids
    call init_wavefunctions
    call init_fields
    call statichf

    ! Writing the single particle states to 'out'
    write(6,*) "Single Particle States | "
    do iq =1,2

      if (iq == 1) then
      write(6,*) "Neutrons | "
      npr = nn
      else
      write(6,*) "Protons  | "
      npr = np
      end if
      do i = 1, npr
            if (sortenergies(i,iq) < - small) then
              j = sortstates(i,2,iq) + spin(sortstates(i,3,iq))
              if (sortstates(i,2,iq) == 0) j = 0
              write (6,*) "n=", sortstates(i,1,iq), &
                          &"l=",sortstates(i,2,iq),&
                          &"is=",sortstates(i,3,iq),&
                          &"Energy=", sortenergies(i,iq)
            end if
      end do
    end do
    do ir=0,nbox
      write (13,*) ir*h, wfr(ir,1,0,1,1), wfr(ir,1,1,1,1)
    end do
    ! Particle by way of density integration
    call totenergy
    nneu = sum(4*pi*h*mesh(:)**2*rho(:,1))
    npro = sum(4*pi*h*mesh(:)**2 *rho(:,2))
    write(6,*) "Total Neutrons =", nneu, &
    "Total Protons =", npro
    write(6,*) "Total Energy =", totalenergy,"Total Kinetic Energy:",totalkinetic
    !do ir = 0,nbox
    !  write(14,*) ir*h, rho(ir,1),rho(ir,2),rho(ir,3),drho(ir,1),ddrho(ir,1),tau(ir,1),tau(ir,2),jsc(ir,1)!,rho(ir,4)
    !end do
    rneu = sqrt(4*pi*h*sum(mesh(:)**4 * rho(:,1))/nneu)
    rpro = sqrt(4*pi*h*sum(mesh(:)**4 * rho(:,2))/npro)
    rch = rpro + 0.8**2
    write (6,*) "Neutron RMS Radius(fm) =", rneu
    write (6,*) "Proton RMS Radius(fm) =", rpro
    write (6,*) "Charge Radius(fm) =", rch

    close(13)
    close(6)
    close(5)
    close(14)
    close(15)

end program main
