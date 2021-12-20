!> Driver module that runs the main iteration procedure
module solver
  use init
  use fields
  use wavefunctions
  implicit none

contains
  !> statichf performs the Hartree-Fock iterations until the convergence criteria
  !! or max number of iterations are met. Within, the Numerov method is used to
  !! solve for the energy of a given state.
  subroutine statichf

    integer :: i, ir, nnodes,l, is, iq, n,iter
    real(wp) :: Etrial, Eupper, Elower, a1, a2, a3, b1, b2, b3, norm, j, diff,oldtotenergy,currentconv
    real(wp), allocatable :: temp1(:)

    oldtotenergy = 0.0
    wfr(:,:,:,:,:) = 0.0
    !Main loop begins here; be careful
    allocate(temp1(0:nbox))
    do iter = 1,itermax
      if(iter>1) then
        call build_fields
        !stop
      end if
      call totenergy
      if(iter==1) totfunct=0.
      oldtotenergy = totfunct

      do iq =1,2
        do is =1,2
          do l =0,lmax
              j = l + spin(is)
              if (l==0) j=0.5
              temp1(:) = -0.25*(dumr(:,iq)/umr(:,iq))**2&
              +(-uc(:,iq) -ucso(:,iq)-udd(:,iq)-uso(:,iq)*(j*(j+1)- l*(l+1) - 0.75) &
              -dumr(:,iq)/mesh(:) + (1-iq)*ucoul(:) - umr(:,iq)*l*(l+1)/mesh(:)**2)/umr(:,iq)&
              -0.5*(d2umr(:,iq)*umr(:,iq) - dumr(:,iq)**2)/(umr(:,iq)**2)
              do n = 1,nnmax
                  Eupper = 1_wp
                  Elower = -100._wp

                  do i=1,1000000
                    ! Trial Energy for Numerov Algorithm
                    Etrial = (Eupper+Elower)/2.0
                    do ir=0,nbox
                      ! Isospin dependent potential using matrices from before
                         potential(ir) = temp1(ir) + Etrial/umr(ir,iq)
                    end do

                    ! Setting initial values for left and right W.F.
                    wfr(nbox,n,l,is,iq) = 0.
                    wfr(nbox-1,n,l,is,iq) = 1.
                    wfl(0,n,l,is,iq) = 0.
                    wfl(1,n,l,is,iq) = mesh(1)**(l+1)

                    nnodes = 0
                    ! Our Friend, Numerov
                    do ir=nbox-1,1,-1
                      ! These are for the right propagating W.F.'s
                      a1 = 2.0 * (1.0 - (5.0/12.0) * h**2 * potential(ir))
                      a2 = (1.0 + (1.0/12.0) * h**2 * potential(ir+1))
                      a3 = (1.0 + (1.0/12.0) * h**2 * potential(ir-1))
                      ! These are for the left propagating W.F.'s
                      b1 = 2.0 * (1.0 - (5.0/12.0) * h**2 * potential(nbox-ir))
                      b2 = (1.0 + (1.0/12.0) * h**2 * potential(nbox - ir-1))
                      b3 = (1.0 + (1.0/12.0) * h**2 * potential(nbox - ir+1))

                      wfr(ir-1,n,l,is,iq) = (a1*wfr(ir,n,l,is,iq) - a2*wfr(ir+1,n,l,is,iq))/a3
                      wfl(nbox - ir+1,n,l,is,iq) = (b1*wfl(nbox - ir,n,l,is,iq) - b2*wfl(nbox - ir-1,n,l,is,iq))/b3
                      ! Setting coefficient for left W.F. scaling W.R.T. right
                      if (ir == njoin) diff = wfr(ir,n,l,is,iq) / wfl(ir,n,l,is,iq)
                      ! Scale all the left W.F. from njoin -> 0
                      if (ir <= njoin) then
                        wfl(ir,n,l,is,iq) = diff * wfl(ir,n,l,is,iq)
                      end if
                      ! Count the nodes
                      if((wfl(nbox-ir,n,l,is,iq)*wfl(nbox-ir+1,n,l,is,iq) < 0)) nnodes = nnodes+1
                    end do
                    ! Unite the left and right
                    wfr(0:njoin,n,l,is,iq) = wfl(0:njoin,n,l,is,iq)

                    if (nnodes > n-1) then
                      Eupper = Etrial
                    else if (nnodes <= n-1) then
                      Elower = Etrial
                    end if
                    ! If the energies converge on something that is not vpb, select that solution
                    if (abs(Eupper - Elower) < conv) then
                      if (Etrial < 0 .AND. Etrial > -100.+.01) then
                            wfr(:,n,l,is,iq) = wfr(:,n,l,is,iq)/sqrt(umr(:,iq))
                            vocc(n,l,is,iq) = 2*l+1
                            energies(n,l,is,iq) = etrial
                            norm = sqrt(sum(h*wfr(:,n,l,is,iq)*wfr(:,n,l,is,iq)))
                            wfr(:,n,l,is,iq) = wfr(:,n,l,is,iq)/norm
                         end if
                      ! Leave this place
                      exit
                    end if
                  end do
                end do
            end do
          end do
        end do
        call energy_sort
        call build_densities
        call totenergy
        currentconv = abs(totfunct - oldtotenergy)
        write (6,*) "Iteration:",iter,"Convergence:",currentconv
        if (currentconv<conv) exit

       end do
  end subroutine statichf

end module solver
