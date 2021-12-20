!> The module where all the functions and subroutines related to the fields
!! and energies are located.
module fields
use init
use wavefunctions
implicit none

contains
  !> build_fields builds the fields and populates the field arrays
  subroutine build_fields
    integer :: ir, iq, ir2
    real(wp), dimension(0:nbox,2) :: ucnew,umrnew,uddnew,usonew,ucsonew,dumrnew,d2umrnew
    real(wp), dimension(0:nbox) :: ucoulnew
    real(wp) :: tot1=0.0d0,tot2=0.0d0
    real(wp) :: xmix, ymix,a,b,c

    xmix = 0.4
    ymix = 1.-xmix

    do iq = 1,2
     do ir = 0,nbox

    !!Central Field U(r)
             ucnew(ir,iq) = 2*(a0r0-a1r1)*rho(ir,3) + 4*a1r1 * rho(ir,iq)  &
                            + (a0tau0-a1tau1) *tau(ir,3)+ 2 *a1tau1*tau(ir,iq) &
                            +2*( a0r0p-a1r1p )*laprho(ir,3) + 4 *a1r1p * laprho(ir,iq)
    !!Part of U(r) coming from so)
             ucsonew(ir,iq) = (cso0-cso1 ) *(djsc(ir,3) + 2 * jsc(ir,3)/mesh(ir) ) &
                            + 2 *cso1 * ( djsc(ir,iq) + 2 * jsc(ir,iq) / mesh(ir) )
    !!Mq(r) contributions
             umrnew(ir,iq) = hbar22m*cmcorr+(a0tau0-a1tau1)*rho(ir,3) + 2 * a1tau1*rho(ir,iq)
             dumrnew(ir,iq) = (a0tau0-a1tau1)*drho(ir,3) + 2 * a1tau1*drho(ir,iq)
             d2umrnew(ir,iq) = (a0tau0-a1tau1)*ddrho(ir,3) + 2 * a1tau1*ddrho(ir,iq)
    !! t3 part of U(r)
              uddnew(ir,iq) = ( 2 + sig ) * (cddr0-cddr1)*rho(ir,3)**(sig+1)  &
                            +2*sig*cddr1*(rho(ir,1)**2+rho(ir,2)**2)*rho(ir,3)**(sig-1.) &
                            + 4 * cddr1 * rho(ir,iq) * rho(ir,3)**sig
     !!spin-orbit part
             usonew(ir,iq) = - (cso0-cso1 )*drho(ir,3)/mesh(ir) &
                             - 2 *cso1 * drho(ir,iq) / mesh(ir)
             if (j2terms) then
             usonew(ir,iq) = usonew(ir,iq)-(a0t0-a1t1) *jsc(ir,3) / mesh(ir) &
                             - 2 *a1t1 * jsc(ir,iq) / mesh(ir)
             end if
     !!coulomb
      if (iq==2) then
            tot1=0.0d0
            tot2=0.0d0
        do ir2=0,ir
         tot1=tot1+rho(ir2,2)*(mesh(ir2)**2)
        end do

        do ir2=ir+1,nbox
         tot2=tot2+rho(ir2,2)*mesh(ir2)
        end do
        ucoulnew(ir)=4.0d0*pi*e2*(tot1/mesh(ir)&
        + tot2)*h - e2*(3./pi)**(1./3.)*rho(ir,2)**(1./3.)
      end if
     end do
    end do


    uc(:,:) = ucnew(:,:)*xmix + uc(:,:)*ymix
    ucso(:,:) = ucsonew(:,:)*xmix + ucso(:,:)*ymix
    umr(:,:) = umrnew(:,:)*xmix + umr(:,:)*ymix
    udd(:,:) = uddnew(:,:)*xmix + udd(:,:)*ymix
    uso(:,:) = usonew(:,:)*xmix + uso(:,:)*ymix
    ucoul(:) = ucoulnew(:)*xmix + ucoul(:)*ymix
    dumr(:,:) = dumrnew(:,:)*xmix + dumr(:,:)*ymix
    d2umr(:,:) = d2umrnew(:,:)*xmix + d2umr(:,:)*ymix

    if(icoul==0) ucoul(:) = 0.
  end subroutine build_fields
  !> totenergy calculates the total energy using both the koopman theorem and
  !! the integration of the energy density
  subroutine totenergy
    real(wp) :: kinetic(1:2),val
    integer :: iq,i,l,is
    !Functional variables
    real(wp) :: ekin,ecould,ecoulex,ecentr,edd, eso
    real(wp) :: tmp(0:nbox),rg(0:nbox)


    !method using the functional energy
    totfunct = 0.0_wp
    !!!Kinetic energy
      ekin = hbar22m * cmcorr * 4._wp*pi*h*sum( tau(:,3)*mesh(:)**2)
    !!!Coulomb energy
      ecould = 0.0_wp
      ecoulex = 0.0_wp
    if (icoul==1) then
    rg = mesh
    tmp = 0.0_wp
    do i = 0, nbox
       rg(0:i) = mesh(i)
       tmp(i) = 4._wp*pi*h*sum( mesh(:)**2 / rg * rho(:,2) )
    end do
    ecould = 4._wp*pi*h*sum( mesh(:)**2 *rho(:,2)* tmp )*e2 / 2._wp
    ecoulex = - (3.0_wp/4.0_wp)*e2*4._wp*pi*h*sum( mesh(:)**2*rho(:,2)**(4._wp/3._wp) )*( 3 / pi)**(1._wp/3._wp)
    end if
    !!!Central energy
    ecentr = 0.0_wp
    ecentr = a0r0 * 4._wp*pi*h*sum( mesh(:)**2*rho(:,3)**2) &
           & + a1r1*4._wp*pi*h*sum(mesh(:)**2*rho(:,4)**2) &
    	   & + a0r0p*4._wp*pi*h*sum(mesh(:)**2*rho(:,3)*laprho(:,3)) &
    	   & + a1r1p*4._wp*pi*h*sum(mesh(:)**2*rho(:,4)*laprho(:,4)) &
    	   & + a0tau0*4._wp*pi*h*sum(mesh(:)**2*rho(:,3)*tau(:,3))  &
    	   & + a1tau1*4._wp*pi*h*sum(mesh(:)**2*rho(:,4)*tau(:,4))
    if (j2terms) then
    ecentr = ecentr - 0.5 * a0t0 *4._wp*pi*h* sum(mesh(:)**2*jsc(:,3)**2)&
    	& - 0.5 * a1t1 * 4._wp*pi*h*sum(mesh(:)**2*jsc(:,4)**2)
    end if
    !!!Density Dependent Energy
    edd = 4._wp*pi*h*sum(mesh(:)**2*( cddr0 * rho(:,3)**2 + cddr1 * rho(:,4)**2 )  &
         * rho(:,3)**sig )
    !!Spin-orbit Energy
    eso = 0.0_wp
    eso = cso0 *4._wp*pi*h* sum(mesh(:)**2*rho(:,3)* &
        &  ( djsc(:,3) + 2 * jsc(:,3) / mesh(:))) &
        & + cso1 *4._wp*pi*h*sum(mesh(:)**2*rho(:,4)* &
        &  ( djsc(:,4) + 2 * jsc(:,4) / mesh(:)))

    !!!Total Energy with the functional
    totfunct = ekin+ecould+ecoulex+ecentr+edd+eso

    !Energy calculation using Koopman's theorem
    totalenergy = 0.
    do iq=1,2
      kinetic(iq) = 4*pi*h*cmcorr*hbar22m*sum(mesh(:)**2 * tau(:,iq))
      do i=1,nn
        l = sortstates(i,2,iq)
        is= sortstates(i,3,iq)
        if (sortenergies(i,iq) < -small) then
          val = (2*(l+spin(is))+1)
          !val = 2*l+1
          if (l==0) val = 2.
          totalenergy = totalenergy + sortenergies(i,iq)*val
        end if
      end do
    end do
    totalenergy = totfunct !(totalenergy + kinetic(1) + kinetic(2))/2.&
                !-4*pi*h*t3*sum(mesh(:)**2 * (rho(:,3)**sig &
                !*(rho(:,3)**2 -(rho(:,1)**2 +rho(:,2)**2)/2.)))/24.&
                !-e2*(3./pi)**(1./3.)*pi*h*sum(mesh(:)**2*rho(:,2)**(4./3.))
                !-4*h*pi*t3*0.125*sum(mesh(:)**2 * rho(:,3)*rho(:,1)*rho(:,2))
    totalkinetic = ekin!sum(kinetic(:))
    !print *,totfunct - totalenergy
  end subroutine totenergy
  !> energy_sort sorts the occupied single particle states
  subroutine energy_sort
    integer :: n, l, iq, is,k,i,nfill,nfull
    real(wp) :: temp,j
    integer, dimension(1:3) :: state
    ! This sorts the energies
    sortenergies = small
    sortstates = 0
    do iq =1,2
     nfull = nn
     if (iq == 2) nfull = np
     k = 1
     nfill = 0
     do while(nfill.lt.nfull)
        temp = 0._wp
        state = 1
       do n = 1,lmax
        do l = 0,lmax
          do is=1,2
            if(temp > energies(n,l,is,iq)) then
            temp = energies(n,l,is,iq)
            state(1) = n
            state(2) = l
            state(3) = is
            end if
          end do
        end do
       end do
       sortenergies(k,iq) = energies(state(1),state(2),state(3),iq)
       do i = 1,3
       sortstates(k,i,iq) = state(i)
       end do
       energies(state(1),state(2),state(3),iq) = 0.0_wp
       if (state(2)==0) energies(state(1),state(2),:,iq) = 0.0_wp
         j = state(2) + spin(state(3))
       if (state(2)==0) j = 0.5
       nfill = nfill + INT(2*j+1)
        if (nfill.gt.nfull) then
         sortenergies(k-1,iq) = sortenergies(k,iq)
         sortenergies(k,iq) = 0.0_wp
         do i = 1,3
          sortstates(k-1,i,iq) = sortstates(k,i,iq)
          sortstates(k,i,iq) = 1
         end do
        end if
       k = k+1
     end do
    !if (nfill.ne.nfull) then
    !write (6,*) 'Not the correct number of particles'
    !stop
    !end if
    end do
  end subroutine energy_sort

end module fields
