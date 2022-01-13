!> The module where all the functions and subroutines related to the wavefunctions,
!! their derivatives, and the densities that are constructed from them
module wavefunctions
use init
implicit none

contains
  !> The densities (rho, tau, J) are calculated here
  subroutine build_densities
    integer :: npr,iq,ir,i,l,n,is
    real(wp) :: j, dwfr(0:nbox,lmax,0:lmax,2,2)
    do iq =1,2

        if (iq == 1) then
        npr = nn
        else
        npr = np
        end if
        rho(:,iq)=small
        tau(:,iq)=small
        jsc(:,iq)=small
        do i = 1, npr
           if (sortenergies(i,iq) < - small) then
             n = sortstates(i,1,iq)
             l = sortstates(i,2,iq)
             is= sortstates(i,3,iq)
             j = sortstates(i,2,iq) + spin(sortstates(i,3,iq))
             if (sortstates(i,2,iq) == 0) j = 0.5
             dwfr = 0.0
             call derivative_sca( nbox, wfr(1:nbox,n,l,is,iq), dwfr(1:nbox,n,l,is,iq), l )
             do ir=1,nbox
               !tau(ir,iq) = tau(ir,iq) + (2*j+1)*((dwavefunction(ir,n,l,is,iq)&
               !-wfr(ir,n,l,is,iq)/mesh(ir))**2+l*(l+1)*(wfr(ir,n,l,is,iq)**2)&
               !/mesh(ir)**2)/(4*pi*mesh(ir)**2)
               tau(ir,iq) = tau(ir,iq) + (2*j+1)*((dwfr(ir,n,l,is,iq)&
               -wfr(ir,n,l,is,iq)/mesh(ir))**2+l*(l+1)*(wfr(ir,n,l,is,iq)**2)&
               /mesh(ir)**2)/(4*pi*mesh(ir)**2)

               rho(ir,iq) = rho(ir,iq) + (2*j+1)*wfr(ir,sortstates(i,1,iq),sortstates(i,2,iq),sortstates(i,3,iq),iq)&
               *wfr(ir,sortstates(i,1,iq),sortstates(i,2,iq),sortstates(i,3,iq),iq) / (4*pi*mesh(ir)**2)

               jsc(ir,iq) = jsc(ir,iq) + (2*j+1)*(j*(j+1)-sortstates(i,2,iq)*(sortstates(i,2,iq)+1)-0.75)&
               *wfr(ir,sortstates(i,1,iq),sortstates(i,2,iq),sortstates(i,3,iq),iq)**2&
               /(4*pi*mesh(ir)**3)
              end do
              rho(0,iq) = rho(1,iq)
              !tau(2,iq) = tau(3,iq)
              tau(1,iq) = tau(2,iq)
              tau(0,iq) = tau(1,iq)
            end if
        end do
     end do
     rho(:,3)=rho(:,1) + rho(:,2)
     rho(0,3) = rho(1,3)
     rho(:,4)=rho(:,1) - rho(:,2)
     tau(:,3)=tau(:,1) + tau(:,2)
     tau(:,4)=tau(:,1) - tau(:,2)
     jsc(:,3)=jsc(:,1) + jsc(:,2)
     jsc(:,4)=jsc(:,1) - jsc(:,2)

     call ddensities
  end subroutine build_densities

  subroutine restart_densities
    integer :: npr,iq,ir,i,l,n,is
    real(wp) :: j

    do iq =1,2

        if (iq == 1) then
        npr = nn
        else
        npr = np
        end if
        rho(:,iq)=small
        tau(:,iq)=small
        jsc(:,iq)=small
        do i = 1, npr
             n = sortstates(i,1,iq)
             l = sortstates(i,2,iq)
             is= sortstates(i,3,iq)
             j = sortstates(i,2,iq) + spin(sortstates(i,3,iq))
             if (sortstates(i,2,iq) == 0) j = 0.5
             do ir=1,nbox
               tau(ir,iq) = tau(ir,iq) + (2*j+1)*((dwavefunction(ir,n,l,is,iq)&
               -wfr(ir,n,l,is,iq)/mesh(ir))**2+l*(l+1)*(wfr(ir,n,l,is,iq)**2)&
               /mesh(ir)**2)/(4*pi*mesh(ir)**2)

               rho(ir,iq) = rho(ir,iq) + (2*j+1)*wfr(ir,sortstates(i,1,iq),sortstates(i,2,iq),sortstates(i,3,iq),iq)&
               *wfr(ir,sortstates(i,1,iq),sortstates(i,2,iq),sortstates(i,3,iq),iq) / (4*pi*mesh(ir)**2)

               jsc(ir,iq) = jsc(ir,iq) + (2*j+1)*(j*(j+1)-sortstates(i,2,iq)*(sortstates(i,2,iq)+1)-0.75)&
               *wfr(ir,sortstates(i,1,iq),sortstates(i,2,iq),sortstates(i,3,iq),iq)**2&
               /(4*pi*mesh(ir)**3)
              end do
              rho(0,iq) = rho(1,iq)
              !tau(2,iq) = tau(3,iq)
              tau(1,iq) = tau(2,iq)
              tau(0,iq) = tau(1,iq)
        end do
     end do
     rho(:,3)=rho(:,1) + rho(:,2)
     rho(0,3) = rho(1,3)
     rho(:,4)=rho(:,1) - rho(:,2)
     tau(:,3)=tau(:,1) + tau(:,2)
     tau(:,4)=tau(:,1) - tau(:,2)
     jsc(:,3)=jsc(:,1) + jsc(:,2)
     jsc(:,4)=jsc(:,1) - jsc(:,2)

     call ddensities
  end subroutine restart_densities

  ! Adapted from derivatives in HFBRAD
  !! Derivative of a...
  !

  ! ...salar function

  subroutine derivative_sca( n, f, df, l )
    implicit none
    integer, intent(in) :: n
    real (wp), intent(in) :: f(n)
    real (wp), intent(inout) :: df(n)
    integer, intent(in) :: l
    real (wp) :: sig, h_12, h_60
    integer :: i, boundary_condition
    boundary_condition = 2
    h_12 = h*12.0_wp
    h_60 = h*60.0_wp
    sig = ( modulo(l,2) - 0.5_wp ) * 2
    df(1) = ( 8.0_wp * f(2) - f(3) + sig * f(1) ) / h_12
    df(2) = ( 45.0_wp * ( f(3) - f(1) ) - 9.0_wp * f(4) &
         + f(5) - sig * f(1) ) / h_60
    df(3) = ( 45.0_wp * ( f(4) - f(2) ) - 9.0_wp * ( f(5) - f(1) ) &
         + f(6) ) / h_60
    !
    if ( boundary_condition == 0 .or.  &
         !boundary_condition == 2 .and. l == 2 * ( l / 2 ) .or. &
         (boundary_condition == 2 .and. mod(l,2) == 0) .or. &
         (boundary_condition == 3 .and. l /= 2 * ( l / 2 )) ) then
       df(n) = ( -8.0_wp * f(n-1) + f(n) + f(n-2) ) / h_12
       df(n-1) = ( 45.0_wp * ( f(n) - f(n-2) ) + 9.0_wp * f(n-3) &
            - f(n) - f(n-4) ) / h_60
       df(n-2) = ( 45.0_wp * ( f(n-1) - f(n-3) ) &
            - 9.0_wp * ( f(n) - f(n-4) ) - f(n-5) ) / h_60
    end if
    if ( boundary_condition == 1 .or.  &
         (boundary_condition == 2 .and. mod(l,2) /= 0) .or. &
         (boundary_condition == 3 .and. l == 2 * ( l / 2 )) ) then
       df(n) = ( - 54._wp * f(n-1) + 45._wp * f(n) + 10._wp * f(n-2) - f(n-3) ) / h_60
       df(n-1) = ( 36._wp * f(n) + f(n-1) - 45._wp * f(n-2) &
            + 9._wp * f(n-3) - f(n-4) ) / h_60
       df(n-2) = ( 45._wp * ( f(n-1) - f(n-3) ) - 8._wp * f(n)  &
            + 9._wp * f(n-4) - f(n-5) ) / h_60
    end if
    !
    do i = 4, n - 3
       df(i) = ( 45.0_wp * ( f(i+1) - f(i-1) )  &
            - 9.0_wp * ( f(i+2) - f(i-2) ) &
            + f(i+3) - f(i-3) ) / h_60
    end do
    !
  end subroutine derivative_sca

  !> dwavefunction calculates the derivative of the wavefunctions and a given
  !! point using second order finite difference
  function dwavefunction(ir,n,l,is,iq) result(derv)
    integer, intent(in) :: ir,n,l,is,iq
    real(wp) :: derv

    if(ir == 0) then
      if(mod(l,2)==1) then
        derv = (-wfr(ir+2,n,l,is,iq) + 8*wfr(ir+1,n,l,is,iq) &
               -8*wfr(ir+1,n,l,is,iq) + wfr(ir+2,n,l,is,iq))/(12*h)
      else
        derv = (-wfr(ir+2,n,l,is,iq) + 8*wfr(ir+1,n,l,is,iq) &
               +8*wfr(ir+1,n,l,is,iq) - wfr(ir+2,n,l,is,iq))/(12*h)
      end if
    else if (ir == 1) then
      if(mod(l,2)==1) then
        derv = (-wfr(ir+2,n,l,is,iq) + 8*wfr(ir+1,n,l,is,iq) &
               -8*wfr(ir-1,n,l,is,iq) + wfr(ir,n,l,is,iq))/(12*h)
      else
        derv = (-wfr(ir+2,n,l,is,iq) + 8*wfr(ir+1,n,l,is,iq) &
               +8*wfr(ir-1,n,l,is,iq) - wfr(ir,n,l,is,iq))/(12*h)
      end if
    else if ((ir >= 2) .AND. (ir <= nbox-2)) then
      derv = (-wfr(ir+2,n,l,is,iq) + 8*wfr(ir+1,n,l,is,iq) &
             -8*wfr(ir-1,n,l,is,iq) + wfr(ir-2,n,l,is,iq))/(12*h)
    else if ((ir > nbox-2) .AND. (ir/=nbox)) then
      derv = 0.!(2*wfr(ir+1,n,l,is,iq) + 3*wfr(ir,n,l,is,iq) &
             !-6*wfr(ir-1,n,l,is,iq) + wfr(ir-2,n,l,is,iq))/(6*h)
    else
      derv = 0.
    end if
  end function

  !> ddensities calculates the derivative of rho and J and the second derivative
  !! of rho and stores the value in an array
  subroutine ddensities
    integer :: ir,iq,i
    do iq=1,2
      do ir=0,nbox
        if(ir < 1) then
          djsc(ir,iq) = (-jsc(ir+2,iq) + 8*jsc(ir+1,iq) &
                 +8*jsc(ir+1,iq) - jsc(ir+2,iq))/(12*h)

          drho(ir,iq) = (-rho(ir+2,iq) + 8*rho(ir+1,iq) &
                 -8*rho(ir+1,iq) + rho(ir+2,iq))/(12*h)

          ddrho(ir,iq) = (-rho(ir+2,iq)+16*rho(ir+1,iq)-30*rho(ir,iq)&
                  +16*rho(ir+1,iq)-rho(ir+2,iq))/(12*h**2)
        else if (ir < 2) then
          djsc(ir,iq) = (-jsc(ir+2,iq) + 8*jsc(ir+1,iq) &
                 +8*jsc(ir-1,iq) - jsc(ir,iq))/(12*h)

          drho(ir,iq) = (-rho(ir+2,iq) + 8*rho(ir+1,iq) &
                 -8*rho(ir-1,iq) + rho(ir,iq))/(12*h)

          ddrho(ir,iq) = (-rho(ir+2,iq)+16*rho(ir+1,iq)-30*rho(ir,iq)&
                  +16*rho(ir-1,iq)-rho(ir,iq))/(12*h**2)
        else if ((ir >= 2) .AND. (ir <= nbox-2)) then
          drho(ir,iq) = (-rho(ir+2,iq) + 8*rho(ir+1,iq) &
                 -8*rho(ir-1,iq) + rho(ir-2,iq))/(12*h)

          djsc(ir,iq) = (-jsc(ir+2,iq) + 8*jsc(ir+1,iq) &
                 -8*jsc(ir-1,iq) + jsc(ir-2,iq))/(12*h)

          ddrho(ir,iq) = (-rho(ir+2,iq)+16*rho(ir+1,iq)-30*rho(ir,iq)&
                  +16*rho(ir-1,iq)-rho(ir-2,iq))/(12*h**2)


        else if ((ir > nbox-2) .AND. (ir/=nbox)) then
          drho(ir,iq) = 0.
          ddrho(ir,iq) = 0.
          djsc(ir,iq) = 0.
        else
          drho(ir,iq) = 0.
          ddrho(ir,iq) = 0.
          djsc(ir,iq) = 0.
        end if
      end do
      ddrho(0:3,iq)=ddrho(4,iq)
    end do

    drho(:,3) = drho(:,1) + drho(:,2)
    drho(:,4) = drho(:,1) - drho(:,2)
    djsc(:,3) = djsc(:,1) + djsc(:,2)
    djsc(:,4) = djsc(:,1) - djsc(:,2)
    ddrho(:,3) = ddrho(:,1) + ddrho(:,2)
    ddrho(:,4) = ddrho(:,1) - ddrho(:,2)
    do i =1,4
    laprho(:,i) = ddrho(:,i) + 2._wp/mesh(:)*drho(:,i)
    end do
  end subroutine ddensities

  subroutine restartwfs
    implicit none
    Open(Unit=20, File="wf_numpy.bin", form='unformatted', access="stream")

    Read(20) h, nbox, lmax, nmax, sortstates, sortenergies, wfr

    Close(Unit=20,Status='keep')

  end subroutine restartwfs

  subroutine printwfs

    Open(Unit=19, File="wf.bin", Status='unknown', Form='unformatted', Position='asis')
    Open(Unit=20, File="wf_numpy.bin", status='replace', access="stream")

    Write(19) sortstates, sortenergies, wfr
    Write(20) h, nbox, lmax, nmax, sortstates, sortenergies, wfr

    Close(Unit=19,Status='keep')
    Close(Unit=20,Status='keep')

  end subroutine printwfs

  subroutine printhpsis

    Open(Unit=21, File="hpsi.bin", status='replace', access="stream")

    Write(21) nbox, lmax, nmax, sortenergies, hpsi

    Close(Unit=21,Status='keep')

  end subroutine printhpsis

end module
