!> The module where all the functions and subroutines related to the wavefunctions,
!! their derivatives, and the densities that are constructed from them
module wavefunctions
use init
implicit none

contains
  !> The densities (rho, tau, J) are calculated here
  subroutine build_densities
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
           if (sortenergies(i,iq) < - small) then
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

end module
