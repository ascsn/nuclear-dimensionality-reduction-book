subroutine gethpsi(h, nbox, lmax, nmax, sortstates, sortenergies, wfr, hpsi, &
  t0,x0,t1,x1,t2,x2,t3,x3,sig,w0,ij2terms,icoul,icm)
  implicit none
  integer, parameter :: wp=kind(1.0d0)
  !real(wp), parameter :: pi = 3.14159265358979_wp
  real(wp), parameter :: pi = 3.141592653589793_wp
  real(wp), parameter :: e2 = 1.4399784_wp
  real(wp), parameter :: hbar = 6.582119E-22_wp
  real(wp), parameter :: a = 0.67_wp
  real(wp), parameter :: vso = 23_wp
  real(wp) :: t0,x0,t1,x1,t2,x2,t3,x3,sig,w0
  real(wp) :: cmcorr, xmix
  real(wp) :: h,conv,hbar22m,v0,nrad,vpb(2),r0,small,spin(2),totalenergy,totfunct,totalkinetic
  real(wp) :: a0r0,a1r1,a0s0,a1s1,a0tau0,a1tau1,a0t0,a1t1,a0r0p,a1r1p,&
             & a0s0p,a1s1p,cddr0,cddr1,cdds0,cdds1,cso0,cso1,tot1,tot2
  real(wp), allocatable,dimension(:) :: mesh,ucoul
  real(wp), allocatable, dimension(:,:) :: rho,tau,jsc,drho,ddrho,dtau,djsc,laprho
  real(wp), allocatable, dimension(:,:) ::uc,umr,dumr, d2umr,udd,uso,ucso
  real(wp), dimension(0:nbox,lmax,0:lmax,2,2) :: wfr,hpsi
  real(wp), allocatable, dimension(:,:,:,:) :: energies
  real(wp), allocatable :: potential(:),vocc(:,:,:,:)
  real(wp) :: j, dwfr(0:nbox,lmax,0:lmax,2,2),ddwfr(0:nbox,lmax,0:lmax,2,2),sortenergies(1:nmax,2)

  integer :: nbox, nodes, radius, lmax,nmax,njoin,itermax,nnmax
  integer :: nt,icoul,icm,npr,iq,ir,i,l,n,is,nn,np,ir2,ij2terms
  logical :: j2terms
  integer :: sortstates(1:nmax,1:3,2)

  if(ij2terms==1) j2terms = .True.
!f2py depend(nbox,lmax,nmax) :: wfr,sortstates,sortenergies,hpsi
!f2py depend(nmax) :: sortstates,sortenergies
!f2py intent(in,out):: hpsi
  if(icoul==1) iq = 1
  hbar22m = 20.73553
  nn = 28
  np = 20
  sig = 1./sig
  nnmax = lmax-2
  if (lmax <= 3) nnmax = lmax 
  nt = np+nn
  nrad = r0 * (nt)**(1._wp/3._wp)

  cmcorr = 1._wp - (1._wp/real(nt))
  if (icm==0) cmcorr = 1._wp
  spin(1) = -0.5
  spin(2) = 0.5
  njoin = 2.0/h
  nmax = nn - np

  if (nmax.ge.0) then
  nmax = nn
  else
  nmax = np
  end if
  !!
  !setting the coupling constants
  !!
  a0r0 = 3._wp/8._wp * t0
  a1r1 = - 1._wp/4._wp * t0 * ( 1._wp/2._wp + x0 )
  a0s0 = - 1._wp/4._wp * t0 * ( 1._wp/2._wp - x0 )
  a1s1 = - 1._wp/8._wp * t0
  !
  a0tau0 = 3._wp/16._wp * t1 + 1._wp/4._wp * t2 * ( 5._wp/4._wp + x2 )
  a1tau1 = - 1._wp/8._wp * t1 * ( 1._wp/2._wp + x1 ) + &
         & 1._wp/8._wp * t2 * ( 1._wp/2._wp + x2 )
  !
  a0t0 = - 1._wp/8._wp *t1* ( 1._wp/2._wp - x1 ) + &
         & 1._wp/8._wp * t2 * ( 1._wp/2._wp + x2 )
  a1t1 =  1._wp/16._wp * (t2-t1)
  !
  a0r0p = - 9._wp/64._wp * t1 + 1._wp/16._wp * t2 *( 5._wp/4._wp + x2 )
  a1r1p = 3._wp/32._wp * t1 * ( 1._wp/2._wp + x1 ) + &
         & 1._wp/32._wp * t2 * ( 1._wp/2._wp + x2 )

  a0s0p = 3._wp/32._wp * t1 * ( 1._wp/2._wp - x1 ) + &
         & 1._wp/32._wp * t2 * ( 1._wp/2._wp + x2 )
  a1s1p = 3._wp/64._wp * t1 + 1._wp/64._wp * t2
!
  cddr0 =  t3 / 16.0_wp
  cddr1 = - 1._wp/24._wp * t3 * ( 1._wp/2._wp + x3 )
  cdds0 = - 1._wp/24._wp * t3 * ( 1._wp/2._wp - x3 )
  cdds1 = - 1._wp/48._wp * t3
  !
  cso0 = - 3._wp/4._wp * w0
  cso1 = - 1._wp/4._wp * w0
  small = 1E-25_wp
  allocate(mesh(0:nbox))
  mesh = (/ (real(i)*h,i=0,nbox) /)
  mesh(0) = small
  allocate(rho(0:nbox,4),drho(0:nbox,4),ddrho(0:nbox,4),&
          tau(0:nbox,4),jsc(0:nbox,4),djsc(0:nbox,4),laprho(0:nbox,4))

  allocate(uc(0:nbox,2),umr(0:nbox,2),udd(0:nbox,2),uso(0:nbox,2),ucoul(0:nbox),ucso(0:nbox,2),dumr(0:nbox,2),d2umr(0:nbox,2))
  allocate(potential(0:nbox),vocc(lmax,0:lmax,2,2),energies(lmax,0:lmax,2,2))
  udd=0._wp
  ucso=0.0_wp
  dumr=0.0_wp
  d2umr=0.0_wp
  dwfr = 0.0
  ddwfr = 0.0
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
         call derivative_sca( nbox, wfr(1:nbox,n,l,is,iq), dwfr(1:nbox,n,l,is,iq), l, h )
         call derivative_sca( nbox, dwfr(1:nbox,n,l,is,iq), ddwfr(1:nbox,n,l,is,iq), l, h )
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

  do iq = 1,2
    do ir = 0,nbox

   !!Central Field U(r)
            uc(ir,iq) = 2*(a0r0-a1r1)*rho(ir,3) + 4*a1r1 * rho(ir,iq)  &
                           + (a0tau0-a1tau1) *tau(ir,3)+ 2 *a1tau1*tau(ir,iq) &
                           +2*( a0r0p-a1r1p )*laprho(ir,3) + 4 *a1r1p * laprho(ir,iq)
   !!Part of U(r) coming from so)
            ucso(ir,iq) = (cso0-cso1 ) *(djsc(ir,3) + 2 * jsc(ir,3)/mesh(ir) ) &
                           + 2 *cso1 * ( djsc(ir,iq) + 2 * jsc(ir,iq) / mesh(ir) )
   !!Mq(r) contributions
            umr(ir,iq) = hbar22m*cmcorr+(a0tau0-a1tau1)*rho(ir,3) + 2 * a1tau1*rho(ir,iq)
            dumr(ir,iq) = (a0tau0-a1tau1)*drho(ir,3) + 2 * a1tau1*drho(ir,iq)
            d2umr(ir,iq) = (a0tau0-a1tau1)*ddrho(ir,3) + 2 * a1tau1*ddrho(ir,iq)
   !! t3 part of U(r)
             udd(ir,iq) = ( 2 + sig ) * (cddr0-cddr1)*rho(ir,3)**(sig+1)  &
                           +2*sig*cddr1*(rho(ir,1)**2+rho(ir,2)**2)*rho(ir,3)**(sig-1.) &
                           + 4 * cddr1 * rho(ir,iq) * rho(ir,3)**sig
    !!spin-orbit part
            uso(ir,iq) = - (cso0-cso1 )*drho(ir,3)/mesh(ir) &
                            - 2 *cso1 * drho(ir,iq) / mesh(ir)
            if (j2terms) then
            uso(ir,iq) = uso(ir,iq)-(a0t0-a1t1) *jsc(ir,3) / mesh(ir) &
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
       ucoul(ir)=4.0d0*pi*e2*(tot1/mesh(ir)&
       + tot2)*h - e2*(3./pi)**(1./3.)*rho(ir,2)**(1./3.)
     end if
    end do
   end do

  
    do iq =1,2

      if (iq == 1) then
      npr = nn
      else
      npr = np
      end if
      do i = 1, npr
         if (sortenergies(i,iq) < - small) then
           n = sortstates(i,1,iq)
           l = sortstates(i,2,iq)
           is= sortstates(i,3,iq)
           j = sortstates(i,2,iq) + spin(sortstates(i,3,iq))
           if (sortstates(i,2,iq) == 0) j = 0.5
           do ir=1,nbox

            hpsi(ir,n,l,is,iq) = +umr(ir,iq)*ddwfr(ir,n,l,is,iq)+dumr(ir,iq)*dwfr(ir,n,l,is,iq) &
            + (-uc(ir,iq) -ucso(ir,iq)-udd(ir,iq) &
            -uso(ir,iq)*(j*(j+1._wp)- l*(l+1._wp) - 0.75_wp)&!/(2.0*mesh(:)) &
            -dumr(ir,iq)/mesh(ir) + (1._wp-iq)*ucoul(ir) - umr(ir,iq)*l*(l+1._wp)/mesh(ir)**2) * wfr(ir,n,l,is,iq)

            end do

          end if
      end do
    end do
  
end subroutine

subroutine derivative_sca( n, f, df, l, h )
  implicit none
  integer, parameter :: wp=kind(1.0d0)

  integer, intent(in) :: n
  real (wp), intent(in) :: f(n)
  real (wp), intent(inout) :: df(n)
  integer, intent(in) :: l
  real (wp) :: sig, h_12, h_60, h
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