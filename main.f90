PROGRAM TRIDIAG
implicit none
integer rc,i,number_of_iterations,j,l,nmat, timestep, mm, kr
real*8, allocatable, dimension(:) :: a, b, c, r, u, t, cap, k, x, w, diag, spdiag, sbdiag, xout, rhs, wgamma, told, tprev
real*8 :: rho=917 !density
real*8 :: c0=2106 !capacity of fresh ice
real*8 :: L0=334000 !latent heat of fusion of fresh ice
real*8 :: mu=0.054 !ratio between the freezing temperature and salinity of brine
real*8 :: k0=2.03 !the conductivity of fresh ice
real*8 :: beta=0.13 !empirical constant
real*8 :: S=0.0 !solenity of ice
real*8 :: dt=1*3600 !time step in seconds
real*8 :: I0=0 !1 !10!2*25*(1-0.65)*(1-0)*0.7/100 !radiation
real*8 :: kappa0=1.4 !radiation
real*8 :: cw=3985 !capacity of water
real*8 :: emissivity=0.95 !emissivity
real*8 :: flw=0.3 !incoming longwave radiation
real*8 :: Tw=0.0 !temperature of water in C
!real*8 :: Ta=-10 !temperature of air in C
real*8 :: vonkar=0.4 ! von Karman constant
real*8 :: rhoa = 1.3    ! air density (kg/m^3)
real*8 :: zref   = 10  ! reference height for stability (m)
real*8 :: iceruf   = 0.0005 ! ice surface roughness (m)
logical :: conv
real*8 x1, x2, h, Fbot, Fcb, q, Tf, dh, h0, wbeta, totalh, totalh0, F0, dfsurf_dT, Ta, rdn, ustar, rh, shcoef
integer :: n=7 !number of points (>2)
open(2,file='lolminus.txt',status='unknown',decimal='comma')
open(3,file='haha5.txt',status='unknown',decimal='comma')
open(4,file='t(2)check.txt',status='unknown',decimal='comma')
allocate(a(n+1))
allocate(b(n+1))
allocate(c(n+1))
allocate(r(n+1))
allocate(u(n+1))
allocate(t(n+2))
allocate(told(n+2))
allocate(tprev(n+2))
allocate(cap(n+2))
allocate(k(n+2))
allocate(x(n+2))
allocate(w(n+1))
allocate(diag(n+1))
allocate(spdiag(n+1))
allocate(sbdiag(n+1))
allocate(xout(n+1))
allocate(rhs(n+1))
allocate(wgamma(n+1))
rdn  = vonkar/log(zref/iceruf)
ustar = rdn*1
!real*8 :: psimhs = -(0.7*hol(ij) + 0.75*(hol(ij)-14.3) * exp(-0.35*hol(ij)) + 10.7)
!real*8 :: psixh(ij)  = psimhs*stable(ij) + (1 - stable(ij))*psixhu(xqq)
rh = rdn !/ (1+rdn/vonkar*(log(1)-psixh(ij)))
shcoef=1 !+rhoa*1*ustar*rh !transfer coefficient for sensible heat
number_of_iterations=4
totalh=2
totalh0=totalh
Tf=-mu*S
do i=2,n
    t(i)=-(i+1)/2
    w(i)=0
end do
t(1)=0
w(1)=0
write(*,*) Tf, "Tf"
write(*,*) shcoef
do timestep=1,50000
    if (t(n-1) >= 0 ) then
        write(*,*) t(n-1)
        exit
    end if
    if (totalh<=0) then
        exit
    endif
    Ta=-10*sin(2*3.14*(timestep-1200/1.0)/(8760/1.0) + 3.14/2)-5
    write(*,*) timestep
  !      write(2,*) "============================================================================="
    write(2,*) "STEP NO", timestep
    write(2,*) "Ta", Ta
  !  write(2,*) "============================================================================="
    h0=totalh0/(n-1)
    do i=1,n
        write(2,*) t(i), "t0(", i, ")"
        told(i)=t(i)
        tprev(i)=t(i)
    end do
    F0=shcoef*(Ta-t(n))+I0
    dfsurf_dT=-shcoef
    Fbot=-0.006*1020*cw*0.0005*(Tw-t(1))
    if (t(n) >= Tf) then
        t(n) = Tf-0.000000000001
    endif
    conv = .false.
    do while (.not. conv)
        do i=1,n
            tprev(i)=t(i)
        end do
        h=totalh/(n-1)
   ! k(1)=k0+beta*S/t(n)
        k(1)=k0
     !   write(2,*) "F0", F0, "Fct", k(1)*(t(n)-t(n-1))/h
      !  if (t(n) >= Tf) then
     !   if (t(n) >= Tf-0.000001) then
     !       t(n) = Tf
     !   endif
        if (t(n) >= Tf) then
          !  write(*,*) t(n), "!!!"
            write(*,*) "0 0 0 0 0 0 0"
            w(n)=k(1)*(t(n)-t(n-1))/h - F0
            F0=k(1)*(t(n)-t(n-1))/h
        !       q=c0*(t(n)-Tf)-L0*(1-Tf/t(n))+cw*Tf
            q=c0*t(n)-L0
            w(n)=w(n)/(rho*q)
            nmat=n-2
            if (nmat == 1) then
             !   diag(1)=2*k0/h0+2*h*rho*c0/(dt*3)
             !   rhs(1)=h*rho*c0*(2*told(2)/3)+k0*told(1)/h0 +k0*told(3)/h0 &
             !   +rho*w(2)*(c0*(told(n)+told(2))/2-c0*(told(1)+told(2))/2)
            else
                diag(1)=2*k0/h0+2*h*rho*c0/(dt*3)
                spdiag(1)=-k0/h0+h*rho*c0/(dt*6)
                rhs(1)=h*rho*c0*(2*told(2)/3+told(3)/6)/dt+k0*told(1)/h0 &
                +rho*w(2)*(c0*(told(n)+told(2))/2-c0*(told(1)+told(2))/2) &
                - I0*((1-exp(kappa0*h*1))/2 + (exp(kappa0*h*1)-exp(kappa0*h*2))/2)
                do kr = 2,nmat-1
                    spdiag(kr)=-k0/h0+h*rho*c0/(dt*6)
                    diag(kr)=2*k0/h0+2*h*rho*c0/(dt*3)
                    sbdiag(kr)=-k0/h0+h*rho*c0/(dt*6)
                    rhs(kr)=h*rho*c0*(told(kr)/6+2*told(kr+1)/3+told(kr+2)/6)/dt &
                    +rho*w(kr+1)*(c0*(told(kr+2)+told(kr+1))/2-c0*(told(kr) +told(kr+1))/2)&
                     - I0*((exp(kappa0*h*(kr-1)) &
                    -exp(kappa0*h*kr))/2 + (exp(kappa0*h*kr)-exp(kappa0*h*(kr+1)))/2)
                enddo
                diag(nmat)=2*k0/h0+2*h*rho*c0/(dt*3)
                sbdiag(nmat)=-k0/h0+h*rho*c0/(dt*6)
                rhs(nmat)=h*rho*c0*(2*told(nmat+1)/3+told(nmat)/6)/dt+k0*told(nmat+1)/h0 &
                +rho*w(nmat+1)*(c0*(told(nmat+2)+ told(nmat+1))/2-c0*(told(nmat)+told(nmat+1))/2) &
                -I0*((exp(kappa0*h*(nmat-1)) -exp(kappa0*h*nmat))/2 + (exp(kappa0*h*nmat)-exp(kappa0*h*(nmat+1)))/2)
            endif
            wbeta=diag(1)
            xout(1)=rhs(1)/wbeta
            do l = 2, nmat
                wgamma(l) = spdiag(l-1) / wbeta
                wbeta = diag(l) - sbdiag(l)*wgamma(l)
                xout(l) = (rhs(l) - sbdiag(l)*xout(l-1)) / wbeta
            enddo
            do l = nmat-1, 1, -1
                xout(l) = xout(l) - wgamma(l+1)*xout(l+1)
            enddo
            do l = 1, nmat
                t(l+1) = xout(l)
            enddo
        else

            nmat=n-1
            diag(1)=2*k0/h0+2*h*rho*c0/(dt*3)
            spdiag(1)=-k0/h0+h*rho*c0/(dt*6)
            rhs(1)=h*rho*c0*(2*told(2)/3+told(3)/6)/dt+k0*told(1)/h0 &
            +rho*w(2)*(c0*(told(n)+told(2))/2-c0*(told(1)+told(2))/2) &
            - I0*((1-exp(kappa0*h*1))/2 + (exp(kappa0*h*1)-exp(kappa0*h*2))/2)
            do kr = 2,nmat-1
                spdiag(kr)=-k0/h0+h*rho*c0/(dt*6)
                diag(kr)=2*k0/h0+2*h*rho*c0/(dt*3)
                sbdiag(kr)=-k0/h0+h*rho*c0/(dt*6)
                rhs(kr)=h*rho*c0*(told(kr)/6+2*told(kr+1)/3+told(kr+2)/6)/dt &
                +rho*w(kr+1)*(c0*(told(kr+2)+told(kr+1))/2-c0*(told(kr)+told(kr+1))/2) &
                -I0*((exp(kappa0*h*(kr-1))-exp(kappa0*h*kr))/2 + (exp(kappa0*h*kr)-exp(kappa0*h*(kr+1)))/2)
            enddo
            diag(nmat)=-dfsurf_dT+k0/h0 +h*rho*c0/(3*dt)
            sbdiag(nmat)=-k0/h0 +h*rho*c0/(6*dt)
            rhs(nmat)=-dfsurf_dT*t(n)+F0 +rho*w(n)*(c0*(told(n)+told(n-1))/2) &
            +h*rho*c0*(told(n-1)/6+told(n)/3)/dt
          !  write(*,*) diag(nmat)
          !  write(*,*) sbdiag(nmat)
          !  write(*,*) rhs(nmat)
            wbeta=diag(1)
            xout(1)=rhs(1)/wbeta
            do l = 2, nmat
                wgamma(l) = spdiag(l-1) / wbeta
                wbeta = diag(l) - sbdiag(l)*wgamma(l)
                xout(l) = (rhs(l) - sbdiag(l)*xout(l-1)) / wbeta
            enddo
            do l = nmat-1, 1, -1
                xout(l) = xout(l) - wgamma(l+1)*xout(l+1)
            enddo
            do l = 1, nmat
                t(l+1) = xout(l)
            enddo
            if (t(n) >= Tf) then
                t(n) = Tf
           ! else
            !j    F0=k(1)*(t(n)-t(n-1))/h
            w(n)=0
            endif
        endif
    !k/h0*(-t(1)+2*t(2)-t(3))+h*rho*c0*(t(1)/6+2*t(2)/3+t(3)/6)=rho*w(2)*(c0*(t(n)+t(2))/2-c0*(t(1)+t(2))/2)+h*rho*c0*(t(1)/6+2*t(2)/3+t(3)/6)
    !k0*(t(n)-t(2))/h=shcoef*(Ta-t(n))
    !    write(2,*) (k(1)*(t(2)-t(1))/h), "(k(1)*(t(2)-t(1))/h)"
  !      write(2,*) Fbot, "Fbot"
        w(1)=k(1)*(t(2)-t(1))/h - Fbot
        Fbot=k(1)*(t(2)-t(1))/h
    !q=c0*(t(1)-Tf)-L0*(1-Tf/t(1))+cw*Tf
        q=c0*t(1)-L0
        w(1)=w(1)/(rho*q)
        totalh=(w(1)-w(n))*dt + totalh
        do l = 2, n-1
            w(n-l+1)=(l-1)*(w(1)-w(n))/(n-1)+w(n)
        enddo
        do i=1,n
            if (ABS(t(i)-tprev(i))>1.e-12) then
                conv = .false.
                exit
            else
                conv = .true.
            end if
        end do
    enddo
    totalh0 = totalh
    write(3,*) timestep*150.0/8760, totalh
    if (timestep > 49990000) then
    write(4,*) timestep
    do mm = 1,n
            write(4,*) t(mm), "t(", mm, ")"
        enddo
        do mm = 1,n
            write(4,*) w(mm), "w(", mm, ")"
        enddo
        write(4,*) totalh, "totalh"
    endif
enddo
close(2)
close(3)
close(4)
  end

