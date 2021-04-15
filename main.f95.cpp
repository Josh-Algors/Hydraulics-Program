!Pipeline may have up to 10 pipes in series, and each pipe may have
!up to 100 computational nodes; this limits can be increased by
!modifying the dimension statement

!Wave velocity is adjusted, if necessary, to avoid interpolation error

!Relative valve opening vs time curve is specified at discrete points and
!values at intermediate times are determined by parabolic interpolation

!SI units are used here.

!====================Notation====================!
!   a       = wave speed (m/s)
!   ar      = pipe c/s area (m2)
!   d       = pipe diameter (m)
!   dt      = computational time interval (s)
!   dxt     = time interval for storing tau vs time curve (s)
!   f       = darcy-weisbach friction factor
!   h       = piexometric head at beginning of time interval (m)
!   hmax    = maximum piezometric head (m)
!   hmin    = minimum piezometrci head (m)
!   hp      = piezometric head at end of the time interval (m)
!   hres    = reservoir level above datum (m)
!   hs      = valve head loss for flow of qs (m)
!   iprint  = number of time interval after which condition are to be printed
!   l       = pipe length (m)
!   m       = number of points on tau vs time curve
!   n       = number of reaches into which pipe is subdivided
!   np      = number of pipes
!   nrlp    = number of reaches on last pipe
!   q       = discharge at beginning of time interval (m3/s)
!   qo      = steady state discharge (m3/s)
!   qp      = discharge at the end of time interval (m3/s)
!   qs      = valve discharge (m3/s)
!   t       = time (s)
!   tau     = relative valve opening
!   tauf    = final valve opening
!   tauo    = initial valve opening
!   tlast   = time upto which conditions are to be computed (s)
!   tv      = valve opening or closing time (s)
!   y       = stored tau values
!================================================!
implicit none
real , dimension(10,100) :: q,h,qp,hp,hmax,hmin
real , dimension(10)     :: ca,f,cf,ar,a,l,n,d,y,an,aunadj,ani
integer                  :: np,nrlp,iprint,m,i,j,k,nn,n1,np1,jp1,jm1,nm1
real                     :: qo,hres,tlast,g,tv,dxt,tauo,tauf,qs,hs,dt,t,tau,cn,cp,cv
!reading and writing the data

    g=9.81
    print *,'Enter np,nrlp,iprint,qo,hres,tlast'
    read *,np,nrlp,iprint,qo,hres,tlast
    write (*,1) np,nrlp
    write (*,2) qo,hres
    write (*,3) tlast
    1 format (8x, 'Number of pipes=',i3/8x,'Number of reaches on last pipe=',i3)
    2 format (8x,'Steady state discharge=',f6.3,' m3/s'/8x,'Reservoir level=',f6.1,' m')
    3 format (8x,'Time for which transients computed=',f6.1,' s'/)

!data for valve

    print *,'Enter m,tv,dxt,tauo,tauf,qs,hs,y'
    read *,m,tv,dxt,tauo,tauf,qs,hs,(y(i),i=1,m)
    write (*,4) m,tv
    write (*,5) dxt,hs,qs
    write (*,6) (y(i),i=1,m)
    4 format (8x, 'Number of points on tau vs time curve=',i2/8x,'Valve operation time=',f6.2,' s')
    5 format (8x, 'Time interval for storing tau curve=',f6.3' s'/8x,'Valve loss=',f6.2,' m for qs=',f6.3,' m3/s')
    6 format (8x, 'Stored tau values:'/8x,15f8.3/)

!data for pipes

    print *,'Enter Length of Pipe,Diameter of Pipe,Wave Velocity,Friction Factor'
    read *,(l(i),d(i),a(i),f(i),i=1,np)
    write (*,7)
    7 format (/8x, 'Pipe no',5x, 'Length',5x , 'Dia',5x, 'Wave vel.',5x, 'Fric factor'/21x,'(m)',7x,'(m)',7x, '(m/s)')
    write (*,8)(i,l(i),d(i),a(i),f(i),i=1,np)
    8 format (10x, i3,6x, f7.1,3x, f5.2,5x, f7.1,11x, f5.3)
    dt=l(np)/(nrlp*a(np))
    write (*,9)
    9 format(/8x,'Pipe no',5x,'Adjusted wave velocity'/27x,'(m/s)')

! calculation of pipe constants
    do i=1,np
        ar(i)=0.7854*d(i)*d(i)
        aunadj(i)=a(i)
        an(i)=l(i)/(dt*a(i))
        n(i)=an(i)
        ani(i)=n(i)
        if (an(i)-ani(i) >= 0.5) n(i)=n(i)+1
        a(i)=l(i)/(dt*n(i))
        write(*,10) i,a(i)
        10 format (10x, i3,12x, f7.1)
        ca(i)=g*ar(i)/a(i)
        cf(i)=f(i)*dt/(2.0*d(i)*ar(i))
        f(i)=f(i)*l(i)/(2.0*g*d(i)*n(i)*ar(i)*ar(i))
        continue
    end do

!Calculation of steady state

    h(1,1)=hres
    do i=1,np
        nn=n(i)+1
        do j=1,nn
            h(i,j)=h(i,1)-(j-1)*f(i)*qo*qo
            q(i,j)=qo
            continue
        end do
        h(i+1,1)=h(i,nn)
        continue
    end do
    nn=n(np)+1
    hs=h(np,nn)
    qs=qo
    do i=1,np
        nn=n(i)+1
        do j=1,nn
            hmax(i,j)=h(i,j)
            hmin(i,j)=h(i,j)
            continue
        end do
        continue
    end do
    t=0.0
    tau=tauo
    write(*,11)
    11 format (/8x,'Time',2x,'Tau',2x,'Pipe',7x,'Head (m)',7x,'Disch. ','(m3/s)'/20x,'No',5x,'(1)',5x,'(n+1)',5x,'(1)',5x,'(n+1)'/)
    128 k=0
        i=1
        nn=n(i)+1
        write(*,12) t,tau,i,h(i,1),h(i,nn),q(i,1),q(i,nn)
    12 format (f12.1,f6.3,i4,2f9.2,f9.3,f10.3)
    if (np == 1) go to 14
    do i=2,np
        nn=n(i)+1
        write(*,13) i,h(i,1),h(i,nn),q(i,1),q(i,nn)
        13 format (20x,i2,2f9.2,f9.3,f10.3)
        continue
    end do
    14 t=t+dt
       k=k+1

    if ( t>tlast) go to 19

! Upstream reservoir

    hp(1,1)=hres
    cn=q(1,2)-h(1,2)*ca(1) - cf(1)*q(1,1)*ABS(q(1,1))
    qp(1,1)= cn+ca(1)*hres

!interior points
    do i=1,np
        nn=n(i)
        do j=2,nn
            jp1=j+1
            jm1=j-1
            cn=q(i,jp1)-ca(i)*h(i,jp1)-cf(i)*q(i,jp1)*ABS(q(i,jp1))
            cp=q(i,jm1)+ca(i)*h(i,jm1)-cf(i)*q(i,jm1)*ABS(q(i,jm1))
            qp(i,j)=0.5*(cp+cn)
            hp(i,j)=(cp-qp(i,j))/ca(i)
            continue
        end do
        continue
    end do

!Series Connection
    np1 = np-1
    if(np==1) go to 178
    do i=1,np1
        n1 = n(i)
        nn = n(i) + 1
        cn = q(i+1,2) - ca(i+1) * h(i+1,2) - cf(i+1) * q(i+1,2) * ABS(q(i+1,2))
        cp = q(i,n1) + ca(i) * h(i,n1) - cf(i) * q(i,n1) * ABS(q(i,n1))
        hp(i,nn) = (cp-cn)/(ca(i)+ca(i+1))
        hp(i+1,1) = hp(i,nn)
        qp(i,nn) = cp-ca(i)*hp(i,nn)
        qp(i+1,1) = cn + ca(i+1) *hp(i+1,1)
        continue
    end do


!valve at downstream end

    178 nn=n(np)+1
    nm1=n(np)
    cp=q(np,nm1)+ca(np)*h(np,nm1)-cf(np)*q(np,nm1)*ABS(q(np,nm1))
    if (t >= tv) go to 15
    call parab(t,dxt,y,tau)
    go to 16
    15 tau=tauf
    if (tau <= 0.0) go to 17
    16 cv=(qs*tau)**2/(hs*ca(np))
       qp(np,nn)=0.5*(-cv+sqrt(cv*cv+4.*cp*cv))
       hp(np,nn)=(cp-qp(np,nn))/ca(np)
       go to 18
    17 qp(np,nn)=0.0
       hp(np,nn)=cp/ca(np)

! storing variables for next time step

    18 do i=1,np
            nn=n(i)+1
            do j=1,nn
                q(i,j)=qp(i,j)
                h(i,j)=hp(i,j)
                if (h(i,j) > hmax(i,j)) hmax(i,j)=h(i,j)
                if (h(i,j) < hmin(i,j)) hmin(i,j)=h(i,j)
                continue
            end do
            continue
       end do

    if (k==iprint) go to 128
    go to  14
    19 write(*,20)
    20 format (/8x,'Pipe no.',3x,'Section no.',3x,'Max press.',3x,'Min press.'/)

    do i=1,np
        nn=n(i)+1
        do j=1,nn
        write(*,21) i,j,hmax(i,j),hmin(i,j)
        21 format (9x,i2,13x,i2,2f13.2)
        end do
    end do
    stop
end
subroutine parab(t,dxt,y,tau)
implicit none
real ,dimension(10) ::y
real                ::t,dxt,tau,r
integer             ::i
    i=t/dxt
    r=(t-i*dxt)/dxt

    if (i==0) r=r-1
      i=i+1
    if (i<2) i=2
      tau=y(i)+0.5*r*(y(i+1)-y(i-1)+r*(y(i+1)+y(i-1)-2.0*y(i)))
end subroutine parab
