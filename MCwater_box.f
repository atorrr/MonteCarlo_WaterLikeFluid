c This program calculates the structure of a 3D system composed of spheres that interact with a water-like potential.
c The diameter of the particles is used as unit length.
c Ramon Castaneda Priego. 04.03.2019
      implicit integer*4(i-n),real*8(a-h,o-z)
      parameter(mp=10000,mr=2**9,nvq=20)
	  dimension x(mp),y(mp),z(mp)
      dimension r(mr),g(mr),q(mr),s(mr)
	  common/box/boxl,rc,np,ky,kz,rcz,rcy
      common/pot/deps,dl,dw,dro,drc
	  common/parameters/diam,rho
      common/sf/qx(mr,nvq),qy(mr,nvq),qz(mr,nvq)
	  pi=4.d0*datan(1.d0)
c number of particles
	  np=15**3
c filling fraction
	  phi=0.364353
c diameter of the particles
      diam=1.d0
c box length
	  ky=1.d0
	  kz=14.d0
	  boxl=(pi*np/(6.*phi*ky*kz))**(1.d0/3.d0)
	  rc=boxl/2.d0
	  rcy=rc*ky
	  rcz=rc*kz
      dr=rc/mr
      dq=pi/rc
c      print*,'rc= ',rc
c	  print*,'dr = ',dr
c      print*,'dq = ',dq
      do i=1,mr
         r(i)=(i-1)*dr
         q(i)=(i-1)*dq
      enddo
      open(121,file='qvectors.dat',status='unknown')
	  ncq=0
	  do i=1,mr
	     do j=1,nvq
            dt=pi*ranf(iseed)
	        dphi=2*pi*ranf(iseed)
	        qx(i,j)=q(i)*dcos(dphi)*dsin(dt)
            qy(i,j)=q(i)*dsin(dphi)*dsin(dt)
	        qz(i,j)=q(i)*dcos(dt)
	        ncq=ncq+1
	        qs=dsqrt(qx(i,j)*qx(i,j)+qy(i,j)*qy(i,j)+qz(i,j)*qz(i,j))
	        write(121,100)dfloat(ncq),qs
	     enddo
	  enddo
 	  close(121)
c reduced density
	  rho=6.d0*phi/pi
      d=rho**(-1./3.)
      print*,'Number density: ',rho
	  print*,'Mean interparticle distance: ',d
c initial configuration
	  call iniconfig(x,y,z)
c write the initial configuration
	  open(10,file='iniconf.dat',status='unknown')
	  do i=1,np
	     write(10,100)x(i),y(i),z(i)
      enddo
	  close(10)
c potential parameters
      deps=1.d0/3.d0
      dl=1.717d0
      dw=5.d0
      dro=1.5
      drc=2.5
      call potential(d,ud)
      print*,d,ud

      call energy(x,y,z,ener)
      print*,'Energy per particle of the initial configuration:',ener/np

c MC cycle to thermalize the system
      open(30,file='energy.dat',status='unknown')
	  del=0.1d0
	  nattemp=0
	  nacc=1
	  iseed=123456789
	  do i=1,1000000
	     call mcmove(x,y,z,ener,nattemp,nacc,del,iseed)
	     call adjust(nattemp,nacc,del)
	     if (mod(i,100) .eq. 0) write(30,100)i*1.d0,ener/np
	     if (mod(i,100000) .eq. 0) then
            print*,i,del,ener/np
         endif
      enddo
	  print*,'The system has thermalized'
c write the final configuration and the energy
	  open(20,file='finalconf.dat',status='unknown')
	  do i=1,np
	     write(20,100)x(i),y(i),z(i)
	  enddo
	  close(20)
	  close(30)

c MC cycle to calculate the g(r)
      nacco=nacc
	  do i=1,mr
	     g(i)=0.d0
	  enddo
      ng=0
      naveg=0
c Averages
	  do i=1,1000000
	     call ave(x,y,z,g,s,ener,nattemp,nacc,ng,naveg,del,dr,iseed)
         call adjust(nattemp,nacc,del)
         if (mod(i,100000) .eq. 0) print*,i,'Averaging'
	  enddo
	  nav=nacc-nacco
	  print*,'Average number for energy: ',nav
      print*,'Average number for g(r): ',naveg

c This is the radial distribution function
	  open(50,file='gr.dat',status='unknown')
	  do i=2,mr
	     dv=(4.d0*pi*r(i)**2.*dr)*rho
	     g(i)=g(i)/(np*naveg*dv)
	     write(50,200)r(i),g(i)
	  enddo
	  close(50)

c This is the structure factor from the definition
      open(51,file='sq.dat',status='unknown')
	  do i=3,mr
	     s(i)=s(i)/naveg
         if (q(i) .lt. 30.d0) then
	        write(51,100)q(i),s(i)
         endif
	  enddo
	  close(51)

100   format(3f15.7)
200   format(2f15.7)
	  end
c
c This subroutine calculates the initial configuration in 3D.		 
      subroutine iniconfig(x,y,z)
      implicit integer*4(i-n),real*8(a-h,o-z)
      parameter(mp=10000)
      dimension x(mp),y(mp),z(mp)
	  common/box/boxl,rc,np,ky,kz,rcz,rcy
	  common/parameters/diam,rho
	  d=rho**(-1.d0/3.d0)
      x(1)=-rc+d/2.d0
      y(1)=-rcy+d/2.d0
	  z(1)=-rcz+d/2.d0
      do i=1,np-1
         x(i+1)=x(i)+d
         y(i+1)=y(i)
	     z(i+1)=z(i)
         if (x(i+1) .gt. rc) then
            x(i+1)=-rc+d/2.
            y(i+1)=y(i+1)+d
	        z(i+1)=z(i)
	        if (y(i+1) .gt. rcy) then
	           x(i+1)=-rc+d/2.
	           y(i+1)=-rcy+d/2.
	           z(i+1)=z(i+1)+d
	        endif
	     endif
	   enddo
       return
       end
c
c This subroutine calculates the energy of a given configuration
      subroutine energy(x,y,z,ener)
	  implicit integer*4(i-n),real*8(a-h,o-z)
	  parameter(mp=10000)
      dimension x(mp),y(mp),z(mp)
	  common/box/boxl,rc,np,ky,kz,rcz,rcy
	  common/parameters/diam,rho
	  pi=4.d0*datan(1.d0)
	  ener=0.d0
      do i=1,np-1
         do j=i+1,np
            xij=x(j)-x(i)
            yij=y(j)-y(i)
	        zij=z(j)-z(i)
            xij=xij-boxl*dnint(xij/boxl)
            yij=yij-boxl*dnint(yij/boxl)
	        zij=zij-boxl*dnint(zij/boxl)
            rij2=xij*xij+yij*yij+zij*zij
            rij=dsqrt(rij2)
            if (rij .lt. rc) then
               call potential(rij,uij)
               ener=ener+uij
	        endif
         enddo      
      enddo
	  return
	  end
c
c This subroutine calculates the difference in energy when a particle is displaced
      subroutine denergy(x,y,z,no,dener)
	  implicit integer*4(i-n),real*8(a-h,o-z)
	  parameter(mp=10000)
      dimension x(mp),y(mp),z(mp)
	  common/box/boxl,rc,np,ky,kz,rcz,rcy
	  common/pot/deps,dl,dw,dro,drc
	  common/parameters/diam,rho
	  pi=4.d0*datan(1.d0)
	  dener=0.d0
      do i=1,np
         if (i .ne. no) then
            xij=x(no)-x(i)
            yij=y(no)-y(i)
	        zij=z(no)-z(i)
            xij=xij-boxl*dnint(xij/boxl)
            yij=yij-boxl*dnint(yij/boxl)
	        zij=zij-boxl*dnint(zij/boxl)
            rij2=xij*xij+yij*yij+zij*zij
            rij=dsqrt(rij2)
            if (rij .lt. rc) then
               dB=dl/(dl-1.d0)
               if (rij .lt. dB) then
                  call potential(rij,uij)
               else
                  uij=0.d0
               endif
               dener=dener+uij
	        endif
         endif
      enddo
	  return
	  end
c
c This subroutine calculates the pair between particles i & j
      subroutine potential(rij,uij)
	  implicit integer*4(i-n),real*8(a-h,o-z)
	  common/pot/deps,dl,dw,dro,drc
      if (rij .lt. drc) then
	     uij=4.d0*deps*(rij**(-12)-rij**(-6))
         uij=uij-dl*deps*dexp(-(dw*(rij-dro))**2)
         uij=uij-4.d0*deps*(drc**(-12)-drc**(-6))
      else
         uij=0.d0
      endif
	  return
	  end
c
c This subroutine displace the system to a new configuration
	  subroutine mcmove(x,y,z,ener,nattemp,nacc,del,iseed)
	  implicit integer*4(i-n),real*8(a-h,o-z)
	  parameter(mp=10000)
      dimension x(mp),y(mp),z(mp)
	  common/box/boxl,rc,np,ky,kz,rcz,rcy
	  nattemp=nattemp+1
	  no=int(ranf(iseed)*np)+1
	  xo=x(no)
	  yo=y(no)
	  zo=z(no)
	  call denergy(x,y,z,no,enero)
	  x(no)=x(no)+(ranf(iseed)-0.5d0)*del
	  y(no)=y(no)+(ranf(iseed)-0.5d0)*del
	  z(no)=z(no)+(ranf(iseed)-0.5d0)*del
c periodic boundary conditions
	  x(no)=x(no)-boxl*dnint(x(no)/boxl)
	  y(no)=y(no)-boxl*dnint(y(no)/boxl)
	  z(no)=z(no)-boxl*dnint(z(no)/boxl)
	  call denergy(x,y,z,no,enern)
      dener=enern-enero
	  if (ranf(iseed) .lt. dexp(-dener)) then
	     ener=ener+dener
	     nacc=nacc+1
	  else
	     x(no)=xo
	     y(no)=yo
	     z(no)=zo
	  endif
	  return
	  end
c
c This subroutine calculates the g(r)
	  subroutine ave(x,y,z,g,s,ener,nat,nacc,ng,nvg,del,dr,iseed)
	implicit integer*4(i-n),real*8(a-h,o-z)
	  parameter(mp=10000,mr=2**9,nvq=20)
	dimension x(mp),y(mp),z(mp),g(mr),s(mr)
	  common/box/boxl,rc,np,ky,kz,rcz,rcy
	  pi=4.d0*datan(1.d0)
	  nat=nat+1
	  no=int(ranf(iseed)*np)+1
	  xo=x(no)
	  yo=y(no)
	  zo=z(no)
	  call denergy(x,y,z,no,enero)
	  x(no)=x(no)+(ranf(iseed)-0.5d0)*del
	  y(no)=y(no)+(ranf(iseed)-0.5d0)*del
	  z(no)=z(no)+(ranf(iseed)-0.5d0)*del
c periodic boundary conditions
	  x(no)=x(no)-boxl*dnint(x(no)/boxl)
	  y(no)=y(no)-boxl*dnint(y(no)/boxl)
	  z(no)=z(no)-boxl*dnint(z(no)/boxl)
	  call denergy(x,y,z,no,enern)
      dener=enern-enero
	  if (ranf(iseed) .lt. dexp(-dener)) then
         ener=ener+dener
         nacc=nacc+1
         ng=ng+1
c calculating the g(r)
         if (mod(ng,np) .eq. 0) then
            nvg=nvg+1
            call rdf(x,y,z,dr,g)
c            call sq(x,y,z,s)
         endif
	  else
	     x(no)=xo
	     y(no)=yo
	     z(no)=zo
	  endif
	  return
	  end
c
c This subroutine calculates the radial distribution function
      subroutine rdf(x,y,z,dr,g)
	  implicit integer*4(i-n),real*8(a-h,o-z)
	  parameter(mp=10000,mr=2**9)
      dimension x(mp),y(mp),z(mp),g(mr)
	  common/box/boxl,rc,np,ky,kz,rcz,rcy
	  do i=1,np-1
         do j=i+1,np
	        xij=x(j)-x(i)
            yij=y(j)-y(i)
	        zij=z(j)-z(i)
            xij=xij-boxl*dnint(xij/boxl)
            yij=yij-boxl*dnint(yij/boxl)
	        zij=zij-boxl*dnint(zij/boxl)
            rij2=xij*xij+yij*yij+zij*zij
            rij=dsqrt(rij2)
	        nbin=int(rij/dr)+1
	        if (nbin .le. mr) then
	           g(nbin)=g(nbin)+2.
	        endif
	     enddo
      enddo
	  return
	  end
c
c This subroutine calculates the S(q)
	  subroutine sq(x,y,z,s)
      implicit integer*4(i-n),real*8(a-h,o-z)
	  parameter(mp=10000,mr=2**9,nvq=20)
      dimension x(mp),y(mp),z(mp),s(mr)
	  dimension auxc(nvq),auxs(nvq)
	  common/box/boxl,rc,np,ky,kz,rcz,rcy
	  common/sf/qx(mr,nvq),qy(mr,nvq),qz(mr,nvq)
      do i=2,mr
         sum=0.d0
	     parti=0.d0
	     do kq=1,nvq
	        auxc(kq)=0.d0
	        auxs(kq)=0.d0
	     enddo
	     do k=1,np
	        xaux=x(k)-boxl*dnint(x(k)/boxl)
	        yaux=y(k)-boxl*dnint(y(k)/boxl)
	        zaux=z(k)-boxl*dnint(z(k)/boxl)
            rij2=xaux*xaux+yaux*yaux+zaux*zaux
	        rij=dsqrt(rij2)
	        if (rij .lt. rc) then
               parti=parti+1.d0
	           do j=1,nvq
		          arg=qx(i,j)*xaux+qy(i,j)*yaux+qz(i,j)*zaux
	              auxc(j)=auxc(j)+dcos(arg)
		          auxs(j)=auxs(j)+dsin(arg)
		       enddo
	        endif
         enddo
	     do kq=1,nvq
	        sum=sum+auxc(kq)*auxc(kq)+auxs(kq)*auxs(kq)
	     enddo
	     auxsq=sum/(nvq*parti)
	     s(i)=s(i)+auxsq
	  enddo
	  return
	  end

c Random generator algorithm
	function ranf(idum)
	implicit integer*4(i-n),real*8(a-h,o-z)
c      integer idum,IA,IM,IQ,IR,MASK
c      real ran0,AM
      parameter (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
     *MASK=123459876)
c      integer k
      idum=ieor(idum,MASK)
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      ranf=AM*idum
      idum=ieor(idum,MASK)
      return
      end
c
c This subroutine adjusts the displacement of particles
	  subroutine adjust(nattemp,nacc,dr)
	  implicit integer*4(i-n),real*8(a-h,o-z)
	  if (mod(nattemp,nacc) .eq. 0) then
	     ratio=real(nacc)/real(nattemp)
	     if (ratio .gt. 0.5) then
	        dr=dr*1.05
	     else
	        dr=dr*0.95
	     endif
	  endif
      return
      end
