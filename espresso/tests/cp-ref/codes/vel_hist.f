        Program Velocity_histogram

        character*150 filename
        character*2 atom
        integer startt,nfr,dt,grid,bin
        real ts,vac0,vac,currt,com(3),gridv,temp,mass,velnorm,mb
        real mb_pre,two_kt,p,meangrid
        real,allocatable :: xyz(:,:),xyz_old(:,:),vel(:,:)
        integer,allocatable :: velhist(:)

        p=3.1415926535
        write(*,*) 'Calculates normalized histogram of velocities for'
        write(*,*) 'a single species of atoms in a simulation cell'
        write(*,*) 'and compares this to the ideal Maxwell-Boltzmann'
        write(*,*) 'distribution'
        write(*,*) 'Output: vel_hist.dat (vel, counts, counts_MB)'
        write(*,*) '[vel in atomic units, others are unitless]'
        write(*,*) ''
        write(*,*) 'Which file? '
        read(*,*) filename
        open(unit=10,file=filename,status='old')
        write(*,*) 'Start analysis from timestep #?'
        read(*,*) startt
        startt=startt+1 !need at least 2 timesteps
        !write(*,*) 'End analysis at timestep #? (0 = end of file)'
        !read(*,*) nfr
        !if(nfr.eq.0) nfr=1e5
        nfr=1e5
        write(*,*) 'Simulation timestep [in au]? '
        read(*,*) ts
        write(*,*) 'Simulation temperature [in K]?'
        read(*,*) temp
        write(*,*) 'Ion mass [in amu]?'
        read(*,*) mass
        mass=mass*1822.8884854 ! convert to atomic units
        write(*,*) 'Velocity bin spacing [in au]?'
        read(*,*) gridv
        ivel=0
        do ifr=1,nfr
         read(10,*,end=10) nion
         read(10,*,end=10) filename
         if(ifr.eq.1) then
          allocate(xyz(3,nion))
          allocate(xyz_old(3,nion))
          allocate(vel(nion,nfr-startt+1))
         endif
         com(:)=0.0
         do iat=1,nion
          read(10,*,end=10) atom,(xyz(ixyz,iat),ixyz=1,3)
          com(:)=com(:)+xyz(:,iat)
         enddo
         com(:)=com(:)/(1.0*nion)
         do iat=1,nion
          xyz(:,iat)=xyz(:,iat)-com(:)
         enddo
         if(ifr.gt.startt) then
          ivel=ivel+1
          do iat=1,nion
           vel(iat,ivel)=sqrt((xyz(1,iat)-xyz_old(1,iat))**2
     &                       +(xyz(2,iat)-xyz_old(2,iat))**2
     &                       +(xyz(3,iat)-xyz_old(3,iat))**2)
     &                       /(ts*0.529177249) !converts to atomic units
          enddo
         endif
         xyz_old=xyz
        enddo
10      continue
        close(10)
        deallocate(xyz)
        deallocate(xyz_old)
        nfr=ivel

        write(*,*) 'vel max', maxval(vel(:,1:nfr))
        grid=int(maxval(vel(:,1:nfr))/gridv)+1
        allocate(velhist(grid))
        velhist=0.0
        do ifr=1,nfr
         do iat=1,nion
          bin=int(vel(iat,ifr)/gridv)+1
          if(bin.le.grid) then
           velhist(bin)=velhist(bin)+1
          endif
         enddo
        enddo

        open(unit=15,file='vel_hist.dat',status='unknown')
        two_kt=2.0*3.1668153e-6*temp
! MB prefactor 4pi*(m/2pikT)^3/2
        mb_pre=4*p*(mass/(two_kt*p))**1.5
        write(*,*) 'mb_pre',mb_pre
        write(*,*) 'two_kt',two_kt
        write(*,*) '1/2 mv2',0.5*mass*(gridv*(grid/2))**2
        write(*,*) 'exp',exp((-1.0*mass*(gridv*(grid/2))**2)/two_kt) 
        do igrid=1,grid
         meangrid=gridv*(igrid-0.5)
! Calculate MB factor for meangrid
! = prefactor * v^2 * e^(-mv^2/2kT)
         mb=mb_pre*meangrid**2*exp((-1.0*mass*meangrid**2)/two_kt) 
         !write(15,'(f14.8,2f11.6)') meangrid,
       !          (1.0*velhist(igrid))/(1.0*nion*nfr), mb
         write(15,*) meangrid,
     &           (1.0*velhist(igrid))/(1.0*nion*nfr), mb*gridv
        enddo

        deallocate(vel)
        deallocate(velhist)
        close(15)
        end
