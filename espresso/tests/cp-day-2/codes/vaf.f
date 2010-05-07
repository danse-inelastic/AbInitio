        Program Velocity_autocorrelation

        character*75 filename
        character*2 atom
        integer startt,nfr,dt
        real ts,cutt,vac0,vac,currt,com(3)
        real,allocatable :: xyz(:,:),xyz_old(:,:),vel(:,:,:)

        write(*,*) 'Calculates unnormalized velocity autocorrelation '
        write(*,*) 'function for a single species of atoms in a '
        write(*,*) 'simulation cell'
        write(*,*) 'Output: vaf.dat (time, vaf)'
        write(*,*) '[time in units of ps, vaf in (atomic units)^2]'
        write(*,*) ''
        write(*,*) 'Which file? '
        read(*,*) filename
        open(unit=10,file=filename,status='old')
        write(*,*) 'Start analysis from timestep #?'
        read(*,*) startt
        !write(*,*) 'End analysis at timestep #? (0 = end of file)'
        !read(*,*) nfr
        !if(nfr.eq.0) nfr=1e5
        nfr=1e5
        write(*,*) 'Simulation timestep [in au]? '
        read(*,*) ts
        !write(*,*) 'VAF cutoff time [in ps, 0 = half of simulation] ?'
        !read(*,*) cutt
        cutt=0
        ivel=0
        do ifr=1,nfr
         read(10,*,end=10) nion
         read(10,*,end=10) filename
         if(ifr.eq.1) then
          allocate(xyz(3,nion))
          allocate(xyz_old(3,nion))
          allocate(vel(3,nion,nfr-startt+1))
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
          vel(:,:,ivel)=(xyz(:,:)-xyz_old(:,:))/
     &                  (ts*.529177249) ! convert to atomic units
         endif
         xyz_old=xyz
        enddo
10      continue
        close(10)
        deallocate(xyz)
        deallocate(xyz_old)
        nfr=ivel
        if(cutt.eq.0) cutt=0.5*2.4189e-5*ts*nfr

        open(unit=15,file='vaf.dat',status='unknown')

        !write(15,'(f9.4,3x,f11.6)') 0.,1.
        vac0=0.0
        do ifr=1,nfr
         do iat=1,nion
          vac0=vac0+vel(1,iat,ifr)**2
     &             +vel(2,iat,ifr)**2
     &             +vel(3,iat,ifr)**2
         enddo
        ! if(ifr.eq.1) write(*,*) 'first',vac0/(1.0*nion)
        enddo
        vac0=vac0/(1.0*nion*nfr)
        write(15,*) 0.,vac0
        do k=1,nfr-1
         if(k.gt.500) then
          dt=500
         else
          dt=k
         endif
         nit=nfr-k
         vac=0.0
         do m=k+1,nfr
          mm=m-k
          do iat=1,nion
           vac=vel(1,iat,m)*vel(1,iat,mm)
     &        +vel(2,iat,m)*vel(2,iat,mm)
     &        +vel(3,iat,m)*vel(3,iat,mm)+vac
          end do
         end do
         currt=ts*k*2.4189e-5
        ! write(15,'(f9.4,3x,f11.6)') currt,vac/(1.0*nion*nit*vac0)
        ! write(15,*) currt,vac/(1.0*nion*nit*vac0)
         write(15,*) currt,vac/(1.0*nion*nit)
         if(currt.gt.cutt) goto 20
        end do
20      continue        
        deallocate(vel)
        close(15)
        end
