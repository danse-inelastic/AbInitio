        Program Mean_square_displacement

        real,allocatable :: xyz(:,:,:)
        integer startt,fr,nfr,dt
        character*150 filename
        character*2 atom
        real ts,msd,com(3),com_orig(3)

        write(*,*) 'Calculates mean square displacement for a single '
        write(*,*) 'species of atoms in a simulation cell '
        write(*,*) 'Averages over all time windows in the simulation '
        write(*,*) 'for maximum statistics (use only in equilibrium!)'
        write(*,*) 'Output: msd.dat (time, distance^2)'
        write(*,*) '[time in units of ps, distance in Ang]'
        write(*,*) ''
        write(*,*) 'Which file? '
        read(*,*) filename
        open(unit=10,file=filename,status='old')
        write(*,*) 'Start analysis from timestep #? '
        read(*,*) startt
        write(*,*) 'End analysis at timestep #? (0 = end of file)'
        read(*,*) nfr
        if(nfr.eq.0) nfr=1e5 ! max number of frames
        !write(*,*) 'Analyze every nth timestep [all=1]? '
        !read(*,*) fr
        fr=1
        write(*,*) 'Simulation timestep (atomic units)? '
        read(*,*) ts
        ts=ts*2.4189e-5

        open(unit=15,file='msd.dat',status='unknown')
        write(15,'(f11.6,f11.5)') 0,0

        do ifr=1,startt
         read(10,*) nion
         read(10,*) filename
         if(ifr.eq.1) allocate(xyz(3,nion,nfr-startt+1))
         do iat=1,nion
          read(10,*) atom,(xyz(ixyz,iat,1),ixyz=1,3)
         enddo
        enddo
        com_orig(:)=0.0
        do iat=1,nion
         com_orig(:)=com_orig(:)+xyz(:,iat,1)
        enddo
        com_orig(:)=com_orig(:)/(1.0*nion)
        do iat=1,nion
         xyz(:,iat,1)=xyz(:,iat,1)-com_orig(:)
        enddo
        nxyz=1
        do ifr=startt+1,nfr
         com(:)=0.0
         read(10,*,end=10) nion
         read(10,*,end=10) filename
         do iat=1,nion
          read(10,*,end=10) atom,(xyz(ixyz,iat,nxyz+1),ixyz=1,3)
          com(:)=com(:)+xyz(:,iat,nxyz+1)
         enddo
         com(:)=com(:)/(1.0*nion)
         do iat=1,nion
          xyz(:,iat,nxyz+1)=xyz(:,iat,nxyz+1)-com(:)
         enddo
         nxyz=nxyz+1
        enddo
10      continue
        nfr=nxyz
        close(10)

        do k=1,nfr-1
         nit=nfr-k
         msd=0.0
         do m=k+1,nfr
          mm=m-k
          do iat=1,nion
           msd=((xyz(1,iat,m)-xyz(1,iat,mm))**2)+
     &         ((xyz(2,iat,m)-xyz(2,iat,mm))**2)+
     &         ((xyz(3,iat,m)-xyz(3,iat,mm))**2)+msd
          enddo
         enddo
         write(15,'(f11.6,f11.5)') ts*k,msd/(1.0*nion*nit)
        enddo
        deallocate(xyz)
        close(15)
        end

