        Program Mean_square_displacement

        real,allocatable :: xyz_orig(:,:),xyz(:,:)
        integer startt,fr,nfr
        character*150 filename
        character*2 atom
        real ts,msd,com_orig(3),com(3)

        write(*,*) 'Calculates mean square displacement for a single '
        write(*,*) 'species of atoms in a simulation cell '
        write(*,*) 'Output: msd.dat (time, distance^2)'
        write(*,*) '[time in units of iteration no, distance in Ang]'
        write(*,*) ''
        write(*,*) 'Which file? '
        read(*,*) filename
        open(unit=10,file=filename,status='old')
        write(*,*) 'Start analysis from timestep #? '
        read(*,*) startt
        !write(*,*) 'End analysis at timestep #? (0 = end of file)'
        !read(*,*) nfr
        !if(nfr.eq.0) nfr=1e5 ! max number of frames
        nfr=1e5
        !write(*,*) 'Analyze every nth timestep [all=1]? '
        !read(*,*) fr
        fr=1
        write(*,*) 'Timestep increment (atomic units)? '
        read(*,*) ts
        ts=ts*2.4189e-5

        open(unit=15,file='msd.dat',status='unknown')
        write(15,'(f11.6,f11.5)') 0,0

        do ifr=1,startt 
         read(10,*) nion
         read(10,*) filename
         if(ifr.eq.1) then
          allocate(xyz_orig(3,nion))
          allocate(xyz(3,nion))
         endif
         do iat=1,nion
          read(10,*) atom,(xyz_orig(ixyz,iat),ixyz=1,3)
         enddo
        enddo
        com_orig(:)=0.0
        do iat=1,nion
         com_orig(:)=com_orig(:)+xyz_orig(:,iat)
        enddo
        com_orig(:)=com_orig(:)/(1.0*nion)
        do iat=1,nion
         xyz_orig(:,iat)=xyz_orig(:,iat)-com_orig(:)
        enddo
        nfr=nfr-startt+1
        do ifr=2,nfr
         msd=0.0
         read(10,*,end=10) nion
         read(10,*,end=10) filename
         com(:)=0.0
         do iat=1,nion
          read(10,*,end=10) atom,(xyz(ixyz,iat),ixyz=1,3)
          com(:)=com(:)+xyz(:,iat)
         enddo
         com(:)=com(:)/(1.0*nion)
         do iat=1,nion
          xyz(:,iat)=xyz(:,iat)-com(:)
          msd=((xyz(1,iat)-xyz_orig(1,iat))**2)+
     &        ((xyz(2,iat)-xyz_orig(2,iat))**2)+
     &        ((xyz(3,iat)-xyz_orig(3,iat))**2)+msd
         enddo
         msd=msd/(1.0*nion)
        ! write(15,'(f11.6,f11.5)') ts*(ifr-1),msd
         write(15,'(i9,f11.5)') ifr,msd
        enddo
10      continue
        deallocate(xyz_orig)
        deallocate(xyz)
        close(10)
        close(15)
        end
