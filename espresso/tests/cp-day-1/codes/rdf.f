        program pair_corr

        real,allocatable :: xyz(:,:),gofr(:)
        character*150 filename
        character*2 atom
        real gridr,a,gfact,p,ri,ro,dist,meangrid
        integer nion,grid,bin,totfr,startt,fr,nfr

        write(*,*) 'Calculates radial pair distribution function '
        write(*,*) 'for a single species of atoms in a cubic cell'
        write(*,*) 'Output: rdf.dat (distance, g(r))'
        write(*,*) '[distance in same units as file]'
        write(*,*) ''
        write(*,*) 'Which file? '
        read(*,*) filename
        open(10,file=filename,status='old')
        p=3.14159265358
        write(*,*) 'Start analysis from timestep #? '
        read(*,*) startt
        !write(*,*) 'End analysis at timestep #? (0 = end of file)'
        !read(*,*) nfr
        !if(nfr.eq.0) nfr=1e5 ! max number of frames
        nfr=1e5
        write(*,*) 'Cubic lattice parameter (in atomic units)? '
        read(*,*) a
        a=a*0.529177249
        write(*,*) 'RDF grid spacing Delta_r (in atomic units)? ' 
        read(*,*) gridr
        gridr=gridr*0.529177249
        grid=int(a/gridr)+1
        allocate(gofr(grid))
        totfr=0
        gofr=0.0
        do ifr=1,nfr
         read(10,*,end=10) nion
         read(10,*,end=10) filename
         if(ifr.eq.1) allocate(xyz(3,nion))
         do iat=1,nion
          read(10,*,end=10) atom,(xyz(ixyz,iat),ixyz=1,3)
         enddo
         if(ifr.ge.startt) then
          totfr=totfr+1
          do ni=1,nion-1
           do ixyz=1,3
            rel=xyz(ixyz,iat)/a
            if(rel.lt.0.0) then
             xyz(ixyz,iat)=a*(rel-int(rel-1))
            else
             xyz(ixyz,iat)=a*(rel-int(rel))
            endif
           enddo
           do ii=-1,1
            do jj=-1,1
             do kk=-1,1
              do nii=ni+1,nion
               dist=sqrt((xyz(1,ni)-xyz(1,nii)+1.0*ii*a)**2+
     &                   (xyz(2,ni)-xyz(2,nii)+1.0*jj*a)**2+
     &                   (xyz(3,ni)-xyz(3,nii)+1.0*kk*a)**2)
               bin=int(dist/gridr)+1
               if(bin.le.grid) then
                gofr(bin)=gofr(bin)+1.0
               endif
              enddo
             enddo
            enddo
           enddo
          enddo
         endif
        enddo
10      continue
        close(10)
        open(15,file='rdf.dat',status='unknown')
        gfact=0.5*nion*(nion-1)/(a**3)*totfr*(4.0*p/3.0)
        do igrid=1,grid
         ro=1.0*igrid*gridr
         ri=1.0*(igrid-1)*gridr
         gofr(igrid)=gofr(igrid)/(ro**3-ri**3)
         gofr(igrid)=gofr(igrid)/gfact
         meangrid=gridr*(igrid-0.5)
         write(15,'(2F11.6)') meangrid,gofr(igrid)
        enddo
        deallocate(gofr)
        deallocate(xyz)
        close(15)
       end
