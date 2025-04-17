module setup
use mpi
use parameters
use flowvariables
use solvers


implicit none

contains
!***********************************************************************
!***********************************************************************
subroutine allocate_variabs()

!....wavenumbers
allocate(kx(0:Nx-1),kzloc(0:Nzloc-1))

!....dealiasing array
allocate(dealiasing(0:Ny,0:Nx-1,0:Nzloc-1))

!....3d wavenumbers
allocate(ikkx(0:Ny,0:Nx-1,0:Nzloc-1),ikkz(0:Ny,0:Nx-1,0:Nzloc-1))

!.... k^2
allocate(kkquad(0:Nx-1,0:Nzloc-1))
allocate(kkquadno0(0:Nx-1,0:Nzloc-1))

		
!....space coordinates		
allocate(x(0:Nx),zloc(0:Nzloc),y(0:Ny)) 
allocate(yp(0:Ny)) 

!....chebyshev diff matrix and diagonalization
allocate(D2(0:Ny,0:Ny),D(0:Ny,0:Ny),D2tilde(1:Ny-1,1:Ny-1))	
allocate(Qmat(1:Ny-1,1:Ny-1),QmatInv(1:Ny-1,1:Ny-1))
allocate(Lambda(1:Ny-1))	
allocate(Lambda2(1:Ny-1))

!....Fields in Physical Space
allocate(u0(0:Ny,0:Nx,0:Nzloc),v0(0:Ny,0:Nx,0:Nzloc),w0(0:Ny,0:Nx,0:Nzloc))
allocate(uMT(0:Ny,0:Nx,0:Nzloc),vMT(0:Ny,0:Nx,0:Nzloc),wMT(0:Ny,0:Nx,0:Nzloc))

!....Turbulent Kinetic Energy
allocate(Kturb(0:Ny))
allocate(ProdTurb(0:Ny))
allocate(Dissip(0:Ny,0:Nx-1,0:Nzloc-1))
allocate(DissipMT(0:Ny,0:Nx-1,0:Nzloc-1))
allocate(TurbTr(0:Ny,0:Nx-1,0:Nzloc-1))
allocate(TurbTrMT(0:Ny,0:Nx-1,0:Nzloc-1))
allocate(PressTr(0:Ny,0:Nx-1,0:Nzloc-1))
allocate(PressTrMT(0:Ny,0:Nx-1,0:Nzloc-1))
allocate(epsilonM(0:Ny))
allocate(DiffTurb(0:Ny))
allocate(DiffMolec(0:Ny))
allocate(DiffPress(0:Ny))


!....Fields in Fourier Space
allocate(uhat(0:Ny,0:Nx-1,0:Nzloc-1),vhat(0:Ny,0:Nx-1,0:Nzloc-1),what(0:Ny,0:Nx-1,0:Nzloc-1))
!allocate(uhatp1(0:Ny,0:Nx-1,0:Nzloc-1),vhatp1(0:Ny,0:Nx-1,0:Nzloc-1),whatp1(0:Ny,0:Nx-1,0:Nzloc-1))
allocate(uhatMT(0:Ny,0:Nx-1,0:Nzloc-1),vhatMT(0:Ny,0:Nx-1,0:Nzloc-1),whatMT(0:Ny,0:Nx-1,0:Nzloc-1))
allocate(ConvStr(0:Ny,0:Nx-1,0:Nzloc-1,1:6))
!allocate(ConvStrp1(0:Ny,0:Nx-1,0:Nzloc-1,1:6))
allocate(ReStress(0:Ny,0:Nx-1,0:Nzloc-1,1:6))
allocate(ReStressPh(0:Ny,0:Nx-1,0:Nzloc-1,6))

! !....Variables for Passive Scalars

! allocate(Theta(0:Ny,0:Nx,0:Nzloc))
! allocate(ThetaMT(0:Ny,0:Nx,0:Nzloc))
! allocate(ThetaSpz(0:Ny,0:Nx,0:Nzloc-1))
! allocate(ThetaSpz2(0:Ny,0:Nzloc-1))
! allocate(ThetaSpzMT(0:Ny,0:Nzloc-1))
! allocate(Thetahat(0:Ny,0:Nx-1,0:Nzloc-1))
! allocate(ThetahatMT(0:Ny,0:Nx-1,0:Nzloc-1))
! allocate(ThetahatSpz(0:Ny,0:Nx-1,0:Nzloc-1))


! allocate(ConvStrSc(0:Ny,0:Nx-1,0:Nzloc-1,1:4))
! allocate(TheStress(0:Ny,0:Nx-1,0:Nzloc-1,1:4))
! allocate(TheStressPh(0:Ny,0:Nx-1,0:Nzloc-1,1:4))
! allocate(dThetady(0:Ny))



!....Indices of Savings
allocate(IndSave(1:Nsavings))
allocate(Ubulk(1:Nsavings))
allocate(uumax(1:Nsavings))

allocate(nneven(0:Ny/2),nnodd(0:Ny/2-1))


!....Grid Spacing in the wall-normal direction
allocate(deltay(0:Ny-1))


!..Velocity incrememts 
allocate(deltux1(0:Ny,0:Nx-1,0:Nzloc-1))
allocate(deltuy1(0:Ny,0:Nx-1,0:Nzloc-1))
allocate(deltuz1(0:Ny,0:Nx-1,0:Nzloc-1))
allocate(deltux2(0:Ny,0:Nx-1,0:Nzloc-1))
allocate(deltuy2(0:Ny,0:Nx-1,0:Nzloc-1))
allocate(deltuz2(0:Ny,0:Nx-1,0:Nzloc-1))

!....Spacing of the velocity increments in the wall-normal directions
allocate(SpacY1Up(0:Ny))
allocate(SpacY1Do(0:Ny))
allocate(SpacY2Up(0:Ny))
allocate(SpacY2Do(0:Ny))
allocate(SpacY1rUp(0:Ny))
allocate(SpacY1rDo(0:Ny))
allocate(SpacY2rUp(0:Ny))
allocate(SpacY2rDo(0:Ny))

allocate(InDom1Up(0:Ny))
allocate(InDom1Do(0:Ny))
allocate(InDom2Up(0:Ny))
allocate(InDom2Do(0:Ny))

allocate(InDom1(0:Ny))
allocate(InDom2(0:Ny))

!!....Allocate Particles Variables
!allocate(Part(1:PartInfo,1:Nplocmax))
!allocate(PartBuffer(1:PartInfo,1:Nplocmax))
!allocate(PartIndices(1:3,1:Nplocmax))
!allocate(Weights(1:8,1:Nplocmax))



!#ifdef SCAL
!deallocate(ConvStrScp1)
!#endif
!deallocate(ConvStrp1)
!allocate(ReStressPh(0:Ny,0:Nx-1,0:Nzloc-1,6))


end subroutine
!***********************************************************************
!***********************************************************************
subroutine deallocate_variabs()

!....wavenumbers
deallocate(kx,kzloc)

!....dealiasing array
deallocate(dealiasing)

!....3d wavenumbers
deallocate(ikkx,ikkz)

!....k^2
deallocate(kkquad)
deallocate(kkquadno0)


!....space coordinates		
deallocate(x,zloc,y) 
deallocate(yp) 

!....chebyshev diff matrix and diagonalization
deallocate(D2,D,D2tilde)	
deallocate(Qmat,QmatInv)
deallocate(Lambda)	
deallocate(Lambda2)	

!....Fields in Physical Space
deallocate(u0,v0,w0)
deallocate(uMT,vMT,wMT)

!....Turbulent Kinetic Energy
deallocate(Kturb)
deallocate(ProdTurb)
deallocate(Dissip)
deallocate(DissipMT)
deallocate(TurbTr)
deallocate(TurbTrMT)
deallocate(PressTr)
deallocate(PressTrMT)
deallocate(epsilonM)
deallocate(DiffTurb)
deallocate(DiffMolec)
deallocate(DiffPress)


!....Fields in Fourier Space
deallocate(uhat,vhat,what)
deallocate(uhatMT,vhatMT,whatMT)
deallocate(ReStress)
deallocate(ReStressPh)
deallocate(ConvStr)


! !....Variables for Passive Scalars
! deallocate(Theta)
! deallocate(ThetaMT)
! deallocate(ThetaSpz)
! deallocate(ThetaSpz2)
! deallocate(ThetaSpzMT)
! deallocate(Thetahat)
! deallocate(ThetahatMT)
! deallocate(ThetahatSpz)
! deallocate(TheStress)
! deallocate(TheStressPh)
! deallocate(ConvStrSc)
! deallocate(dThetady)


!....Indices of Savings
deallocate(IndSave)

deallocate(Ubulk)
deallocate(uumax)

deallocate(nneven,nnodd)

!!....Deallocate Particles Variables
!deallocate(Part)
!deallocate(PartBuffer)
!deallocate(PartIndices)
!deallocate(Weights)

!....Grid Spacing in the wall-normal direction
deallocate(deltay)



!..Velocity incrememts 
deallocate(deltux1)
deallocate(deltuy1)
deallocate(deltuz1)
deallocate(deltux2)
deallocate(deltuy2)
deallocate(deltuz2)


!....Spacing of the velocity increments in the wall-normal directions
deallocate(SpacY1Up)
deallocate(SpacY1Do)
deallocate(SpacY2Up)
deallocate(SpacY2Do)
deallocate(SpacY1rUp)
deallocate(SpacY1rDo)
deallocate(SpacY2rUp)
deallocate(SpacY2rDo)


deallocate(InDom1Up)
deallocate(InDom1Do)
deallocate(InDom2Up)
deallocate(InDom2Do)

deallocate(InDom1)
deallocate(InDom2)

end subroutine
!***********************************************************************
!***********************************************************************
subroutine InitializeFlow()
  implicit none
  integer::ii,jj,ll
  real(kind=prec)::y_target1,y_target2

  
  jj=1
  do ii=Nmin,Nmax,Ipasso
  IndSave(jj)=ii
  jj=jj+1
  end do
  
  
  !....Some initial definitions
  !...Domain Size
  Lx=Lx_over_hpi*pi
  Lz=Lz_over_hpi*pi
  Lzloc=Lz/noprocs
  Ly=2.0d0
  
  !...Spacing in x and z
  deltax=Lx/dfloat(Nx)
  deltaz=Lz/dfloat(Nz)
  
  
  
  !....Processors id
  if(nid.eq.0)then
  nidp1=nid+1
  nidm1=noprocs-1
  elseif(nid.eq.noprocs-1)then
  nidp1=0
  nidm1=noprocs-2
  else
  nidp1=nid+1
  nidm1=nid-1
  endif
  
  
  !....wavenumber calculation
  allocate(kz(0:Nz-1))
  do ii=0,Nx/2
  kx(ii)=dfloat(ii)
  end do
  do ii=1,Nx/2-1
  kx(Nx-ii)=-dfloat(ii)
  end do  
  
  do ii=0,Nz/2
  kz(ii)=dfloat(ii)
  end do
  do ii=1,Nz/2-1
  kz(Nz-ii)=-dfloat(ii)
  end do  
  
  do ii=0,Nzloc-1
  kzloc(ii)=kz(nid*Nzloc+ii)
  end do
  
  deallocate(kz)
  
  !....dealiasing array calculation
  do ll=0,Nzloc-1
  do jj=0,Nx-1
  do ii=0,Ny
  dealiasing(ii,jj,ll)=.false.
  end do
  end do
  end do
  
  do ll=0,Nzloc-1
  do jj=0,Nx-1
  do ii=0,Ny
  if((dabs(kx(jj)/Nx)>(1.0d0/3.0d0)).OR.(dabs(kzloc(ll)/Nz)>(1.0d0/3.0d0))) then
  dealiasing(ii,jj,ll)=.true.
  end if
  end do
  end do
  end do
  
  
  call MPI_BARRIER(MPI_Comm_World,ierror) 
  
  
  kx(Nx/2)=0.0d0
  	
  
  do ii=0,Nzloc-1
  if((nid*Nzloc+ii).eq.(Nz/2)) then
  kzloc(ii)=0.0d0
  end if
  end do
  
  
  !....Adapt the wavenumbers to the domain sizes
  kx=kx*2.0d0*pi/Lx
  kzloc=kzloc*2.0d0*pi/Lz
  
  
  
  do ll=0,Nzloc-1
  do ii=0,Nx-1
  do jj=0,Ny
  ikkx(jj,ii,ll)=imu*kx(ii)
  ikkz(jj,ii,ll)=imu*kzloc(ll)
  end do
  kkquad(ii,ll)=(kx(ii))**2+(kzloc(ll))**2
  end do
  end do
  
  
  kkquadno0=kkquad
  do ll=0,Nzloc-1
  do jj=0,Nx-1
  if(kkquadno0(jj,ll).eq.0.0d0)then
  kkquadno0(jj,ll)=1.0d0
  end if
  end do
  end do
  
  
  allocate(z(0:Nz)) 
  do ii=0,Nx
  x(ii)=Lx*(ii*1.0d0)/(Nx)
  end do	
  do ii=0,Nz
  z(ii)=Lz*(ii*1.0d0)/(Nz)
  end do		
  do ii=0,(Ny/2-1)
  y(ii)=dcos(pi*ii/Ny)
  y(Ny-ii)=-y(ii)
  end do
  y(Ny/2)=0.0d0
  
  
  do ii=0,Ny
  yp(ii)=(1.0d0-y(ii))*Re
  end do
  
  do ii=0,Nzloc
  zloc(ii)=z(nid*Nzloc+ii)
  end do
  
  deallocate(z)
  
  do ii=0,Ny/2
  nneven(ii)=ii*2
  end do
  do ii=0,Ny/2-1
  nnodd(ii)=ii*2+1
  end do
  
  
  !....Chebyshev Differentiation Matrices calculation and diagonlization
  
  call ChebyshevDiffMatrixV4Peyret()
  
  ! !....Chebyshev-tau Differentiation Matrices calculation
  ! call ChebyshevTauDiffMatrix()
  
  ! do ii=0,Ny
  ! Tnp1(ii)=(1.0d0)**(ii)
  ! Tnm1(ii)=(-1.0d0)**(ii)
  ! end do
  
  
  ! temp1=0.0d0
  
  ! do ii=0,Ny-1
  ! temp1(ii)=-DCOS(ii*pi/Ny+pi/2.0d0/Ny)*DCOS(pi/2/Ny)
  ! end do
  
  ! do ll=0,Nzloc-1
  ! do jj=0,Nx-1
  ! do ii=1,Ny-1
  ! Volumes(ii,jj,ll)=(temp1(ii)-temp1(ii-1))*deltax*deltaz
  ! end do
  ! end do
  ! end do
  
  ! Volumes(0,0:Nx-1,0:Nzloc-1)=deltax*deltaz*(1.0d0+temp1(0))
  ! Volumes(Ny,0:Nx-1,0:Nzloc-1)=Volumes(0,0:Nx-1,0:Nzloc-1)
  
  
  
  ! deallocate(temp1)
  
  !Grid Spacing in the wall-normal direction
  do ii=0,Ny-1
      deltay(ii)=DABS(y(ii+1)-y(ii))
  end do
  

  !...Velocity increments 
  SpacX1 = FLOOR(sep1/deltax)
  SpacX2 = FLOOR(sep2/deltax)
  SpacZ1 = FLOOR(sep1/deltaz)
  SpacZ2 = FLOOR(sep2/deltaz)
  
  !...Remainder
  SpacX1r = sep1/deltax - FLOOR(sep1/deltax)
  SpacX2r = sep2/deltax - FLOOR(sep2/deltax)
  SpacZ1r = sep1/deltaz - FLOOR(sep1/deltaz)
  SpacZ2r = sep2/deltaz - FLOOR(sep2/deltaz)
  
  
  
  ! Do not consider the boundaries (walls)
  ! Initialize everything
  InDom1Up= .false.
  InDom1Do= .false.
  InDom2Up= .false.
  InDom2Do= .false.

  InDom1 = .false.
  InDom2 = .false.

  SpacY1Do = -1
  SpacY2Do = -1
  SpacY1rDo = 0.0d0
  SpacY2rDo = 0.0d0

  SpacY1Up = -1
  SpacY2Up = -1
  SpacY1rUp = 0.0d0
  SpacY2rUp = 0.0d0

  do ii = 1, Ny-1  ! Skip walls
      ! Downward target positions (y decreases downward, so we subtract)
      y_target1 = y(ii) - sep1
      y_target2 = y(ii) - sep2
      
      ! Search for y_target1 in grid
      do jj = ii, Ny-1  ! move downward (increasing jj)
        if ((y(jj) >= y_target1) .and. (y(jj+1) <= y_target1)) then
          SpacY1Do(ii) = jj - ii
          SpacY1rDo(ii) = DABS((y_target1 - y(jj)) / (y(jj+1) - y(jj)))
          InDom1Do(ii) = .true.
          exit
        endif
      end do
      
      ! Same for sep2
      do jj = ii, Ny-1
        if ((y(jj) >= y_target2) .and. (y(jj+1) <= y_target2)) then
          SpacY2Do(ii) = jj - ii
          SpacY2rDo(ii) = DABS((y_target2 - y(jj)) / (y(jj+1) - y(jj)))
          InDom2Do(ii) = .true.
          exit
        endif
      end do    
     
         
     ! Upward target positions (physically upward means +sep, index decreases)
     y_target1 = y(ii) + sep1
     y_target2 = y(ii) + sep2
  
     ! Search for y_target1 in grid (moving upward: decreasing jj)
     do jj = ii, 1, -1
        if (jj > 0) then
          if ((y(jj) <= y_target1) .and. (y_target1 <= y(jj-1))) then
            SpacY1Up(ii) = jj - ii
            SpacY1rUp(ii) = DABS((y_target1 - y(jj)) / (y(jj-1) - y(jj)))
            InDom1Up(ii) = .true.
            exit
          end if
        end if
      end do

    
     ! Same for y_target2
      do jj = ii, 1, -1
        if (jj > 0) then
          if ((y(jj) <= y_target2) .and. (y_target2 <= y(jj-1))) then
            SpacY2Up(ii) = jj - ii
            SpacY2rUp(ii) = DABS((y_target2 - y(jj)) / (y(jj-1) - y(jj)))
            InDom2Up(ii) = .true.
            exit
          end if
        end if
      end do
  end do

  do ii = 1, Ny-1
    InDom1(ii) = InDom1Up(ii) .and. InDom1Do(ii)
    InDom2(ii) = InDom2Up(ii) .and. InDom2Do(ii)
  end do






return
end subroutine
!***********************************************************************
subroutine ChebyshevDiffMatrixV4Peyret()
  integer::ii,jj
  real(kind=prec),allocatable,dimension(:)::temp
  
  allocate(temp(0:Ny))
  
  
  !....First order derivative
  D=0.0d0
  
  !....Off-Diagonal Entries
  do ii=1,Ny-1
  do jj=1,Ny-1
  if (ii.ne.jj)then
  D(ii,jj)=((-1.0d0)**(ii+jj))/(y(ii)-y(jj))
  end if 
  end do
  end do
  
  do jj=1,Ny-1
  D(0,jj)=2.0d0*((-1.0d0)**(jj))/(1-y(jj))
  end do
  
  D(0,Ny)=0.5d0*(-1.0d0)**(Ny)
  
  do ii=1,Ny-1
  D(ii,Ny)=0.5d0*((-1.0d0)**(Ny+ii))/(1+y(ii))
  end do
  
  do ii=1,Ny-1
  D(ii,0)=-0.5d0*((-1.0d0)**(ii))/(1-y(ii))
  end do
  
  D(Ny,0)=-D(0,Ny)
  
  do jj=1,Ny-1
  D(Ny,jj)=-2.0d0*((-1.0d0)**(Ny+jj))/(1+y(jj))
  end do
  
  !....Diagonal Entries
  do jj=0,Ny
  temp=0.0d0
  do ii=0,Ny
  if(ii.ne.jj)then
  temp(ii)=D(jj,ii)
  end if
  end do
  
  D(jj,jj)=-sum(temp)
  end do
  
  !....Second order derivative
  D2=0.0d0
  
  D2=matmul(D,D)
  
  !....Removing Diagonal Entries
  do jj=0,Ny
  D2(jj,jj)=0.0d0
  end do
  
  !....Diagonal Entries
  do jj=0,Ny
  temp=0.0d0
  do ii=0,Ny
  if(ii.ne.jj)then
  temp(ii)=D2(jj,ii)
  end if
  end do
  
  D2(jj,jj)=-sum(temp)
  end do
  
  
  
  do ii=1,Ny-1
  do jj=1,Ny-1
  D2tilde(jj,ii)=D2(jj,ii)	
  end do
  end do
  
  
  !....Matrix Diagonalization
  !Remember: MatrixDiagonalization(Eigenv,Pmat,PmatInv,Matrice,NN)
  call MatrixDiagonalization(Lambda,Qmat,QmatInv,D2tilde,Ny-1)
  
  	
  deallocate(temp)
  
  return
end subroutine
!***********************************************************************
!***********************************************************************
end module setup
