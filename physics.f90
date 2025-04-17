module physics
use mpi
use parameters
use transforms

implicit none


contains
!***********************************************************************
!***********************************************************************
subroutine ConvectiveStress(ConvSt,u,v,w)
implicit none 
complex(kind=prec)::ConvSt(0:Ny,0:Nx-1,0:Nzloc-1,1:6)
real(kind=prec)::u(0:Ny,0:Nx-1,0:Nzloc-1)
real(kind=prec)::v(0:Ny,0:Nx-1,0:Nzloc-1)
real(kind=prec)::w(0:Ny,0:Nx-1,0:Nzloc-1)
integer::ii,jj,ll
real(kind=prec),allocatable,dimension(:,:,:):: tempr1



allocate(tempr1(0:Ny,0:Nx-1,0:Nzloc-1))


!....uu
tempr1=u*u
call transform_rows_batched(ConvSt(0:Ny,0:Nx-1,0:Nzloc-1,1),tempr1)!(u*u)_hat

!....vv
tempr1=v*v
call transform_rows_batched(ConvSt(0:Ny,0:Nx-1,0:Nzloc-1,2),tempr1)!(v*v)_hat

!....ww
tempr1=w*w
call transform_rows_batched(ConvSt(0:Ny,0:Nx-1,0:Nzloc-1,3),tempr1)!(w*w)_hat

!....uv
tempr1=u*v
call transform_rows_batched(ConvSt(0:Ny,0:Nx-1,0:Nzloc-1,4),tempr1)!(u*v)_hat

!....uw
tempr1=u*w
call transform_rows_batched(ConvSt(0:Ny,0:Nx-1,0:Nzloc-1,5),tempr1)!(u*w)_hat


!....vw
tempr1=v*w
call transform_rows_batched(ConvSt(0:Ny,0:Nx-1,0:Nzloc-1,6),tempr1)!(v*w)_hat


!....dealiasing the convective stress
where(dealiasing) 
ConvSt(0:Ny,0:Nx-1,0:Nzloc-1,1)=cpx0
ConvSt(0:Ny,0:Nx-1,0:Nzloc-1,2)=cpx0
ConvSt(0:Ny,0:Nx-1,0:Nzloc-1,3)=cpx0
ConvSt(0:Ny,0:Nx-1,0:Nzloc-1,4)=cpx0
ConvSt(0:Ny,0:Nx-1,0:Nzloc-1,5)=cpx0
ConvSt(0:Ny,0:Nx-1,0:Nzloc-1,6)=cpx0
end where




deallocate(tempr1)

return
end subroutine
!***********************************************************************
!***********************************************************************
subroutine EvalSquaredVelIncrementsNormal()
    ! Compute the velocity increments squared in the three directions
    ! these are normal to the chosen direction
    implicit none 
    integer(kind=prec):: ii,jj,ll,hhLef,hhRig
    real(kind=prec):: deltuxL,deltuyL,deltuzL
    real(kind=prec):: deltuxLfwd,deltuxLbwd

    !....Compute velocity increments 
    ! in streamwise (x)
    do ll=0,Nzloc-1
      do jj=0,Nx-1
        do ii=0,Ny
            
         !....First Separation
            !...Indices of the separated point 
            hhLef = jj + SpacX1
            hhRig = jj + SpacX1 + 1 

            


            ! deltuxLfwd = 

         !....Second Separation

            

        end do
      end do
    end do
            


    return
end subroutine
!***********************************************************************
subroutine DissipationTerm(Disspdmy,uhatdmy,vhatdmy,whatdmy)
implicit none
real(kind=prec)::Disspdmy(0:Ny,0:Nx-1,0:Nzloc-1) 
complex(kind=prec)::uhatdmy(0:Ny,0:Nx-1,0:Nzloc-1)
complex(kind=prec)::vhatdmy(0:Ny,0:Nx-1,0:Nzloc-1)
complex(kind=prec)::whatdmy(0:Ny,0:Nx-1,0:Nzloc-1)
complex(kind=prec),allocatable,dimension(:,:,:)::tempcp
real(kind=prec),allocatable,dimension(:,:,:)::tempr
integer::jj,ll


allocate(tempcp(0:Ny,0:Nx-1,0:Nzloc-1))
allocate(tempr(0:Ny,0:Nx-1,0:Nzloc-1))

!...Zeroing the Output
Disspdmy=0.0d0

!Nine Terms have to be evaluated (dui/dxj*dui/dxj)


!....DERIVATIVES for u
!...(du/dx*du/dx)
tempcp=ikkx*uhatdmy!(du/dx)_hat
call transformInv_rows_batched(tempr,tempcp)! du/dx
Disspdmy=tempr*tempr

!...(du/dy*du/dy)
do ll=0,Nzloc-1
do jj=0,Nx-1
tempcp(0:Ny,jj,ll)=matmul(D(0:Ny,0:Ny),uhatdmy(0:Ny,jj,ll))!(du/dy)_hat
end do
end do
 
call transformInv_rows_batched(tempr,tempcp)! du/dy
Disspdmy=Disspdmy+tempr*tempr


!...(du/dz*du/dz)
tempcp=ikkz*uhatdmy!(du/dz)_hat 
call transformInv_rows_batched(tempr,tempcp)! du/dz
Disspdmy=Disspdmy+tempr*tempr



!....DERIVATIVES for v
!...(dv/dx*dv/dx)
tempcp=ikkx*vhatdmy!(dv/dx)_hat
call transformInv_rows_batched(tempr,tempcp)! dv/dx
Disspdmy=Disspdmy+tempr*tempr

!...(dv/dy*dv/dy)
do ll=0,Nzloc-1
do jj=0,Nx-1
tempcp(0:Ny,jj,ll)=matmul(D(0:Ny,0:Ny),vhatdmy(0:Ny,jj,ll))!(dv/dy)_hat
end do
end do
 
call transformInv_rows_batched(tempr,tempcp)! dv/dy
Disspdmy=Disspdmy+tempr*tempr


!...(dv/dz*dv/dz)
tempcp=ikkz*vhatdmy!(dv/dz)_hat 
call transformInv_rows_batched(tempr,tempcp)! dv/dz
Disspdmy=Disspdmy+tempr*tempr


!....DERIVATIVES for w
!...(dv/dx*dv/dx)
tempcp=ikkx*whatdmy!(dw/dx)_hat
call transformInv_rows_batched(tempr,tempcp)! dw/dx
Disspdmy=Disspdmy+tempr*tempr

!...(dw/dy*dw/dy)
do ll=0,Nzloc-1
do jj=0,Nx-1
tempcp(0:Ny,jj,ll)=matmul(D(0:Ny,0:Ny),whatdmy(0:Ny,jj,ll))!(dw/dy)_hat
end do
end do
 
call transformInv_rows_batched(tempr,tempcp)! dw/dy
Disspdmy=Disspdmy+tempr*tempr


!...(dw/dz*dw/dz)
tempcp=ikkz*whatdmy!(dw/dz)_hat 
call transformInv_rows_batched(tempr,tempcp)! dw/dz
Disspdmy=Disspdmy+tempr*tempr


!....NOW:  Disspdmy=dui/dxj*dui/dxj

deallocate(tempr)
deallocate(tempcp)


return
end subroutine
!***********************************************************************
!***********************************************************************
subroutine TurbulentTransportTerm(TurbTrdmy,udmy,vdmy,wdmy)
implicit none 
real(kind=prec)::TurbTrdmy(0:Ny,0:Nx-1,0:Nzloc-1) 
real(kind=prec)::udmy(0:Ny,0:Nx-1,0:Nzloc-1) 
real(kind=prec)::vdmy(0:Ny,0:Nx-1,0:Nzloc-1)
real(kind=prec)::wdmy(0:Ny,0:Nx-1,0:Nzloc-1)
real(kind=prec),allocatable,dimension(:,:,:)::tempr

allocate(tempr(0:Ny,0:Nx-1,0:Nzloc-1))

TurbTrdmy=0.0d0
tempr=0.0d0


tempr=udmy*udmy+vdmy*vdmy+wdmy*wdmy

TurbTrdmy=0.5d0*(tempr*vdmy)


deallocate(tempr)

return
end subroutine
!***********************************************************************
!***********************************************************************
subroutine MeanXZ(averval,Sfield)
implicit none
real(kind=prec)::averval(0:Ny) 
real(kind=prec)::Sfield(0:Ny,0:Nx-1,0:Nzloc-1) 
complex(kind=prec),allocatable,dimension(:,:,:)::tempcp1,tempcp2
real(kind=prec),allocatable,dimension(:,:,:)::tempr

allocate(tempcp1(0:Ny,0:Nx-1,0:Nzloc-1))
allocate(tempcp2(0:Ny,0:Nx-1,0:Nzloc-1))
allocate(tempr(0:Ny,0:Nx-1,0:Nzloc-1))

averval=0.0d0
tempr=0.0d0
tempcp2=cpx0


call transform_rows_batched(tempcp1,Sfield)

if(nid.eq.0)then
tempcp2(0:Ny,0,0)=tempcp1(0:Ny,0,0)
end if

call transformInv_rows_batched(tempr,tempcp2)


averval(0:Ny)=tempr(0:Ny,0,0)




deallocate(tempcp1)
deallocate(tempcp2)
deallocate(tempr)

return
end subroutine
!***********************************************************************
!***********************************************************************
subroutine EvalMeanTurbTranspor(avervaldy,Sfield)
implicit none
real(kind=prec)::avervaldy(0:Ny) 
real(kind=prec)::Sfield(0:Ny,0:Nx-1,0:Nzloc-1) 
complex(kind=prec),allocatable,dimension(:,:,:)::tempcp1,tempcp2
real(kind=prec),allocatable,dimension(:,:,:)::tempr


allocate(tempcp1(0:Ny,0:Nx-1,0:Nzloc-1))
allocate(tempcp2(0:Ny,0:Nx-1,0:Nzloc-1))
allocate(tempr(0:Ny,0:Nx-1,0:Nzloc-1))

avervaldy=0.0d0
tempr=0.0d0
tempcp2=cpx0


call transform_rows_batched(tempcp1,Sfield)

if(nid.eq.0)then
tempcp2(0:Ny,0,0)=tempcp1(0:Ny,0,0)
end if

call transformInv_rows_batched(tempr,tempcp2)


avervaldy(0:Ny)=matmul(D(0:Ny,0:Ny),tempr(0:Ny,0,0))




deallocate(tempcp1)
deallocate(tempcp2)
deallocate(tempr)

return
end subroutine
!***********************************************************************
!***********************************************************************
subroutine EvalUbulk(Ubulkdmy,uhatdmy)
implicit none
real(kind=prec)::Ubulkdmy
complex(kind=prec)::uhatdmy(0:Ny) 
complex(kind=prec),allocatable, dimension(:,:,:)::uhatm
real(kind=prec),allocatable, dimension(:,:,:)::urm
real(kind=prec),allocatable,dimension(:)::temp1
real(kind=prec),allocatable, dimension(:,:,:)::temp2
integer::ii


allocate(uhatm(0:Ny,0:Nx-1,0:Nzloc-1))
allocate(urm(0:Ny,0:Nx-1,0:Nzloc-1))
allocate(temp2(0:Ny,0:Nx-1,0:Nzloc-1))
allocate(temp1(0:Ny))

uhatm=cpx0
urm=0.0d0
temp2=0.0d0

uhatm(0:Ny,0,0)=uhatdmy(0:Ny)

call transformInv_rows_batched(urm,uhatm)

!...Chebyshev Transform in y
call DCT_1D_FWD(temp2(0:Ny,0:Nx-1,0:Nzloc-1),urm(0:Ny,0:Nx-1,0:Nzloc-1))


! call ChebTransform1Dy(temp1,urm(0:Ny,0,0))

temp1(0:Ny)=temp2(0:Ny,0,0)


temp1(nnodd)=0.0d0

temp1(0)=-temp1(0)
do ii=2,Ny
temp1(ii)=temp1(ii)/((dfloat(ii)**2)-1.0d0)
end do

temp1(nnodd)=0.0d0

Ubulkdmy=-1.0d0*sum(temp1(0:Ny))




deallocate(uhatm)
deallocate(urm)
deallocate(temp1)
deallocate(temp2)
return
end subroutine
!***********************************************************************
!***********************************************************************
subroutine EvaluateMaxuu(maxdmy,Cfield)
real(kind=prec)::maxdmy
complex(kind=prec)::Cfield(0:Ny,0:Nx-1,0:Nzloc-1)
real(kind=prec),allocatable,dimension(:,:,:)::Rfield
real(kind=prec)::maxlocale,maxglob

allocate(Rfield(0:Ny,0:Nx-1,0:Nzloc-1))

Rfield=0.0d0

call transformInv_rows_batched(Rfield,Cfield)


!...Ogni processore calcola il massimo LOCALE
maxlocale=maxval(Rfield)

call MPI_Allreduce(maxlocale,maxglob,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_Comm_World,ierror)

maxdmy=maxglob




deallocate(Rfield)
return
end subroutine
!***********************************************************************
!***********************************************************************
!***********************************************************************
!***********************************************************************
subroutine MeanX(averval,Sfield)
implicit none
real(kind=prec)::averval(0:Ny,0:Nzloc-1)
real(kind=prec)::Sfield(0:Ny,0:Nx-1,0:Nzloc-1)
real(kind=prec),allocatable, dimension(:)::uaus
integer::ii,jj

allocate(uaus(0:Nx-1))

averval=0.0d0



do jj=0,Nzloc-1
do ii=0,Ny

uaus(0:Nx-1)=Sfield(ii,0:Nx-1,jj)
averval(ii,jj)=(sum(uaus(0:Nx-1)))/Nx

end do 
end do





deallocate(uaus)

return
end subroutine
!***********************************************************************
!***********************************************************************

end module physics

