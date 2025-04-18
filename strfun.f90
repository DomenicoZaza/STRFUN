program strfun
use, intrinsic :: iso_c_binding
use mpi
use parameters
use physics
use flowvariables
use in_out
use setup
use transforms


implicit none

integer:: ii,jj,ll

!.....Variables for parallelization:
integer status(MPI_Status_size)

!.....Parallel initialization:
call MPI_Init(ierror)
call MPI_Comm_size(MPI_COMM_WORLD,noprocs,ierror)!noprocs=number of processors
call MPI_Comm_rank(MPI_COMM_WORLD,nid,ierror)    !nid=rank of each processor
call fftw_mpi_init()

!.....End of initialization...
print *,'nid',nid
call MPI_BARRIER(MPI_Comm_World,ierror)    

!....Allocate namefiles
allocate(nomefile(1:3))
allocate(nomefile1(1:3))
allocate(resfile(1:3))
! allocate(gacc(1:2))


!.....Parameter are read by processor 0 from file in_dns_tc.txt
!and spread to all processors
call read_inputs()

call MPI_BARRIER(MPI_Comm_World,ierror) 

!...Allocate Variables
call allocate_variabs()

call MPI_BARRIER(MPI_Comm_World,ierror) 

!...Initialize calculation
call InitializeFlow()

! !...Initialize calculation
! call initialize()

call MPI_BARRIER(MPI_Comm_World,ierror) 

call transforms_init()

call MPI_BARRIER(MPI_Comm_World,ierror) 

!....Zeroing the mean fields
uhatMT=cpx0
vhatMT=cpx0
whatMT=cpx0
! ThetahatMT=cpx0
! ThetaSpzMT=0.0d0

!Zeroing variables for statistical stationary
Ubulk=0.0d0
uumax=0.0d0



if(nid.eq.0)then
WRITE(*,*)'Evaluating Mean Profiles.'
WRITE(*,*)'Processing step: ',IndSave(1)
end if

!....Generate names of the files to be loaded
write(snapshot,'(i8)') (IndSave(1))
nomefile(1)='u_'//trim(adjustl(snapshot))//'.bin'
nomefile(2)='v_'//trim(adjustl(snapshot))//'.bin'
nomefile(3)='w_'//trim(adjustl(snapshot))//'.bin'


!....Read solutions at the time-step 1
call LoadScalarField(u0(0:Ny,0:Nx-1,0:Nzloc-1),nomefile(1))
call LoadScalarField(v0(0:Ny,0:Nx-1,0:Nzloc-1),nomefile(2))
call LoadScalarField(w0(0:Ny,0:Nx-1,0:Nzloc-1),nomefile(3))


!....Compute the Velocity Increments
call EvalSquaredVelIncrementsNormal()

dux1MT = dux1/2.0d0
dux2MT = dux2/2.0d0
duy1MT = duy1/2.0d0
duy2MT = duy2/2.0d0


call MPI_BARRIER(MPI_Comm_World,ierror) 


!....Loop to calculate means
do Iext=2,Nsavings-1

if(nid.eq.0)then
WRITE(*,*)'Processing step: ',IndSave(Iext)
end if

!....Generate names of the files to be loaded
write(snapshot,'(i8)') (IndSave(Iext))
nomefile(1)='u_'//trim(adjustl(snapshot))//'.bin'
nomefile(2)='v_'//trim(adjustl(snapshot))//'.bin'
nomefile(3)='w_'//trim(adjustl(snapshot))//'.bin'



!....Read solutions at the time-step n
call LoadScalarField(u0(0:Ny,0:Nx-1,0:Nzloc-1),nomefile(1))
call LoadScalarField(v0(0:Ny,0:Nx-1,0:Nzloc-1),nomefile(2))
call LoadScalarField(w0(0:Ny,0:Nx-1,0:Nzloc-1),nomefile(3))



!....Compute the Velocity Increments
call EvalSquaredVelIncrementsNormal()

dux1MT = dux1MT + dux1
dux2MT = dux2MT + dux2
duy1MT = duy1MT + duy1
duy2MT = duy2MT + duy2



call MPI_BARRIER(MPI_Comm_World,ierror) 


end do
!End of loop for central savings

call MPI_BARRIER(MPI_Comm_World,ierror) 

if(nid.eq.0)then
WRITE(*,*)'Processing step: ',IndSave(Nsavings)
end if

!....Generate names of the files to be loaded
write(snapshot,'(i8)') (IndSave(Nsavings))
nomefile(1)='u_'//trim(adjustl(snapshot))//'.bin'
nomefile(2)='v_'//trim(adjustl(snapshot))//'.bin'
nomefile(3)='w_'//trim(adjustl(snapshot))//'.bin'


! fileScal='PScal_'//trim(adjustl(snapshot))//'.bin'


!....Read solutions at the time-step n
call LoadScalarField(u0(0:Ny,0:Nx-1,0:Nzloc-1),nomefile(1))
call LoadScalarField(v0(0:Ny,0:Nx-1,0:Nzloc-1),nomefile(2))
call LoadScalarField(w0(0:Ny,0:Nx-1,0:Nzloc-1),nomefile(3))



!....Compute the Velocity Increments
call EvalSquaredVelIncrementsNormal()

dux1MT = dux1MT + dux1/2.0d0
dux2MT = dux2MT + dux2/2.0d0
duy1MT = duy1MT + duy1/2.0d0
duy2MT = duy2MT + duy2/2.0d0


dux1MT = dux1MT / Ninterv
dux2MT = dux2MT / Ninterv
duy1MT = duy1MT / Ninterv
duy2MT = duy2MT / Ninterv

call MPI_BARRIER(MPI_Comm_World,ierror) 


!Average om wall-parallel planes
call MeanInXZ(dux1MT(0:Ny,0:Nx-1,0:Nzloc-1),dux1MTY(0:Ny))
call MeanInXZ(dux2MT(0:Ny,0:Nx-1,0:Nzloc-1),dux2MTY(0:Ny))
call MeanInXZ(duy1MT(0:Ny,0:Nx-1,0:Nzloc-1),duy1MTY(0:Ny))
call MeanInXZ(duy2MT(0:Ny,0:Nx-1,0:Nzloc-1),duy2MTY(0:Ny))
!call EvalUbulk(Ubulk(Nsavings),uhat(0:Ny,0,0))

call MPI_BARRIER(MPI_Comm_World,ierror) 
!.....End of means evaluation


!....Saving Mean Profiles
if(nid.eq.0)then
WRITE(*,*)'**********************************************************'
WRITE(*,*)'**********************************************************' 
WRITE(*,*)'Saving the mean profiles...'
WRITE(*,*)
end if

resfile(1)='VIncr_X.txt'
! resfile(2)='v_mean.bin'
! resfile(3)='w_mean.bin'


call MPI_BARRIER(MPI_Comm_World,ierror) 


FMT1="(2X,E12.5,2X,E12.5,2X,E12.5,2X,E12.5,2X,E12.5,2X,E12.5)"

if(nid.eq.0)then
open(2,file=resfile(1))
write(2,*) '%       y        y+      du_normal_x_sep1     du_normal_y_sep1      du_normal_x_sep2      du_normal_y_sep2     '
do ii=0,Ny
write(2,FMT1) y(Ny-ii),yp(ii),dux1MTY(Ny-ii),duy1MTY(Ny-ii),dux2MTY(Ny-ii),duy2MTY(Ny-ii) 
end do

close(2)

end if





if(nid.eq.0)then
WRITE(*,*)
WRITE(*,*)'Mean profiles Saved.'
end if

call MPI_BARRIER(MPI_Comm_World,ierror)  
!....End of saving mean profiles

 
call transforms_finalization()


call MPI_BARRIER(MPI_Comm_World,ierror) 


!....Deallocate Variables  
call deallocate_variabs()

! deallocate(gacc)
deallocate(nomefile,nomefile1)
deallocate(resfile)


!.....Closing
write(*,*) 'Stop!!! Processore:',nid
call MPI_Finalize(ierror)
      
stop

end program strfun

