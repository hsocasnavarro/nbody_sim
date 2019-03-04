  Implicit None

  Integer, Parameter :: RealBytes=2 
  Integer, Parameter :: dp=8
  Real (kind=8), Parameter :: G_grav=6.67408d-11, cc=2.99792458d8
!  Real (kind=8), Parameter :: G_grav=1.0, cc=2.99792458d8
  Integer :: NN
  Real (kind=8), Dimension(:), Allocatable :: m, record
  Real (kind=8), Dimension(:,:), Allocatable :: x, v, v_0
  Real (kind=8), Dimension(:,:), Allocatable :: x_slice, v_slice, a_slice
  Real (kind=8) :: eps, eps2, dt=1e6*365.25*24.*3600 ! Myr in sec
  Real (kind=8) :: Norm_x=1e9*365.25*24.*3600.*cc
  Real (kind=8), Dimension(3) :: CubeSize ! Glyr in m (half-length of cube side)
  Real (kind=8) :: Norm_t=1e6*365.25*24.*3600. ! Myr in s
  Real (kind=8) :: Norm_m=1d9*1.989d30 ! 1e9 solar masses in kg
  Real (kind=8), Dimension(3) :: xp
  Integer :: myrank, nprocs, nslaves, nperproc, islave
  Integer, Dimension(:), Allocatable :: i0, i1, NN_slice
  Integer :: ii, i, j, status, idir, it, nt, ierr, TaskComing, recl
  Integer :: iframe, iframeskip, iperiodic, myi0, myi1
  Real (kind=8), Dimension(25) :: Param_vector ! length must match MPI_SEND command
  Real (kind=8) :: kk, t, Geff, a, F, dist2, dist32, maxvel

  Include 'mpif.h'
  Integer :: stat(MPI_STATUS_SIZE)
  Integer :: My_MPI_Real
  !
  myrank=0
  nprocs=1
  !
  Call MPI_INIT(status)
  Call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, status)
  Call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, status)
  My_MPI_Real=MPI_Real
  If (Kind(1.0d0) .eq. 8) then
     My_MPI_Real=MPI_Real8
  Else
     Print *,'Something is wrong. Please check MPI_Real'
     Stop
  End if
  !
  If(myrank.eq.0) then ! I'm the master. Do all initializations
     If (nprocs .lt. 2) then
        Print *,'At least 2 processes are required for parallel run'
        Stop
     Endif
     !
     ! Simulation parameters (read from file)
     !
     Open (30,File='input.bin',Form='Unformatted',Access='Direct', Recl=25*RealBytes)
     Allocate(record(25))
     Read (30, rec=1) record
     NN=record(1)+0.5; CubeSize(1)=record(2); Norm_t=record(3); Norm_m=record(4)
     eps=record(5); dt=record(6); nt=record(7)+0.5; iframeskip=record(8)+0.5
     iperiodic=record(9)+0.5; CubeSize(2)=record(10)*CubeSize(1)
     CubeSize(3)=record(11)*CubeSize(1)
     Norm_x=CubeSize(1)
     If (iframeskip .eq. 0) iframeskip = 1
     Deallocate(record)
     recl=NN+NN*3+NN*3
     Allocate(record(recl))
     Close(30)
     Print *,'Input file read. Particles:',NN,' Frames:',nt
     !
     Geff=G_grav*Norm_t/Norm_x*Norm_t/Norm_x*Norm_m/Norm_x ! Gravity in normalized units
     Allocate (x(NN,3), v(NN,3), m(NN), v_0(NN,3))
     !
     ! Simulation starting conditions (read from file)
     !
     Open (30,File='input.bin',Form='Unformatted',Access='Direct', Recl=recl*RealBytes)
     Read (30, rec=2) record
     ii=1
     m(1:NN)=record(ii:ii+NN-1)
     ii=ii+NN
     Do idir=1, 3
        x(1:NN,idir)=record(ii:ii+NN-1)
        ii=ii+NN
     End do
     Do idir=1, 3
        v(1:NN,idir)=record(ii:ii+NN-1)
        ii=ii+NN
     End do
     Read (30, rec=1) record ! keep header record to write it on outputfile
     Close(30)
     !
     nslaves=nprocs-1
     nperproc=NN/nslaves
     if (nperproc*nslaves .lt. NN) nperproc=nperproc+1
     Allocate(i0(nslaves),i1(nslaves),NN_slice(nslaves))
     !
     ii=1
     Do islave=1, nslaves
        i0(islave)=ii
        i1(islave)=Min(ii+nperproc-1, NN)
        NN_slice(islave)=i1(islave)-i0(islave)+1
        ii=i1(islave)+1
        Print *,'Proc:',islave,'Number of particles:',NN_slice(islave),nperproc
     End do
     Allocate (x_slice(nperproc,3), v_slice(nperproc,3))
     !
     ! Prepare output file
     !
     open(30,file='output.bin',form='unformatted',access='direct',recl=recl*RealBytes)
     Write (30, rec=1) record ! write header record
  End if ! if myrank .eq. 0

  !
  ! Start simulation
  !
  If(myrank.eq.0) then ! I'm the master. Start simulation
     ! Use normalized units
     dt=dt/Norm_t 
     x=x/Norm_x 
     v=v/Norm_x*Norm_t 
     m=m/Norm_m 
     eps=eps/Norm_x 
     CubeSize(:)=CubeSize(:)/Norm_x
     ! Remember initial velocity (will need it for leapfrog kick)
     v_0(:,:)=v(:,:)
     ! Pack parameters to transfer to slaves
     Param_vector(1)=NN
     Param_vector(4)=dt
     Param_vector(5)=Geff
     Param_vector(6)=nperproc
     Param_vector(7)=eps**2
     Param_vector(8)=iperiodic
     Param_vector(9:11)=CubeSize(1:3)
     !
     iframe=1
     Do it=1, nt
        Do islave=1, nslaves ! Send job to slaves
           If (i0(islave) .le. NN) then ! If this slave is active
              Param_vector(2)=i0(islave)
              Param_vector(3)=NN_slice(islave)
              !
              x_slice(1:i1(islave)-i0(islave)+1,:)=x(i0(islave):i1(islave),:)
              v_slice(1:i1(islave)-i0(islave)+1,:)=v(i0(islave):i1(islave),:)
              TaskComing=1
              Call MPI_Send(TaskComing, 1, MPI_INTEGER, islave, 1, MPI_COMM_WORLD, ierr) ! Task is coming
              Call MPI_Send(Param_vector, 25, My_MPI_REAL, islave, 2, MPI_COMM_WORLD, ierr) ! Vector of parameters
              Call MPI_Send(m, NN, My_MPI_REAL, islave, 3, MPI_COMM_WORLD, ierr) 
              Call MPI_Send(x, NN*3, My_MPI_REAL, islave, 4, MPI_COMM_WORLD, ierr)
              Call MPI_Send(x_slice, nperproc*3, My_MPI_REAL, islave, 5, MPI_COMM_WORLD, ierr) 
              Call MPI_Send(v_slice, nperproc*3, My_MPI_REAL, islave, 6, MPI_COMM_WORLD, ierr)
           End if
        End do ! do in islave
           !
        Do islave=1, nslaves ! Collect results from slaves
           If (i0(islave) .le. NN) then ! If this slave is active
              Call MPI_Recv(x_slice, nperproc*3, My_MPI_REAL, islave, 7, MPI_COMM_WORLD, stat, ierr)
              Call MPI_Recv(v_slice, nperproc*3, My_MPI_REAL, islave, 8, MPI_COMM_WORLD, stat, ierr)
              x(i0(islave):i1(islave),:)=x_slice(1:i1(islave)-i0(islave)+1,:)
              v(i0(islave):i1(islave),:)=v_slice(1:i1(islave)-i0(islave)+1,:)
           End if
        End do
        ! If this is the first iteration, give initial "kick" to velocity
        ! to make it v_1/2=v_0+a(x_0)*dt/2 and don't update positions
        If (it .eq. 1) then
           v(:,:)=v_0(:,:)+(v(:,:)-v_0(:,:))/2.
        Else
        ! Update positions
           x(:,:)=x(:,:)+v(:,:)*dt
        End if
        ! Write to file (convert back to physical units)
        If (MODULO(it,iframeskip) .eq. 0) then
           ii=1
           record(ii:ii+NN-1)=m(1:NN)*Norm_m
           ii=ii+NN
           Do idir=1, 3
              record(ii:ii+NN-1)=x(1:NN,idir)*Norm_x
              ii=ii+NN
           End do
           Do idir=1, 3
              record(ii:ii+NN-1)=v(1:NN,idir)*Norm_x/Norm_t
              ii=ii+NN
           End do
           Write (30, rec=iframe+1) record
           iframe=iframe+1
        End if
        maxvel=maxval(abs(v))
        print *,'it=',it,' max vel=',maxvel
!        Print *,'Warning! Simulation done with G == 1'
     End do ! do in it
     ! Finished. Notify slaves that job is done
     Do islave=1, nslaves
        TaskComing=0
        Call MPI_Send(TaskComing, 1, MPI_INTEGER, islave, 1, MPI_COMM_WORLD, ierr)
     End do
     ! Close file and finish stuff
     Call MPI_FINALIZE(status)
     Close(30)
     Print *,'DONE'
     Stop
  End if ! if myrank .eq. 0
  !
  If (myrank .gt. 0) then ! If I'm a slave
     Do while (.True.) ! loop waiting for instructions
        Call MPI_Recv(TaskComing, 1, MPI_INTEGER, 0, 1, MPI_COMM_WORLD, stat, ierr)
        If (Taskcoming .eq. 0) Exit ! Master says no more tasks. Break loop
        ! Otherwise collect data
        Call MPI_Recv(Param_vector, 25, My_MPI_REAL, 0, 2, MPI_COMM_WORLD, stat, ierr)
        NN=Param_vector(1)+.5 ! the .5 are added to properly convert float to int
        myi0=Param_vector(2)+.5
        myi1=myi0+Param_vector(3)-1+.5
        dt=Param_vector(4)
        Geff=Param_vector(5)
        nperproc=Param_vector(6)+.5
        eps2=Param_vector(7)
        iperiodic=Param_vector(8)+.5
        CubeSize(1:3)=Param_vector(9:11)
        If (.not. Allocated(m)) & ! If first time, intialize
             Allocate (m(NN), x(NN,3), x_slice(nperproc,3), v_slice(nperproc,3), a_slice(nperproc,3))
        
        Call MPI_Recv(m, NN, My_MPI_REAL, 0, 3, MPI_COMM_WORLD, stat, ierr)
        Call MPI_Recv(x, NN*3, My_MPI_REAL, 0, 4, MPI_COMM_WORLD, stat, ierr)
        Call MPI_Recv(x_slice, nperproc*3, My_MPI_REAL, 0, 5, MPI_COMM_WORLD, stat, ierr)
        Call MPI_Recv(v_slice, nperproc*3, My_MPI_REAL, 0, 6, MPI_COMM_WORLD, stat, ierr)
        
        ! Do processing
        if (iperiodic .eq. 0) then
           Call do_task ! No boundary conditions
        Else if (iperiodic .eq. 1) then
           Call do_task_periodic ! Periodic boundary conditions
        Else
           Print *,'Incorrect switch for boundary conditions:',iperiodic
           Stop
        End if
        
        ! Return updated v and x (but x is only updated if it hits boundaries)
        Call MPI_Send(x_slice, nperproc*3, My_MPI_REAL, 0, 7, MPI_COMM_WORLD, ierr)
        Call MPI_Send(v_slice, nperproc*3, My_MPI_REAL, 0, 8, MPI_COMM_WORLD, ierr)
        !
     End do ! While True wait for instructions
  End if ! my rank .gt. 0
  !
  !
Contains
  Subroutine do_task
    ! Leapfrog algorithm
    ! https://en.wikipedia.org/wiki/Leapfrog_integration
    !
    !  array v is assumed to be v(i-1/2)
    !
    Do ii=1, myi1-myi0+1
       i=ii+myi0-1 ! index in global array
       Do j=1, i-1
          dist2=(x(i,1)-x(j,1))**2+(x(i,2)-x(j,2))**2+(x(i,3)-x(j,3))**2
          dist32=(dist2+eps2)**(3./2.)
          Do idir=1, 3 ! 3 dimensions
             a=-Geff*m(j)*(x(i,idir)-x(j,idir))/dist32 ! F/m(i)
             v_slice(ii,idir)=v_slice(ii,idir)+a*dt ! update v. Now it's v(i+1/2)
          End do
       End do
       Do j=i+1, NN
          dist2=(x(i,1)-x(j,1))**2+(x(i,2)-x(j,2))**2+(x(i,3)-x(j,3))**2
          dist32=(dist2+eps2)**(3./2.)
          Do idir=1, 3 ! 3 dimensions
             a=-Geff*m(j)*(x(i,idir)-x(j,idir))/dist32 ! F/m(i)
             v_slice(ii,idir)=v_slice(ii,idir)+a*dt ! update v. Now it's v(i+1/2)
          End do
       End do
    End do
    
    Return
  End Subroutine do_task
  
  Subroutine do_task_periodic
    ! Leapfrog algorithm
    ! https://en.wikipedia.org/wiki/Leapfrog_integration
    !
    !  array v is assumed to be v(i-1/2)
    !
    ! Periodic boundary condition. Shift particle ii to center of box
    !  (i.e., if a distance is greater than CubeSize, put it in symmetric
    !   position)
    Do ii=1, myi1-myi0+1
       Do idir=1, 3 ! If particle leaves the box, come in from the other side
          If (x_slice(ii,idir) .gt. CubeSize(idir)) x_slice(ii,idir)=x_slice(ii,idir)-2.*CubeSize(idir)
          If (x_slice(ii,idir) .lt. -CubeSize(idir)) x_slice(ii,idir)=x_slice(ii,idir)+2.*CubeSize(idir)
       End do
       i=ii+myi0-1 ! index in global array
       Do j=1, i-1
          xp=x(j,:)
          Do idir=1, 3 ! Reverse j-particle if further than CubeSize
             If (xp(idir)-x(i,idir) .gt. CubeSize(idir)) &
                  xp(idir)=x(i,idir)-2*CubeSize(idir)+(xp(idir)-x(i,idir))
             If (x(i,idir)-xp(idir) .gt. CubeSize(idir)) &
                  xp(idir)=x(i,idir)+2*CubeSize(idir)-(x(i,idir)-xp(idir))
          End do
          dist2=(x(i,1)-xp(1))**2+(x(i,2)-xp(2))**2+(x(i,3)-xp(3))**2
          dist32=(dist2+eps2)**(3./2.)
          Do idir=1, 3 ! 3 dimensions
             a=-Geff*m(j)*(x(i,idir)-xp(idir))/dist32 ! F/m(i)
             v_slice(ii,idir)=v_slice(ii,idir)+a*dt ! update v. Now it's v(i+1/2)
          End do
       End do
       Do j=i+1, NN
          xp=x(j,:)
          Do idir=1, 3 ! Reverse j-particle if further than 1.
             If (xp(idir)-x(i,idir) .gt. CubeSize(idir)) &
                  xp(idir)=x(i,idir)-2*CubeSize(idir)+(xp(idir)-x(i,idir))
             If (x(i,idir)-xp(idir) .gt. CubeSize(idir)) &
                  xp(idir)=x(i,idir)+2*CubeSize(idir)-(x(i,idir)-xp(idir))
          End do
          dist2=(x(i,1)-xp(1))**2+(x(i,2)-xp(2))**2+(x(i,3)-xp(3))**2
          dist32=(dist2+eps2)**(3./2.)
          Do idir=1, 3 ! 3 dimensions
             a=-Geff*m(j)*(x(i,idir)-xp(idir))/dist32 ! F/m(i)
             v_slice(ii,idir)=v_slice(ii,idir)+a*dt ! update v. Now it's v(i+1/2)
          End do
       End do
    End do

    Return
  End Subroutine do_task_periodic
End program
