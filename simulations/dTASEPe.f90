!==============================================================================
!
!   MODEL
!   -----
!   Totally asymmetric simple exclusion process with open boundary conditions,
!   elongated particles and site-dependent hopping rates.
!
!   PROGRAM
!   -------
!   The program computes elongation rates relative to the initiation rate from 
!   the Ribo-seq A-site density profile using NLOPT package and the third order 
!   of the power series method. 
!
!==============================================================================

program dTASEPe
  
  use mt19937
  implicit none

  interface
    subroutine mc(L,ll,alpha,omega,tau,iter,rho,current)
      integer(kind=4),intent(in) :: L,ll
      real(kind=8),intent(in) :: alpha,omega(1:L)
      integer(kind=4),intent(inout) :: tau(1:L)
      integer(kind=8),intent(in) :: iter     
      real(kind=8),intent(out) :: rho(1:L),current(0:L)
    end subroutine mc
  end interface

  ! program parameters
  integer(kind=4) :: err11,err22,err33,err44,err55,stat22,i,L,ll
  integer(kind=8) :: iter
  integer(kind=4),allocatable :: tau(:)
  character(len=80) :: inputfile_rates, outputfile_current, &
    & outputfile_density, outputfile_log
  logical :: eof,found
  real(kind=8),allocatable :: rho(:),current(:),omega(:)
  real(kind=8) :: t1,t2,omegai,tden1,tden2,rel,alpha
  external fnc3

!======================================================================
! INPUT
!======================================================================

  call cpu_time(t1)

  open(unit=11,file='dTASEPe.dat',status='OLD',iostat=err11)
  read(11,*) ll                 ! particle size
  read(11,*) alpha              ! initiation rate
  read(11,*) iter               ! number of iterations for MC (Gillespie algorithm)
  read(11,*) inputfile_rates    ! file with hopping rates omega(1), ... , omega(L)
  read(11,*) outputfile_density ! file with local density rho(1), ... , rho(L)
  read(11,*) outputfile_current ! file with current current(0), ... , current(L)
  read(11,*) outputfile_log     ! log file 
  close (11)

!======================================================================
! MAIN PROGRAM
!======================================================================

  ! begin measuring time
  call cpu_time(t1)
  
  ! write to log file
  open(55,file=outputfile_log,status='REPLACE',iostat=err55)
  write(55,'(A22)') 'INFO: Program started.'
  close(55)
  
  ! initial seed for random number generator
  call init_genrand(4357)
  
  ! read file with hopping rates (to get lattice size)
  !---------------------------------------------------
  open(22,file=inputfile_rates,status='OLD',iostat=err22)
  i=0
  eof=.FALSE.
  do while (eof.EQV..FALSE.)
    read(22,*,iostat=stat22) omegai
    if (stat22.EQ.0) then       
      i=i+1
    else
      eof=.TRUE.
    endif
  enddo
  close(22)

  ! set lattice size
  !--------------
  L=i
  
  ! read file with hopping rates (to get hopping rates)
  !----------------------------------------------------
  allocate(omega(1:L))
  open(22,file=inputfile_rates,status='OLD',iostat=err22)
  do i=1,L
    read(22,*) omega(i)
  enddo
  close(22)
  
  ! run Monte Carlo simulation (to reach the steady state)
  !-------------------------------------------------------
  allocate(tau(1:L))
  tau=0
  allocate(rho(1:L),current(0:L))
  rho=0.0d0
  current=0.0d0
  found=.FALSE.
  do while (found.EQV..FALSE.)
    tden1=sum(rho)/dble(L)
    call mc(L,ll,alpha,omega,tau,100_8,rho,current)
    tden2=sum(rho)/dble(L)
    if (tden2.gt.epsilon(1.0d0)) then
      rel=abs(tden2-tden1)/tden2
      if (rel.lt.1.d-3) then
        found=.TRUE.
      endif
    endif
  enddo
    
  ! run Monte Carlo simulation (in the steady state)
  !--------------------------------------------------
  rho=0.0d0
  current=0.0d0
  call mc(L,ll,alpha,omega,tau,iter,rho,current)
  deallocate(tau)

  ! write output
  !--------------
  open(33,file=outputfile_density,status='REPLACE',iostat=err33)
  open(44,file=outputfile_current,status='REPLACE',iostat=err44)
  write(44,*) current(0)
  do i=1,L
    write(33,*) rho(i)
    write(44,*) current(i)
  enddo
  close(33)
  close(44)
  deallocate(rho,current)
  
  ! stop measuring time
  call cpu_time(t2)
  
  ! write to log file
  open(55,file=outputfile_log,status='OLD',position='APPEND',iostat=err55)
  write(55,'(A39,F9.3,A9)') 'INFO: Program finished successfully in ',t2-t1,' seconds.'
  close(55)

end program dTASEPe

!======================================================================
! SUBROUTINES
!======================================================================

subroutine mc(L,ll,alpha,omega,tau,iter,rho,current)

  use mt19937
  implicit none

  integer(kind=4),intent(in) :: L,ll
  real(kind=8),intent(in) :: alpha,omega(1:L)
  integer(kind=4),intent(inout) :: tau(1:L)
  integer(kind=8),intent(in) :: iter
  real(kind=8),intent(out) :: rho(1:L),current(0:L)

  ! local variables
  integer(kind=8) :: i,k
  real(kind=8) :: r1,r2,a0,atest,t,dt
  real(kind=8),allocatable :: a(:)

  ! initial propensity function
  allocate(a(0:L))
  if (sum(tau(1:min(ll,L))).EQ.0) then
    a(0)=alpha
  else
    a(0)=0.0d0
  endif
  do i=1,L
    if (i+ll.le.L) then
      a(i)=omega(i)*dble(tau(i))*dble(1-tau(i+ll))
    else
      a(i)=omega(i)*dble(tau(i))
    endif
  enddo

  t=0.0d0
  rho=0.0d0
  current=0.0d0

  do k=1,iter*L

    a0=sum(a)
   
    r1=grnd()
    r2=grnd()

    ! chooses time until the next event
    dt=-dlog(dble(r1))/a0
    t=t+dt

    ! density profile
    rho=rho+dble(tau)*dt

    ! chooses which event happens next
    i=0
    atest=a(0)
    do while (r2.GT.atest/a0)
      i=i+1
      atest=atest+a(i)
    enddo

    ! updates the propensity function
    a(i)=0.0d0
    
    ! new particles is placed at site 1
    if (i.EQ.0) then
    
      if (1+ll.le.L) then
        a(1)=omega(1)*dble(1-tau(1+ll))
      else
        a(1)=omega(1)
      endif
      tau(1)=1
    
    ! moves particle to the right
    elseif ((i.ge.1).and.(i.le.L-1)) then
    
      ! checks if the moved particle can move again
      if (i+ll+1.le.L) then
        a(i+1)=omega(i+1)*dble(1-tau(i+1+ll))
      else
        a(i+1)=omega(i+1)
      endif
      
      ! checks if new particle can come onto the lattice
      if (i.eq.ll) then
        a(0)=alpha
      endif
        
      ! checks if the trailing particle can move    
      if (i.ge.ll+1) then
        a(i-ll)=omega(i-ll)*dble(tau(i-ll))
      endif
        
      tau(i)=0
      tau(i+1)=1
      
    ! removes particle from the lattice
    else
    
      ! checks if the trailing particle can move    
      if (i.ge.ll+1) then
        a(i-ll)=omega(i-ll)*dble(tau(i-ll))
      endif
      tau(L)=0
      
    endif
    current(i)=current(i)+1.0d0
  enddo
  deallocate(a)
  rho=rho/t
  current=current/t
  return
end subroutine mc