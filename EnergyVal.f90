program energy
implicit none
  integer(8) n, i
  real(8) beta, lambda, kappa, D, a, En
  real(8), parameter :: hbar=1.05457180d-34
  real(8), parameter :: mu=1.62747d-27
  real(8), parameter :: Ese=1.60217662d-19
  real(8), parameter :: Lse=1.0d-9
  open(1,file='in.dat',action='read')
   open(2,file='eigval.dat',action='write')
  read(1,*) D
  read(1,*) a
   close(1)
  D=D*Ese
  a=a*Lse
  beta=sqrt(2.0d0*mu*D)
  lambda=a*beta
  kappa=hbar/lambda
  n=floor(1.0d0/kappa-0.50d0)     ! determine the number of the energy states
  do i=0, n
    En=-D*(1.0d0-kappa*(i+0.50d0))**2  ! determine the eigenvalue
    En=En/Ese;  ! energy unit conversion, back to eV
  write(2,*) En
  end do
close(2)
end program energy
