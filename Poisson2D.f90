program main
implicit none

real:: xsize, ysize, dx, dy
integer:: Nx, Ny, N, i, j, gi, k
real, allocatable:: A(:,:), b(:)
real, allocatable:: Minv(:,:);
real, allocatable:: r0(:), z0(:), x0(:), p0(:)
real, allocatable:: r1(:), z1(:), x1(:), p1(:)
real:: eps
integer:: maxsteps
real:: alpha
real:: beta
real, allocatable:: phi(:,:)
real:: error

xsize = 10
ysize = 10
Nx = 101
Ny = 51
dx = xsize/(Nx-1)
dy = ysize/(Ny-1)
N = Nx*Ny
eps = 1e-3
maxsteps = 1e3

allocate(A(N,N),b(N))
allocate(Minv(N,N))
allocate(r0(N),z0(N),x0(N),p0(N))
allocate(r1(N),z1(N),x1(N),p1(N))

allocate(phi(Ny,Nx))

do j=1,Nx
do i=1,Ny
gi=(j-1)*Ny+i
if(j==1.or.j==Nx.or.i==1.or.i==Ny) then
    A(gi,gi)=1
    b(gi)=0
else
    A(gi,gi-Ny)=1.0/(dx**2)
    A(gi,gi-1)=1.0/(dy**2)
    A(gi,gi)=-2.0/(dx**2)-2.0/(dy**2)
    A(gi,gi+1)=1.0/(dy**2)
    A(gi,gi+Ny)=1.0/(dx**2)
    
    b(gi)=-1
end if
end do
end do

Minv = 0
do i=1,N
Minv(i,i)=1.0/A(i,i)
end do

x0 = 0

r0 = b-MATMUL(A,x0)
z0 = MATMUL(Minv,r0)
p0 = z0

open(2,file="StepErrorPCG.log",status="replace")
do k=1,maxsteps
alpha = dot_product(r0,z0)/dot_product(p0,MATMUL(A,p0))
x1 = x0+alpha*p0
r1 = r0-alpha*MATMUL(A,p0)
error = sqrt(dot_product(r1,r1)/N)
if(error<eps) then
    print*, "Reach tolerance: ", eps
    exit
else
    print*, k, error
    write(2,'(1000(1pe12.3))'), error
end if
z1 = MATMUL(Minv,r1)
beta = dot_product(z1,r1)/dot_product(z0,r0)
p1 = z1 + beta*p0

r0 = r1
z0 = z1
x0 = x1
p0 = p1
end do
close(2)

do j=1,Nx
do i=1,Ny
gi = (j-1)*Ny+i
phi(i,j)=x1(gi)
end do
end do

open(1,file='Poisson2D.dat',status='replace')
do i=1,Ny
write(1,'(1000(1pe12.3))'),phi(i,:)
end do
close(1)

deallocate(A,b,Minv,r0,r1,z0,z1,x0,x1,p0,p1)

end program
