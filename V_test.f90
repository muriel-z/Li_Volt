program V_test
use precision
implicit none
integer                                 ::nvx,nvy,nvz !numero de grilla
integer                                 ::pnum,N_iter !numero de particulas, numero de iteraciones
real(np),dimension(:,:,:),allocatable   ::V,V0        !voltaje nuevo y anterior
real(np),dimension(:,:),allocatable     ::V1,V2,V3,V4 !Distribuciones 2D de las paredes
real(np)                                ::dv,Vtop     !Vtop es el valor del voltaje en la tapa
real(np),dimension(:),allocatable       ::Vp,Vpp,m        !Voltaje asociado a cada particula,masa
character,dimension(:),allocatable      ::sym         !coso que esta en el archivo de posiciones
real(np)                                ::d1,d2,d3,dist2,div
real(np)                                ::Rcut, Rcut2
real(np)                                ::boxmin(3),boxmax(3),h(3),r3(3)
integer                                 ::i,j,k,l
real(np),dimension(:,:),allocatable     :: r       ! Posiciones
real(np), dimension(:), allocatable     ::hist
real(np)                                ::dx,x
integer                                 ::ibin,nbin
character(len=6)                        ::aux
character(len=15)                       ::filename

!Este programa tiene que agarrar una configuracion inicial y calcular el potencial
!Del programa de python falta agregar el suavizado del campo y la funcion de radio para calcular el potencial en el voxel

!Dimensiones de la caja
boxmin(:)=[0._np,0._np,0._np]
boxmax(:)=[100._np,100._np,400._np]

!Variables enteras
nvx = 10
nvy = 10
nvz = 40
N_iter = 2000

h(1)=(boxmax(1)-boxmin(1))/(nvx-1)
h(2)=(boxmax(2)-boxmin(2))/(nvy-1)
h(3)=(boxmax(3)-boxmin(3))/(nvz-1)

!Variables del potencial
allocate(V(nvx,nvy,nvz))
allocate(V0(0:nvx+1,0:nvy+1,0:nvz+1))

V0(:,:,:) = 0._np
V(:,:,:) = 0._np

!Condicion de borde superior
Vtop = 10._np
V0(:,:,1) = 0._np
V0(:,:,nvz) = Vtop
          
!Doy el valor de la matriz del voltaje en t=0.Resuelta por Laplace

!Itero para encontrar la solucion suave

write(*,*) 'Empezamo'

do l= 1,N_iter
    write(*,*) real(l,np)/real(N_iter,np)*100._np
  
    do i = 1,nvx
        do j = 1,nvy
            do k = 1,nvz-1 
                v(i,j,k) = (V0(i+1,j,k)+V0(i-1,j,k)+V0(i,j+1,k)+V0(i,j-1,k)+V0(i,j,k+1)+V0(i,j,k-1))/6._np
            enddo
        enddo
    enddo
 
    ! Contorno
    ! TODO
    ! do l=1,npart ! sobre las particulas metalicas
    !   ri = int((r(l,1) - boxmin(1)) / h(1)) + 1
    !   rj = int((r(l,2) - boxmin(2)) / h(2)) + 1
    !   rk = int((r(l,3) - boxmin(3)) / h(3)) + 1
    !   V(ri,rj,rk) = 0._np
    !   V(ri+1,rj+1,rk+1) = 0._np
    !   V(ri+1,rj,rk) = 0._np
    !   V(ri+1,rj+1,rk) = 0._np
    ! enddo
    V(:,:,1) = 0._np
    V(:,:,nvz) = Vtop

    ! Avance
    V0(1:nvx,1:nvy,1:nvz)=v(:,:,:)

    ! PBC
    V0(nvx+1,:,:) = V(1,:,:)
    V0(:,nvy+1,:) = V(:,1,:)
    V0(0,:,:) = V(nvx,:,:)
    V0(:,0,:) = V(:,nvy,:)

    ! No creo que haga falta
    V0(:,:,0) = 0._np
    ! V0(:,:,nvz+1) = Vtop
        
enddo


!Paso los valores del potencial a un archivo para graficarlos
open(12,file='Vyz.dat',status='replace')
write(12,10) 'j','k','Vyz'
do j = 1,nvy
    do k = 1,nvz
        write(12,11)real(j,np),real(k,np),V0(nvx/2,j,k)
    enddo
    write(12,7)''
enddo
close(12)

open(13,file='Vxz.dat',status='replace')
write(13,10) 'i','k','Vxz'
do i = 1,nvx
    do k = 1,nvz
        write(13,11)real(i,np),real(k,np),V0(i,nvy/2,k)
    enddo
    write(13,7)''
enddo
close(13)

write(*,*) 'Ya terminó de calcular el potencial inicial V0'

!Llamo a las posiciones iniciales

pnum = 2000

allocate(r(pnum,3))
allocate(m(pnum))
allocate(sym(pnum))

open(14,File='pos_inicial.xyz')
read(14,*)
read(14,*)
do i=1,pnum
    read(14,*) sym(i),r(i,:),m(i)
end do
close(14)

write(*,*) 'Ya leyó las posiciones iniciales'

!Ahora tengo que asignar a cada posición un valor de voltaje
allocate(Vp(0:pnum))
allocate(vpp(0:pnum))


!uso la función de distribución radial
Rcut = sqrt(2._np)
Rcut2 = Rcut*Rcut

open(16,file='Vp.dat',status='replace')
write(16,8) 'l','r(l,1)','r(l,2)','r(l,3)','Vp'


Vp(:) = 0._np
Vpp(:) = 0._np
do l = 1,pnum
    div = 0._np
    write(*,*) real(l,np)/real(pnum,np)*100._np
    
    r3(:)=r(l,:)
    Vpp(l)=trilinear_interpolation(r3, V, boxmin, h)

    write(16,9)real(l,np),r(l,1),r(l,2),r(l,3),Vpp(l)
enddo
close(16)



!************* Armo el histograma del voltaje de cada particula***********

! nbin = 11
! write(aux,'(F5.2)') Vtop
! write(*,*) 'aux=', aux,Vtop
! filename = 'hist_Vpp' //trim(adjustl(aux))// '.d'
! write(*,*) 'filename=',filename
!   
! allocate( hist(nbin) )
! dx = 10._np/real(nbin,np)
! hist(:) = 0
! open(21,file='Vpp_1.d',status='replace')
! do l = 1,pnum
!     write(21,*) l, Vpp(l)
!     ibin = int(Vpp(l)) + 1
!     hist(ibin) = hist(ibin) + 1
! enddo
! close(21)
! open(22,file=filename,status='replace')
! do i = 1,nbin
!     write(22,*) (real(i-1,np) + 0.5_np)*dx ,hist(i)
! enddo
! close(22)
!
! write(*,*) 'Ya asignó un valor de voltaje a cada partícula'

7 FORMAT (A3)
8 FORMAT (5(3X,A15))
9 FORMAT (5(2X,ES17.9)) 
                
10 FORMAT (3(3X,A15))
11 FORMAT (3(2X,ES17.9)) 
   

contains

function trilinear_interpolation(r, f, boxmin, h) result(interp_val)
  real(np), intent(in) :: r(3), f(:,:,:)
  real(np), intent(in) :: boxmin(3), h(3)
  real(np) :: interp_val
  
  integer :: i, j, k, i1, j1, k1
  real(np) :: x, y, z, dx, dy, dz
  
  i = int((r(1) - boxmin(1)) / h(1)) + 1
  j = int((r(2) - boxmin(2)) / h(2)) + 1
  k = int((r(3) - boxmin(3)) / h(3)) + 1
  
  ! Esto deberia devolver error al salirse de la grilla
  ! if (i < 1) then
  !   i = 1
  ! elseif (i > size(f, 1) - 1) then
  !   i = size(f, 1) - 1
  ! endif
  !
  ! if (j < 1) then
  !   j = 1
  ! elseif (j > size(f, 2) - 1) then
  !   j = size(f, 2) - 1
  ! endif
  !
  ! if (k < 1) then
  !   k = 1
  ! elseif (k > size(f, 3) - 1) then
  !   k = size(f, 3) - 1
  ! endif
  
  i1 = i + 1
  j1 = j + 1
  k1 = k + 1
  
  x = (r(1) - boxmin(1)) / h(1) - float(i - 1)
  y = (r(2) - boxmin(2)) / h(2) - float(j - 1)
  z = (r(3) - boxmin(3)) / h(3) - float(k - 1)
  
  dx = 1.0 - x
  dy = 1.0 - y
  dz = 1.0 - z
  
  interp_val = dx*dy*dz*f(i,j,k) + x*dy*dz*f(i1,j,k) &
             + dx*y*dz*f(i,j1,k) + x*y*dz*f(i1,j1,k) &
             + dx*dy*z*f(i,j,k1) + x*dy*z*f(i1,j,k1) &
             + dx*y*z*f(i,j1,k1) + x*y*z*f(i1,j1,k1)
end function trilinear_interpolation


end program V_test
