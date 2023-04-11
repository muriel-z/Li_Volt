program V_test
use precision
implicit none
integer                                 ::nvx,nvy,nvz !numero de voxels
integer                                 ::pnum,N_iter !numero de particulas, numero de iteraciones
real(np),dimension(:,:,:),allocatable   ::V,V0        !voltaje nuevo y anterior
real(np),dimension(:,:),allocatable     ::V1,V2,V3,V4 !Distribuciones 2D de las paredes
real(np)                                ::dv,Vtop     !Vtop es el valor del voltaje en la tapa
real(np),dimension(:),allocatable       ::Vp,Vpp,m        !Voltaje asociado a cada particula,masa
character,dimension(:),allocatable      ::sym         !coso que esta en el archivo de posiciones
real(np)                                ::d1,d2,d3,dist2,div
real(np)                                ::Rcut, Rcut2
integer                                 ::i,j,k,l
real(np),dimension(:,:),allocatable     :: r       ! Posiciones
real(np), dimension(:), allocatable     ::hist
real(np)                                ::dx,x
integer                                 ::ibin,nbin
character(len=6)                        ::aux
character(len=15)                       ::filename

!Este programa tiene que agarrar una configuracion inicial y calcular el potencial
!Del programa de python falta agregar el suavizado del campo y la funcion de radio para calcular el potencial en el voxel

!Variables enteras
nvx = 100
nvy = 100
nvz = 400
N_iter = 2000

!Variables del potencial
allocate(V(0:nvx,0:nvy,0:nvz))
allocate(V0(0:nvx,0:nvy,0:nvz))

V0(:,:,:) = 0._np
V(:,:,:) = 0._np

!Doy el valor de la matriz del voltaje en t=0.Resuelta por Laplace

!Condicion de borde superior
Vtop = 10._np
V0(:,:,0) = Vtop

!Itero para encontrar la solucion suave

write(*,*) 'Empezamo'

do l= 1,N_iter
    write(*,*) real(l,np)/real(N_iter,np)*100._np
    do i = 1,nvx-1
        do j = 1,nvy-1
            do k = 1,nvz-1 
                v0(i,j,k) = (V0(i+1,j,k)+V0(i-1,j,k)+V0(i,j+1,k)+V0(i,j-1,k)+V0(i,j,k+1)+V0(i,j,k-1))/6._np
            enddo
        enddo
    enddo
enddo


!Faltaria revisar las condiciones de borde

!Paso los valores del potencial a un archivo para graficarlos
open(12,file='Vyz.dat',status='replace')
write(12,10) 'j','k','Vyz'
do j = 0,nvy
    do k = 0,nvz
        write(12,11)real(j,np),real(k,np),V0(nvx/2,j,k)
    enddo
    write(12,7)''
enddo
close(12)

open(13,file='Vxz.dat',status='replace')
write(13,10) 'i','k','Vxz'
do i = 0,nvx
    do k = 0,nvz
        write(13,11)real(i,np),real(k,np),V0(i,nvy/2,k)
    enddo
    write(13,7)''
enddo
close(13)

V(:,:,:) = v0(:,:,:)

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
    do i = 1,nvx
        do j = 1,nvy
            do k = 1,nvz
                !Ver esto
                d1 = r(l,1)-real(i,np)
                d2 = r(l,2)-real(j,np)
                d3 = r(l,3)-real(k,np)
                dist2 = d1*d1+d2*d2+d3*d3
                if (dist2.le.Rcut2) then
                !Guardo los puntos de la malla que contribuyen al voltaje de la partícula
                    div = div + 1._np
                    Vp(l) =Vp(l) + V(i,j,k)
                endif
            enddo
        enddo
    enddo
    Vpp(l) = Vp(l)/div
    write(16,9)real(l,np),r(l,1),r(l,2),r(l,3),Vpp(l)
enddo
close(16)



!************* Armo el histograma del voltaje de cada particula***********

nbin = 11
write(aux,'(F5.2)') Vtop
write(*,*) 'aux=', aux,Vtop
filename = 'hist_Vpp' //trim(adjustl(aux))// '.d'
write(*,*) 'filename=',filename
  
allocate( hist(nbin) )
dx = 10._np/real(nbin,np)
hist(:) = 0
open(21,file='Vpp_1.d',status='replace')
do l = 1,pnum
    write(21,*) l, Vpp(l)
    ibin = int(Vpp(l)) + 1
    hist(ibin) = hist(ibin) + 1
enddo
close(21)
open(22,file=filename,status='replace')
do i = 1,nbin
    write(22,*) (real(i-1,np) + 0.5_np)*dx ,hist(i)
enddo
close(22)

write(*,*) 'Ya asignó un valor de voltaje a cada partícula'

7 FORMAT (A3)
8 FORMAT (5(3X,A15))
9 FORMAT (5(2X,ES17.9)) 
                
10 FORMAT (3(3X,A15))
11 FORMAT (3(2X,ES17.9)) 
                
end program V_test
