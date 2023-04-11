program dendritas_03
use precision
use mzranmod
implicit none

!Definicion de variables
real(np)                                ::ti,tf,t,dt !Variables temporales
real(np)                                ::D,q,Long !coeficientes
real(np)                                ::rli0,dx,dy,Lx,Ly, x,y,datt, distx,disty,dist
integer,dimension(:),allocatable       ::ex,ey,e0x,e0y !vectores espaciales de los iones
integer                                ::exs,eys !escalares
integer,dimension(:),allocatable       ::li_xd,li_yd !Arreglo del litio depositado sobre el anodo
integer,dimension(:),allocatable       ::li0x,li0y !li0
integer,dimension(:),allocatable       ::li_aux_x, li_aux_y
integer                                 ::gx,gy !vector unitario aleatorio.
integer                                 ::n0,n0max,nm !numeros de particulas 
integer                                 ::nn !numero de voxels
integer                                 ::nt !numeros de pasos temporales 
integer                                 ::i,j,k,m,l
integer                                 ::nq,nrx,nry
integer                                 ::ndistx,ndisty,ndist,ndatt
real(np),dimension(:,:),allocatable     ::v,v0
real(np)                                ::dv,vm
real(np),dimension(:,:),allocatable     ::E0_x,E0_y,E_x,E_y
real(np),dimension(:),allocatable       ::rx,ry
real(np)                                ::u

!Definicion de variables ******************************************

dt = 1.0e-6_np
n0 = 77                     !comienzo con este numero de li0 depositados sobre la sup
n0max = 600                 !numero maximo de particulas li0
nm = 50                     !numero de iones
nn = 154                    !numero de puntos en la malla
nt = 100                    !numero de pasos temporales
rli0 = 1.67e-10_np          !radio atomico del li0
D = 1.4e-14_np              !coef de dif del Li+ en el electrolito
u = 5.6e-13_np              !movilidad del li+
Long = 16.7e-9_np           !Longitud de la caja simulada
datt = 1.3_np*rli0/Long     !distancia tal que li+->li0 = rli0/(pi/4)
q = sqrt(2._np*D*dt)/Long   !desplazamiento medio debido a la difusion
vm = 85e-3_np               !magnitud del voltaje
dv = vm/real(nn,np)         !diferencial de voltaje

dx = 0.65_np*rli0/Long !1.3*rli0/2
dy = 0.65_np*rli0/Long

!Conversion de estas variables a enteros

nq=anint(q*154._np) !cantidad de voxels que mueve q
ndatt = 2           !cantidad de voxels que determinan la distancia minima tl que li-->li0

!Inicializacion de los valores
gx = 0
gy = 0
x=0._np
y=0._np
ndistx = 0
ndisty = 0
ndist = 0
m = n0 !numero inicial de li0


allocate(ex(1:nm),ey(1:nm),e0x(1:nm),e0y(1:nm)) !son las variales espaciales que guardan informacion de la posciion de todas las particulas a un dado tiempo t !OJO....si se hace esto el numero nm se debe mantener constante en cada paso temporal. Entonces se debe mantener constante el numero de iones
allocate(li_xd(1:n0),li_yd(1:n0))
allocate(li_aux_x(1:600),li_aux_y(1:600))

!creacion de los Li0 depositados uniformemente sobre el anodo--------------------------------------------------------------
open(21,file='init_li0_03.dat',status='replace')
do i = 1,n0
    !vector que guarda la posicion del litio ordenado a lo largo de x normalizado por Long
    li_xd(i) = 2*real(i,np)
    li_yd(i) = 0
    write(21,*)li_xd(i),li_yd(i)
enddo

!Creacion del electrolito con iones*******************************

!Creacion de la posicion inicial de los 50 iones-------------------
open(22,file='init_ion_03.dat',status='replace')
do i = 1,nm
!Posiciones iniciales normalizadas por la longitud Long
    e0x(i) = anint(rmzran()*154._np)
    e0y(i) = anint(rmzran()*154._np)
    write(22,*)e0x(i),e0y(i)
enddo

!Variables del potencial
allocate(v(0:nn,0:nn))
allocate(v0(0:nn,0:nn))
v0(:,:) = 0._np
v(:,:) = 0._np

!Doy el valor de la matriz del voltaje en t=0. Es una disrubucion uniforme de Vm a -Vm

do j = 1,nn
    do i = 1,nn
        v0(i,j) = vm - real(j,np)*dv
    enddo
enddo

v(:,:) = v0(:,:)

!Variables del campo
allocate(E0_x(0:nn,0:nn),E0_y(0:nn,0:nn),E_x(0:nn,0:nn),E_y(0:nn,0:nn))
allocate(rx(1:nm),ry(1:nm))
!Calculo el campo correspondiente a la distribucion de potencial
call calculate_E(nn,dx,dy,v0,v,vm,E_x,E_y)

open(21,file='rx_ry.dat')

do i = 1,nm
    !Guardo la posicion del li+ para calcular el campo en esa coordenada
    k = e0x(i)
    l = e0y(i)
    rx(i) = u*E_x(k,l)*dt/Long
    ry(i) = u*E_y(k,l)*dt/Long
    write(21,*)i,rx(i),ry(i)
enddo
close(21)
li_aux_x(:) = 0
li_aux_y(:) = 0
ex(:) = 0
ey(:) = 0
exs = 0
eys = 0
nrx = 0
nry = 0


open(21,file='rx_ry_evol.dat')

open(23,file='evol_li0_03.dat',status='replace')
open(24,file='evol_lim_03.dat',status='replace')
!Definicion de las ecuaciones de movimiento browniano
!Evoluion temporal
do j = 1,nt
    do i = 1,nm
        !Definicion del vector unitario aleatorio g
        gx = anint(rmzran()*3._np)-1
        gy = anint(rmzran()*3._np)-1
        !Calculo los enteros de voxels que mueve el campo electrico
        !Vuelvo a calcular los coeficientes con el nuevo campo
        !Acordate que entro con campo E0_x y E0_y
        nrx = anint(rx(i)*154._np)
        nry = anint(ry(i)*154._np)
        write(*,*)nrx
        write(*,*)nry
        !Evolucion de la posicion de los iones
        ex(i) = e0x(i) + nq*gx + nrx
        ey(i) = e0y(i) + nq*gy + nry
        !Meto las PBC en una caja de tamaño nnxnn
        ex(i) = ex(i) - nn*anint(real(ex(i),np)/nn-0.5_np)
        ey(i) = ey(i) - nn*anint(real(ey(i),np)/nn+0.5_np) + nn
        write(*,*)ex(i),ey(i)
        !Definicion de la condicion Li+-->Li0
        do k = 1,m
            if ( m.le.n0) then
                ndistx = ex(i)-li_xd(k) 
                ndisty = ey(i)-li_yd(k)
                ndist = anint(sqrt( real(ndistx*ndistx + ndisty*ndisty ,np)))
            elseif ( m.gt.n0) then
                ndistx = ex(i)-li0x(k) 
                ndisty = ey(i)-li0y(k)
                ndist = anint(sqrt( real(ndistx*ndistx + ndisty*ndisty ,np)))
            endif
            if (ndist.le.ndatt) then 
                m = m+1 !Aumenta la dimension del arreglo del li0
                exs = ex(i)
                eys = ey(i)
                call save_li0(m,li_xd,li_yd,exs,eys,li0x,li0y,li_aux_x,li_aux_y) !Guardo la nueva posicion del li0
                !Tengo que reponer un ion en el espacio de forma aleatoria
                ex(i) = anint(rmzran()*154)
                ey(i) = anint(rmzran()*154)
                !Meto las PBC
                ex(i) = ex(i) - nn*anint(real(ex(i),np)/nn-0.5_np)
                ey(i) = ey(i) - nn*anint(real(ey(i),np)/nn+0.5_np) + nn
                !Se llama a la subrutina que calcula la nueva distrbucion de potencial debido al nuevo li0
                call calculate_v(nn,exs,eys,dx,dy,v0,v,vm)
                !Luego sepuede llamar al  !Renombro la variable para tener campo en t=0campo
                !write(*,*) v(exs,eys),'paso li+ a li0'
                call calculate_E(nn,dx,dy,v0,v,vm,E_x,E_y)
                write(*,*)E_x(exs,eys)
            else
                !Necesito poner un else porque de todas formas tengo que allocatear el li0x y li0y mas alla de si se le pega un ion o no !se supone que la linea siguiente llama a la subrutina pero conserva el numero m del litio depositado entonces no deberia cambiar
                call save_li0(m,li_xd,li_yd,exs,eys,li0x,li0y,li_aux_x,li_aux_y)
                !No hay que llamar a calcular v porque el potencial va a ser el mismo del paso anterior (si no paso que li+--->li0), lo mismo con el campo
            endif
        enddo
        !Calculo el vector rx y ry de nuevo con el nuevo campo,hay que poner un if para no calcularlo si no cambio
        do l = 1, nm
            rx(l) = u*E_x(ex(l),ey(l))*dt/Long
            ry(l) = u*E_y(ex(l),ey(l))*dt/Long
            write(21,*)l,rx(l),ry(l),E_x(ex(l),ey(l)),E_y(ex(l),ey(l))
        enddo
        v0(:,:) = v(:,:) 
    enddo
    !Renombre de posiciones para t+dt
    t = real(j,np)*dt
    e0x = ex
    e0y = ey
    write(*,*) real(j,np)/real(nt,np)*100._np
enddo

close(21)

write(*,*)m
write(*,*)nm

do l = 1,nm
    write(24,*)ex(l),ey(l)
enddo

do l = 1,m
    write(23,*)li0x(l),li0y(l)
enddo
    
deallocate(ex,ey)
close(21)
close(22)
close(23)
close(24)

!*******************************************************************
contains

Subroutine save_li0(mm,li_xxd,li_yyd,exx,eyy,li0xx,li0yy,li_auxx,li_auxy)
integer,intent(in)                             ::mm
integer,intent(in)                             ::exx,eyy
integer,dimension(1:77),intent(in)             ::li_xxd, li_yyd
integer,dimension(:),allocatable,intent(out)   ::li0xx, li0yy
integer,dimension(1:600),intent(inout)         ::li_auxx, li_auxy
integer                                        ::i,nn
allocate(li0xx(1:mm),li0yy(1:mm))

!En los primeros 100 lugares guardo el litio depositado sobre el anodo
do i = 1,77
    li_auxx(i) = li_xxd(i)
    li_auxy(i) = li_yyd(i)
enddo
!En los siguientes lugares tengo que guardar las nuevas posiciones de los li+ que pasan a ser li0 

!me creo un vector auxiliar li_auxx,li_auxy. El vector auxiliar tiene en total 600 entradas, las cuales son nulas y las va llenando a medida que m>100 y se va guardando justamente las posiciones de cada li+ que paso a li0. De esta forma tengo una especie de lista donde guardo estos valores. y en li0xx y li0yy me guardo las coordenadas no nulas :D

!Se llena la entrada m
li_auxx(mm) = exx 
li_auxy(mm) = eyy 

do i = 1,mm
    li0xx(i) = li_auxx(i)
    li0yy(i) = li_auxy(i)
enddo

!Doy las PBC
nn=154
do i = 1,mm
    li0xx(i) = li0xx(i) - 154*anint(real(li0xx(i),np)/154._np-0.5_np)
    li0yy(i) = li0yy(i) - 154*anint(real(li0yy(i),np)/154._np+0.5_np)+154
enddo
    
end Subroutine save_li0

Subroutine calculate_v(nnn,exss,eyss,dxx,dyy,v00,vv,vmm)
integer,intent(in)                      ::nnn,exss,eyss
real(np), intent(in)                    ::vmm,dxx,dyy
real(np),dimension(0:nnn,0:nnn),intent(inout) ::v00,vv
integer                                 ::i,j
real(np)                                ::f
!VV se tiene que calcular a partir de v00. Tienen que ser entradas enteras las posiciones de los nuevos li0 adquiridos para renovar los valores -Vm en la matriz para dicha entrada

f = 2._np/(dxx*dxx)+2._np/(dyy*dyy)


!vv=v00 solo modifico la entrada (exss,eyss) donde se encuentra el li0
v00(exss,eyss)=vmm
!Ahora con esta nueva modificacion, se calcula el potencial resultante en la malla

!Primero lleno los valores de los vertices de la malla
vv(1,1) = vmm
vv(nnn,1) = vmm
vv(1,nnn) = 0._np
vv(nnn,nnn) = 0._np
!Luego separo los bordes en x con su ciclicidad por las PBC
do j=2,nnn-1
    vv(1,j) = ((v00(2,j)+v00(nnn,j))/(dxx*dxx)+(v00(1,j+1)+v00(1,j-1))/(dy*dy))*f
    vv(nnn,j) = ((v00(1,j)+v00(nnn-1,j))/(dxx*dxx)+(v00(1,j+1)+v00(1,j-1))/(dy*dy))*f
enddo
!Ahora completo el resto de valores en la malla
do j=2,nnn-1
    do i =2,nnn-1
        vv(i,j) = ((v00(i+1,j)+v00(i-1,j))/(dxx*dxx)+(v00(i,j+1)+v00(i,j-1))/(dyy*dyy))*f
    enddo
enddo

!Tengo que volver a imponer el potencial v en el sitio que se encuentra el li0

vv(exss,eyss)=vm

end Subroutine calculate_v

Subroutine calculate_E(nnn,dxx,dyy,v00,vv,vmm,E_xx,E_yy)
integer,intent(in)                      ::nnn
real(np), intent(in)                    ::dxx,dyy,vmm
real(np),dimension(0:nnn,0:nnn),intent(in) ::v00,vv
real(np),dimension(0:nnn,0:nnn),intent(inout) ::E_xx,E_yy
integer                                 ::i,j

!Escribo el interior
do j = 1,nnn
    do i = 2,nnn-1
        E_xx(i,j) = -(vv(i+1,j)-vv(i-1,j))/dxx
    enddo
enddo
!Escribo los bordes considerando las PBC
do j = 1,nnn
    E_xx(1,j) = -(vv(2,j)-vv(nnn,j))/dxx
    E_xx(nnn,j) = -(vv(1,j)-vv(nnn-1,j))/dxx
enddo


do j = 2,nnn-1
    do i = 1,nnn
        E_yy(i,1) = -(vv(i,2)-vmm)/dyy
        E_yy(i,j) = -(vv(i,j+1)-vv(i,j-1))/dyy
        E_yy(i,nnn) = vv(i,nnn-1)/dyy
    enddo
enddo

end Subroutine calculate_E

end program dendritas_03
