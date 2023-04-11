program pos_inicial
  use gems_groups, only: group, atom, atom_dclist, sys, gindex
  use gems_neighbor, only: test_update, ngroup, ngindex, nn_vlist, nupd_vlist
  use gems_elements, only: set_z, elements_init,elements
  use gems_errors, only: timer_dump, timer_start, sclock_t1, sclock_t2, wstd, sclock_max, sclock_rate,logunit
  use gems_constants, only: time1, dp, sp

  implicit none



subroutine pos_inic()
! Create initial positions. For 1 M [Li+]
use gems_program_types, only: distance 
use gems_program_types, only: tbox, box_setvars  
integer             :: n 
real(dp),allocatable:: r(:,:) !posic.
real(dp)            :: Mol, v1(3),v2(3),dif_pos(3), alto, dif2 ! molaridad, diferencia entre vectores posic., alto caja sim.
real(dp),parameter  :: r0=3.2_dp, mLi= 6.94_dp
integer             :: i,j,l,k
logical,parameter   :: pbc(3)=[.true.,.true.,.false.]

! Set the box
tbox(:,:)= 0._dp
tbox(1,1)= xi
tbox(2,2)= yi
tbox(3,3)= zmax
call box_setvars()

! Calculo nro. de partículas
! n= Molaridad * Volumen(A^3) * 1e-27 L/A^3 * 6.022e23 (Nro Avogadro)
! Molaridad usada = 1 M
Mol = 1._dp 
alto = zmax - o
n = Mol * xi * yi * alto * 6.022e-4 
allocate(r(n,3))

open(12,File='pos_inic.xyz')
write(12,*) n  ! nro. total de átomos
write(12,*)

! Recorro todas las particulas a crear
do i=1,n

  ! Recorro los intentos 
  intento: do k=1,10000
    !Posic. de Li
    r(i,1)=ran(idum)*xi
    r(i,2)=ran(idum)*yi
    r(i,3)=(ran(idum)*alto) + o ! en z 

    ! Recorro las particulas ya creadas 
    do j=1,i-1 

      !Veo distancia Li-Li
      v1(:)=r(i,:)
      v2(:)=r(j,:)
      dif_pos(:)=distance(v2,v1,pbc)
      dif2=dot_product(dif_pos,dif_pos)

      if (dif2<r0*r0) cycle intento

    enddo

    exit
  enddo intento

  if(k==10001) then
      print *, 'Maximo numero de intentos alcanzado'
      stop
  endif

write(12,'(a,4(x,f25.12))') 'Li',r(i,1),r(i,2),r(i,3),mLi

enddo
close(12)

end subroutine pos_inic

end program pos_inicial
