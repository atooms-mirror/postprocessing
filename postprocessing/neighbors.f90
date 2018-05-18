! module neighbors_module

!   implicit none

! contains

  !! Find nearest neighbors 
  !!
  !! Neighbors list is in C-order for easier handling in python.
  !! Also we give up III law because this avoids sorting the neighbors
  !! indices later and overall this is more efficient
  !! 
  subroutine neighbors(offset,box,pos,ids,rcut,nn,neigh)
    character(1), intent(in)  :: offset
    real(8), intent(in)  :: box(:)
    real(8), intent(in)  :: pos(:,:)
    real(8), intent(in)  :: rcut(:,:)
    integer, intent(in)  :: ids(:)
    integer, intent(inout) :: nn(:)
    integer, intent(inout) :: neigh(:,:)
    real(8)              :: rij(size(pos,1)), rijsq, hbox(size(pos,1))
    real(8)              :: rcutsq(size(rcut,1),size(rcut,2))
    integer              :: i, j, isp, jsp, delta
    if (offset == 'C') then
       delta = -1
    else
       delta = 0
    end if    
    nn = 0
    hbox = box / 2
    rcutsq = rcut**2
    do i = 1,size(pos,2)
       isp = ids(i)
       do j = i+1,size(pos,2)
          jsp = ids(j)
          rij = pos(:,i) - pos(:,j)
          where (abs(rij) > hbox)
             rij = rij - sign(box,rij)
          end where
          rijsq = dot_product(rij,rij)
          if (rijsq < rcutsq(isp,jsp)) then
             nn(i) = nn(i) + 1
             nn(j) = nn(j) + 1
             neigh(i,nn(i)) = j+delta
             neigh(j,nn(j)) = i+delta
          end if
       end do
    end do
  end subroutine neighbors

!end module neighbors_module
