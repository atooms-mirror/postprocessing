module neighbors_module

  implicit none

contains

  !! Find nearest neighbors 
  !!
  !! Neighbors list is in C-order for easier handling in python.
  !! Also we give up III law because this avoids sorting the neighbors
  !! indices later and overall this is more efficient
  !! 
  subroutine neighbors(box,pos,ids,rcut,nn,neigh)
    real(8), intent(in)  :: box(:)
    real(8), intent(in)  :: pos(:,:)
    real(8), intent(in)  :: rcut(:,:)
    integer, intent(in)  :: ids(:)
    integer, intent(out) :: nn(:)
    integer, intent(out) :: neigh(:,:)
    real(8)              :: rij(size(pos,1)), rijsq, hbox(size(pos,1))
    real(8)              :: rcutsq(size(rcut,1),size(rcut,2))
    integer              :: i, j, isp, jsp
    nn = 0
    hbox = box / 2
    rcutsq = rcut**2
    do i = 1,size(pos,2)
       isp = ids(i)
       do j = 1,size(pos,2)
          if (i==j) cycle
          jsp = ids(j)
          rij = pos(:,i) - pos(:,j)
          where (abs(rij) > hbox)
             rij = rij - sign(box,rij)
          end where
          rijsq = dot_product(rij,rij)
          if (rijsq < rcutsq(isp,jsp)) then
             nn(i) = nn(i) + 1
             neigh(i,nn(i)) = j
          end if
       end do
    end do
  end subroutine neighbors

end module neighbors_module
