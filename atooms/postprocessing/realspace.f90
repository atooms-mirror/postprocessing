module realspace_module

   implicit none

contains

  subroutine self_overlap(pos1, pos2, a_square, q)
    real(8),intent(in) :: pos1(:,:), pos2(:,:), a_square
    real(8),intent(out) :: q
    real(8) :: tmp(size(pos1,1))
    integer :: i, iq
    do i = 1,size(pos1,1)
       tmp(i) = sum((pos1(i,:) - pos2(i,:))**2)
    end do
    q = count(tmp .lt. a_square)
  end subroutine self_overlap

end module realspace_module
