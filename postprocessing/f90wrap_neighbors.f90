! Module neighbors_module defined in file neighbors.f90

subroutine f90wrap_neighbors(box, pos, ids, rcut, nn, neigh, n0, n1, n2, n3, n4, &
    n5, n6, n7, n8)
    use neighbors_module, only: neighbors
    implicit none
    
    real(8), intent(in), dimension(n0) :: box
    real(8), intent(in), dimension(n1,n2) :: pos
    integer, intent(in), dimension(n3) :: ids
    real(8), intent(in), dimension(n4,n5) :: rcut
    integer, intent(inout), dimension(n6) :: nn
    integer, intent(inout), dimension(n7,n8) :: neigh
    integer :: n0
    !f2py intent(hide), depend(box) :: n0 = shape(box,0)
    integer :: n1
    !f2py intent(hide), depend(pos) :: n1 = shape(pos,0)
    integer :: n2
    !f2py intent(hide), depend(pos) :: n2 = shape(pos,1)
    integer :: n3
    !f2py intent(hide), depend(ids) :: n3 = shape(ids,0)
    integer :: n4
    !f2py intent(hide), depend(rcut) :: n4 = shape(rcut,0)
    integer :: n5
    !f2py intent(hide), depend(rcut) :: n5 = shape(rcut,1)
    integer :: n6
    !f2py intent(hide), depend(nn) :: n6 = shape(nn,0)
    integer :: n7
    !f2py intent(hide), depend(neigh) :: n7 = shape(neigh,0)
    integer :: n8
    !f2py intent(hide), depend(neigh) :: n8 = shape(neigh,1)
    call neighbors(box=box, pos=pos, ids=ids, rcut=rcut, nn=nn, neigh=neigh)
end subroutine f90wrap_neighbors

! End of module neighbors_module defined in file neighbors.f90

