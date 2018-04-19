module fourierspace_module

   implicit none

contains

  subroutine setup_expo(k0,nmax,positions,expo)
    real(8), intent(in) :: k0
    integer, intent(in) :: nmax(:) !,kmax
    real(8),intent(in) :: positions(:,:)
    complex(8),intent(out) :: expo(:,:,:)
    integer :: i,j,ndims,npart
    complex(8), parameter :: IM = (0.d0, 1.d0)
    return
    npart = size(positions,1)
    ndims = size(positions,2)
!    if (associated(expo) .and. size(expo,1)/=size(positions,1)) then
!       deallocate(expo)
!    end if
!    if (.not.associated(expo)) then
!       allocate(expo(npart,ndims,-kmax:kmax))
!    end if
    expo(:,:,0)  = (1.D0,0.D0)
    do j = 1,ndims
       expo(:,j,1) = exp(IM*k0*positions(:,j))
       do i=2,nmax(j)
          expo(:,j,i)=expo(:,j,i-1)*expo(:,j,1)
       end do
    end do
    do i=-nmax(2),-1
       expo(:,2,i)=conjg(expo(:,2,-i))
    end do
    do i=-nmax(3),-1
       expo(:,3,i)=conjg(expo(:,3,-i))
    end do
  end subroutine setup_expo

  subroutine sk_one(position,k0,kmax,ikvec,ikbin,sk)
    real(8),intent(in)        :: position(:,:)  ! (ndim, npart)
    real(8),intent(in)        :: k0
    integer, intent(in)       :: ikvec(:,:), ikbin(:)  ! (ndim, nvec), (nvec)
!    complex(8), intent(inout) :: sk(:)
    real(8), intent(inout) :: sk(:)
    integer                   :: nmax(3)  !?
    complex(8),     pointer   :: expo(:,:,:)
    integer     :: i1, i2, i3, ii, npart, ndim
    integer     :: kbin, kmax
    complex(8) :: rho
    ndim = size(position,1)
    npart = size(position,2)
    allocate(expo(npart,ndim,-kmax:kmax))
    call setup_expo(k0,nmax,position,expo)
    sk = (0.d0, 0.d0)
    print*, '----', ndim, npart, size(ikbin)
    do ii = 1,size(ikvec,2)
       i1   = ikvec(1,ii)
       i2   = ikvec(2,ii)
       i3   = ikvec(3,ii)
       kbin = ikbin(ii)
       rho = sum(expo(:,1,i1) * expo(:,2,i2) * expo(:,3,i3))
       !sk(kbin) = sk(kbin) + rho*conjg(rho)
       sk(kbin) = sk(kbin) + real(rho*conjg(rho),8)
    end do
  end subroutine sk_one

  subroutine sk_bare(expo,ikvec,rho)
    complex(8),intent(in)        :: expo(:,:,:)  ! (ndim, npart)
!    real(8),intent(in)        :: k0
!    integer, intent(in)       :: ikbin(:),ikvec(:,:) ! (ndim, nvec), (nvec)
    integer, intent(in)       :: ikvec(:,:) ! (ndim, nvec), (nvec)
!    complex(8), intent(inout) :: sk(:)
!    real(8), intent(inout) :: sk(:)
    integer                   :: nmax(3)  !?
!    complex(8),     pointer   :: expo(:,:,:)
    integer     :: i1, i2, i3, ii, npart, ndim
    integer     :: kbin, kmax
    complex(8), intent(inout) :: rho(:)
    ! ndim = size(position,1)
    ! npart = size(position,2)
    ! allocate(expo(npart,ndim,-kmax:kmax))
    ! call setup_expo(k0,nmax,position,expo)
    ! sk = (0.d0, 0.d0)
    do ii = 1,size(ikvec,2)
       i1   = ikvec(1,ii)
       i2   = ikvec(2,ii)
       i3   = ikvec(3,ii)
       !kbin = ikbin(ii)
       !if (ii==4) print*, expo(:,1,i1) * expo(:,2,i2) * expo(:,3,i3)
       rho(ii) = sum(expo(:,1,i1) * expo(:,2,i2) * expo(:,3,i3))
       !sk(kbin) = sk(kbin) + rho*conjg(rho)
       !sk(kbin) = sk(kbin) + real(rho*conjg(rho),8)
    end do
  end subroutine sk_bare

end module fourierspace_module
