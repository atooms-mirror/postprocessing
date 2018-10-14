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

!   SUBROUTINE setup_expo_sphere(k0,kmax,positions,expo)
!     REAL(8), INTENT(in) :: k0
!     INTEGER,  INTENT(in) :: kmax
!     REAL(8), INTENT(in) :: positions(:,:)
!     COMPLEX(8),POINTER  :: expo(:,:,:)
!     INTEGER              :: i,j,ndims,npart
!     npart = SIZE(positions,1)
!     ndims = SIZE(positions,2)
!     IF (.NOT.ASSOCIATED(expo)) ALLOCATE(expo(npart,ndims,-kmax:kmax))
!     expo(:,:,0)  = (1.0_8,0.0_8)
!     DO j = 1,ndims
!        expo(:,j,1)  = EXP(IM*k0*positions(:,j))
!        expo(:,j,-1) = CONJG(expo(:,j,1))
!     END DO
!     DO i=2,kmax
!        expo(:,:,i)  = expo(:,:,i-1) * expo(:,:,1)
!        expo(:,:,-i) = CONJG(expo(:,:,i))
!     END DO
!   END SUBROUTINE setup_expo_sphere

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
       rho(ii) = sum(expo(:,1,i1) * expo(:,2,i2) * expo(:,3,i3))
       !sk(kbin) = sk(kbin) + rho*conjg(rho)
       !sk(kbin) = sk(kbin) + real(rho*conjg(rho),8)
    end do
  end subroutine sk_bare


  ! subroutine fskt(position,k0,kmax,ikvec,ikbin,fkt)
  !   real(8),intent(in)        :: position(:,:,:)  ! (ndim, npart, nsteps)
  !   real(8),intent(in)        :: k0
  !   integer, intent(in)       :: ikvec(:,:), ikbin(:)  ! (ndim, nvec), (nvec)
  !   real(8), intent(inout)    :: fkt(:,:)
  !   integer                   :: nmax(3) 
  !   complex(8),     pointer   :: expo(:,:)
  !   integer     :: i1, i2, i3, ii, npart, ndim
  !   integer     :: kbin, kmax
  !   complex(8) :: rho
  !   ndim = size(position,1)
  !   npart = size(position,2)
  !   nstep = size(position,3)

  !   ! Determine n. of time origins for each time separation
  !   nt = 0
  !   DO it = 1,SIZE(t_delta)
  !      t = t_delta(it)
  !      off = t_offset(it)
  !      DO s = off, nstep-t, tskip
  !         nt = nt + 1
  !         tt = t_bin(nt)  ! indirect reference needed here
  !         nori(tt) = nori(tt) + 1
  !      END DO
  !   END DO

  !   fkt = 0.0_8
  !   DO i = 1,npart
  !      !CALL setup_expo(kmin,nmax,kmax,TRANSPOSE(position(:,i,:)),expo)
  !      CALL setup_expo(k0,kmax,position(:,i,:),expo)
  !      DO ik = 1,SIZE(ik_vec)
  !         nt = 0
  !         DO it = 1,SIZE(t_delta)
  !            t = t_delta(it)
  !            off = t_offset(it)
  !            DO s = off, nstep-t, tskip
  !               nt = nt + 1
  !               tt = t_bin(nt)  ! indirect reference needed here
  !               cor(tt) = cor(tt) + &
  !                    expo(s+t,1,ik_vec(1,ik)) * CONJG(expo(s,1,ik_vec(1,ik))) * &
  !                    expo(s+t,2,ik_vec(2,ik)) * CONJG(expo(s,2,ik_vec(2,ik))) * &
  !                    expo(s+t,3,ik_vec(3,ik)) * CONJG(expo(s,3,ik_vec(3,ik))) 
  !            END DO
  !         END DO
  !         fkt(ik_norm(ik),:) = fkt(ik_norm(ik),:) + REAL(cor(:) / nori(:),8)
  !      END DO
  !   END DO    
  ! end subroutine fskt

end module fourierspace_module
