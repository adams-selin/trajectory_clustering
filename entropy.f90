SUBROUTINE entropy(n, epsilon, h, distance, num_segments, opt_MinLns)
!  Entropy formula as defined by Lee et al. Eq. 10.
    use cluster_module, only: calc_Neps
    use distance_functions, only: LOG2

    IMPLICIT NONE
    INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(14, 60)

    INTEGER, INTENT(IN)    :: n
    REAL (dp), INTENT(IN)  :: epsilon(:)
    INTEGER(8), INTENT(IN) :: num_segments
    REAL, INTENT(IN)       :: distance(:)
    REAL (dp), INTENT(OUT) :: h
    REAL, INTENT(OUT)      :: opt_MinLns !avg |Ne(L)|, see text discussion in
                                         !Lee et al. after eq. 10

    ! Local variables
    INTEGER(8)   :: i
    !array containg the number of line segments within a neighborhood of
    !epsilon of this line segment
    INTEGER(8), DIMENSION(num_segments) :: num_neigh_eps
    !total number of neighborhood line segments (sum of num_neigh_eps array)
    INTEGER(8) :: sum_neigh_eps
    !Variables needed to pass to calc_Neps subroutine
    integer(8) :: num_distances
    integer(8), dimension(100000) :: Neps_segID
    real, dimension(100000) :: cluster_distance
    integer :: num_segIDs !number of segments within epsilon neighborhood
    real :: p_x

    h = 0.0_dp
    num_distances = num_segments * (num_segments-1) / 2
    num_neigh_eps(:) = 0

    if ( n .ne. 1 ) then
        stop "Trying to minimize too many variables; n must be set to 1"
    end if
    

    !$OMP PARALLEL DO
    do i = 1, num_segments
        CALL calc_Neps(i, epsilon(n), distance, num_distances, num_segments, &
                num_segIDs, Neps_segID, cluster_distance)
        num_neigh_eps(i) = num_segIDs
    end do
    !$OMP END PARALLEL DO

    sum_neigh_eps = SUM(num_neigh_eps)

    do i = 1, num_segments
        p_x = REAL(num_neigh_eps(i))/REAL(sum_neigh_eps)
        h = h - DBLE(p_x) * DBLE(LOG2(p_x))
    end do

    opt_MinLns = REAL(sum_neigh_eps) / REAL(num_segments)
    !print *, '      opt_MinLns: ', opt_MinLns


    RETURN
END SUBROUTINE entropy
