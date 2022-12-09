module MDL
    IMPLICIT NONE


CONTAINS

    SUBROUTINE MDLpar(npts, psub, MDLpar_length)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! DESCRIPTION: Find the len between p(startIndex) and p(currIndex) if the
    ! two points are both characteristic points. Minimum descriptive length
    ! distance.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        use distance_functions, only: LOG2, LENGTH, dperp_angle

        implicit none

        INTEGER, INTENT(IN) :: npts
        REAL, DIMENSION(3,npts), INTENT(IN) :: psub !x,y,z of points
        REAL, INTENT(OUT) :: MDLpar_length

        !internal variables
        INTEGER inBetween, idx
        REAL, DIMENSION(4) :: totald !total perpendicular angle distances
            ! between each pair of points within partition.  Note x axis
            ! is parallel to one of the vectors.
        REAL l_h, l_d_h

        !number of points in between start and end of partition
        inBetween = npts-1
        !l_h = Euclidean distance between the two points
        l_h = LOG2(LENGTH(psub(:,1), psub(:,npts)))

        if ( inBetween .eq. 0) then
            l_d_h = 0
        else
            !find perpendicular and angle distance between all in between pts
            idx = 2
            l_d_h = 0.
            totald(:) = 0.
            do while ( idx .le. npts )
                CALL dperp_angle(psub(:, 1),psub(:,npts),psub(:, idx-1),&
                    psub(:, idx),totald)
                idx = idx + 1
            end do
        end if
        !sum everything
        l_d_h = LOG2(totald(1)+totald(2)) + LOG2(totald(3)+totald(4))
        MDLpar_length = l_h + l_d_h

    END SUBROUTINE MDLpar

    SUBROUTINE MDLnopar(npts, psub, MDLnopar_length)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! DESCRIPTION:  Find the distance between p(startIndex) and p(currIndex)
    ! if none of the points in between are characteristic points.  Basically
    ! the length of the trajectory if you trace right along the trajectory.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        use distance_functions, only: LOG2, LENGTH

        implicit none

        integer, intent(in) :: npts
        REAL, DIMENSION(3,npts), INTENT(IN) :: psub !x,y,z of points
        REAL, INTENT(OUT) :: MDLnopar_length

        !internal variables
        integer :: idx
        real :: total_length, l_h, l_d_h

        total_length = 0
        do idx = 1, npts-1
            total_length = total_length + LENGTH(psub(:,idx),psub(:,idx+1))
        end do
        l_h = LOG2(total_length)
        l_d_h = 0
        MDLnopar_length = l_h + l_d_h

    END SUBROUTINE MDLnopar

end module
