module distance_functions

CONTAINS

    subroutine distance_matrix(segments,num_segments,distance,num_distances)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! DESCRIPTION:  create a distance matrix for the segments array containing
    ! the distance between each pair of segments. Use indices_1d distance
    ! to flatten indices to a 1d array to save space.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        implicit none
        integer(8), intent(in) :: num_segments, num_distances
        real, intent(in) :: segments(7, num_segments)
        real, intent(out) :: distance(num_distances)
        !internal variables
        integer(8) :: i, j
        integer(8) :: ind1d
        real, dimension(4) :: perp_angles
        real :: parallel

        !find the upper right half of the matrix, sans diagonal
        !$OMP PARALLEL DO
        do i = 1, num_segments-1
            print *, 'i, total: ', i, num_segments
            do j = i+1, num_segments
                perp_angles(:) = 0.
                CALL dperp_angle(segments(1:3,i),segments(4:6,i),&
                                 segments(1:3,j),segments(4:6,j),perp_angles)
                parallel = 0.
                CALL dparallel(segments(1:3,i),segments(4:6,i),&
                               segments(1:3,j),segments(4:6,j),parallel)
                ! convert our 2d indices to a 1d flattened array
                CALL indices_1d(i,j,num_segments,ind1d)
                !write(*,'(a30,i6,1x,i6,1x,i11,1x,i11)') &
                !    ' i, j, ind1d, size(distance): ', i, j, ind1d, &
                !    size(distance, dim=1)
                distance(ind1d) = parallel + SUM(perp_angles)
                ! if ( i .eq. 836 ) then
                !     print *, 'i, j, parallel, perp_angles, total: '
                !     print *, i, j, parallel, perp_angles, distance(ind1d)
                ! end if
            end do
        end do
        !$OMP END PARALLEL DO

    end subroutine distance_matrix


    subroutine indices_1d(i,j,n,ind1d)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! DESCRIPTION:  convert 2d num_segments * num_segments indices to 1d
    ! flattened array of size (num_segments)*(num_segments - 1) / 2
    ! Note that the sum of all numbers between/including m and n is equivalent to
    !   n*(n+1)/2 - m*(m-1)/2
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        implicit none
        integer(8), intent(in) :: i, j
        integer(8), intent(in) :: n !num_segments
        integer(8), intent(out) :: ind1d
        integer(8) :: dumi, dumj

        !symmetric matrix, so make sure i is the smaller number
        if ( i .gt. j ) then
            dumi = j
            dumj = i
        else
            dumj = j
            dumi = i
        end if

        ind1d = (n-1)*n/2 - (n-dumi+1)*(n-dumi)/2 + dumj-dumi

    end subroutine indices_1d


    subroutine dparallel(ps,pe,qs,qe,distance)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! DESCRIPTION:  calculate “parallel distance” along the x axix, between two
    ! trajectories as designed by Lee et al. for two dimensions and modified by
    ! RAS to 3D.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        implicit none

        !starting, ending points of vectors in x,y,z space
        REAL, DIMENSION(3), INTENT(IN) :: ps, pe, qs, qe
        REAL, INTENT(INOUT) :: distance

        !internal variables
        real :: l1, l2
        !start, ending points of vectors in rotated space
        real, dimension(3) :: psPrime, qsPrime, pePrime, qePrime

        !Transform line q to the vector space of p.  Assume the vector p is
        !defined by the equation p = ax + by + cz.  After the transformation,
        !p should align with the x’ axis in the new x’, y’, z’ space, starting
        !at the origin.
        CALL rotate_vector(ps,pe,qs,qe,psPrime,pePrime,qsPrime,qePrime)


        !parallel distances along new x' axis.
        l1 = ABS(qsPrime(1))  !psPrime(1) is 0
        l2 = ABS(qePrime(1) - pePrime(1))
        distance = MIN(l1, l2)

    end subroutine dparallel



    subroutine dperp_angle(ps,pe,qs,qe,distance)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! DESCRIPTION:  calculate “perpendicular distances” along two axes and
    ! “angle distances”, again along two axes, between two trajectories as
    ! designed by Lee et al. for two dimensions and modified by RAS to 3D.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        implicit none
        !starting, ending points of vectors in x,y,z space
        REAL, DIMENSION(3), INTENT(IN) :: ps, pe, qs, qe
        REAL, DIMENSION(4), INTENT(INOUT) :: distance

        !internal variables
        real :: L1y, L1z, L2y, L2z, dperpy, dperpz
        real :: qa, qb, qc, lenq, thetaz, thetay, lenx1, dthetaz, dthetay
        !start, ending points of vectors in rotated space
        real, dimension(3) :: psPrime, qsPrime, pePrime, qePrime
        real, parameter :: pi = 3.14159

        !Transform line q to the vector space of p.  Assume the vector p is
        !defined by the equation p = ax + by + cz.  After the transformation,
        !p should align with the x’ axis in the new x’, y’, z’ space, starting
        !at the origin.
        CALL rotate_vector(ps,pe,qs,qe,psPrime,pePrime,qsPrime,qePrime)

        !Calculate perpendicular distances along new y, z axes.
        L1y = ABS( qsPrime(2)-psPrime(2) )
        L1z = ABS( qsPrime(3)-psPrime(3) )
        L2y = ABS( qePrime(2)-pePrime(2) )
        L2z = ABS( qePrime(3)-pePrime(3) )
        dperpy = (L1y**2 + L2y**2) / (L1y+L2y)
        dperpz = (L1z**2 + L2z**2) / (L1z+L2z)

        !Now the angle distances.
        qa = (qePrime(1)-qsPrime(1))
        qb = (qePrime(2)-qsPrime(2))
        qc = (qePrime(3)-qsPrime(3))
        lenq = LENGTH(qsPrime,qePrime)
        thetaz = ABS( ASIN(qc / lenq))

        if ( thetaz .gt. pi ) then
            thetaz = 2.*pi-thetaz
        end if
        if ( thetaz .le. 0.5*pi ) then
            dthetaz = lenq*SIN(thetaz)
        else
            dthetaz = lenq
        endif
        lenx1 = SQRT(qa**2. + qb**2)
        thetay = ABS( ACOS( qa/lenx1))
        if ( thetay .gt. pi ) then
            thetay = 2*pi-thetay
        end if
        if ( thetay .le. 0.5*pi ) then
            dthetay = lenx1*SIN(thetay)
        else
            dthetay = lenx1
        endif

        !Return the perpendicular and angle distances as an array
        distance(1) = distance(1) + dperpy
        distance(2) = distance(2) + dperpz
        distance(3) = distance(3) + dthetay
        distance(4) = distance(4) + dthetaz

    end subroutine dperp_angle


    subroutine rotate_vector(ps,pe,qs,qe,psPrime,pePrime,qsPrime,qePrime)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! DESCRIPTION:   Given two vectors, p and q, rotate q to the vector space
    ! of p where p is aligned along the x axis starting at the origin.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            implicit none
            real, dimension(3), intent(in) :: ps, pe, qs, qe
            real, dimension(3), intent(out) :: psPrime, qsPrime, pePrime, qePrime

            !internal variables
            real :: a, b, c, phi, lenm, theta
            real, dimension(3) :: qsabc, qeabc
            real, dimension(3,3) :: rotation_matrix

            !Transform line q to the vector space of p.  Assume the vector p
            !is defined by the equation p = ax + by + cz.  After the
            !transformation, p should align with the x’ axis in the new
            !x’, y’, z’ space.
            a = pe(1) - ps(1)
            b = pe(2) - ps(2)
            c = pe(3) - ps(3)

            !Find a, b, c values for qs, qe using ps as the origin.
            qsabc(1) = qs(1) - ps(1)
            qsabc(2) = qs(2) - ps(2)
            qsabc(3) = qs(3) - ps(3)
            qeabc(1) = qe(1) - ps(1)
            qeabc(2) = qe(2) - ps(2)
            qeabc(3) = qe(3) - ps(3)

            !Solve for phi, theta, psi using the inverse cosine equations with
            ! a, b, c, specifically:
            phi = ACOS(a / (a**2. + b**2.)**0.5)
            if (b .lt. 0) then
                phi = -phi
            endif
            lenm = LENGTH(pe,ps)
            theta = -ASIN( c / lenm)

            !Declare the rotation matrix. Transpose so correct in row-order space
            !This rotation matrix rotates the vector -phi about the z axis, and
            !then -theta about the y axis.
            rotation_matrix = transpose(reshape((/&
                cos(theta)*cos(phi), cos(theta)*sin(phi), -sin(theta), &
                -sin(phi),           cos(phi),            0.,           &
                sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)/), &
                shape(rotation_matrix)))

            !Solve for the new p, q start and endpoints in the prime matrix
            psPrime = (/0.,0.,0./)
            pePrime = (/lenm,0.,0./)
            !matrix multiplication to get the answer
            qsPrime = MATMUL(rotation_matrix,qsabc)
            qePrime = MATMUL(rotation_matrix,qeabc)

    end subroutine rotate_vector


    subroutine inverse_rotate(ps,pe,qs,qe,psPrime,pePrime,qsPrime,qePrime)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! DESCRIPTION:   Given two vectors, p and q, where q is in the vector space
    ! of p where p is aligned along the x axis starting at the origin, return
    ! q to cartesian coordinates.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            implicit none
            real, dimension(3), intent(in) :: ps, pe, qs, qe
            real, dimension(3), intent(out) :: psPrime, qsPrime, pePrime, qePrime

            !internal variables
            real :: a, b, c, phi, lenm, theta
            real, dimension(3) :: qsabc, qeabc
            real, dimension(3,3) :: inverse_matrix

            !Transform line q to the vector space of p.  Assume the vector p
            !is defined by the equation p = ax + by + cz.  After the
            !transformation, p should align with the x’ axis in the new
            !x’, y’, z’ space.
            a = pe(1) - ps(1)
            b = pe(2) - ps(2)
            c = pe(3) - ps(3)

            !Find a, b, c values for qs, qe using ps as the origin.
            qsabc(1) = qs(1) - ps(1)
            qsabc(2) = qs(2) - ps(2)
            qsabc(3) = qs(3) - ps(3)
            qeabc(1) = qe(1) - ps(1)
            qeabc(2) = qe(2) - ps(2)
            qeabc(3) = qe(3) - ps(3)

            !Solve for phi, theta, psi using the inverse cosine equations with
            ! a, b, c, specifically:
            phi = ACOS(a / (a**2. + b**2.)**0.5)
            if (b .lt. 0) then
                phi = -phi
            endif
            lenm = LENGTH(pe,ps)
            theta = -ASIN( c / lenm)

            !Declare the rotation matrix. Transpose so correct in row-order space
            !This rotation matrix rotates the vector theta about the y axis, and
            !then phi about the z axis.
            inverse_matrix = transpose(reshape((/&
                cos(theta)*cos(phi),   -sin(phi),  sin(theta)*cos(phi), &
                cos(theta)*sin(phi),    cos(phi),  sin(theta)*sin(phi), &
                -sin(theta),                  0.,  cos(theta)/), &
                shape(inverse_matrix)))

            !Solve for the new p, q start and endpoints in the prime matrix
            psPrime = (/0.,0.,0./)
            pePrime = (/lenm,0.,0./)
            !matrix multiplication to get the answer
            qsPrime = MATMUL(inverse_matrix,qsabc)
            qePrime = MATMUL(inverse_matrix,qeabc)

    end subroutine inverse_rotate


    real function LENGTH(ps,pe)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! DESCRIPTION:  Euclidean distance between two points.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        implicit none
        real, DIMENSION(3):: ps, pe

        length =  SQRT( (pe(1) - ps(1))**2 + (pe(2) - ps(2))**2 + &
                        (pe(3) - ps(3))**2 )
    end function LENGTH

    real function LOG2(x)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! DESCRIPTION:  Log base 2.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        implicit none
        real, intent(in) :: x

        LOG2 = log(x) / log(2.)
    end function

end module
