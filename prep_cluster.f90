program prep_cluster
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DESCRIPTION:  Program designed to prepare all necessary data for clustering.
! Specifically, calculate the distance matrix, and determine a range of
! optimal epsilon and MinLns values. Optimal epsilon is determined via
! the Corana et al. (1987) simulated annealing algorithm minimizing Eq. (10)
! of Lee et al.
!
! How to run: ./prep_cluster traj_part_output_directory traj_part_output_filename_prefix
!
! OUTPUT: binary file containing the distance_matrix array, and a text file
! containing recommended starting epsilon and MinLns values for a starting point.
!
! gfortran -fopenmp -o prep_cluster distance_functions.f90 cluster_module.f90 entropy.f90 simulated_anneal.f90 prep_cluster.f90
! gfortran -fopenmp -g -fcheck=all -Wall -fbacktrace -o prep_cluster distance_functions.f90 cluster_module.f90 entropy.f90 simulated_anneal.f90 prep_cluster.f90
! gfortran -g -frecursive -fcheck=all -Wall -fbacktrace -o prep_cluster distance_functions.f90 cluster_module.f90 entropy.f90 simulated_anneal.f90 prep_cluster.f90
! ifort -qopenmp -o prep_cluster_omp.exe cluster_module.f90 entropy.f90 distance_functions.f90 simulated_anneal.f90 prep_cluster.f90
! ifort -o prep_cluster.exe cluster_module.f90 entropy.f90 distance_functions.f90 simulated_anneal.f90 prep_cluster.f90
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use distance_functions, only: distance_matrix
    use simulated_anneal

    implicit none

    real, allocatable :: segments(:,:), distance(:)
    integer(8) :: num_segments
    real(8) :: opt_epsilon
    real :: opt_MinLns
    !character(len=100) :: char_epsilon, char_MinLns

    !distance matrix variables
    character(len=100) :: file_prefix, dir_prefix
    character(len=200) :: prefix
    logical :: file_exists, file_exists_again
    integer :: ios, i
    integer(8) :: idbl
    integer(8):: num_distances, half_num_distances

    !simulated annealing variables (see simulated_anneal module for documentation)
    INTEGER, PARAMETER :: n = 1, neps = 4 !shouldn't be changed
    REAL (dp)   :: lb(n), ub(n), x(n), xopt(n), c(n), vm(n), t, eps, rt, fopt
    INTEGER     :: ns, nt, nfcnev, ier, iseed1, iseed2, maxevl, iprint,  &
                   nacc, nobds
    LOGICAL     :: max




    !read in filename and directory prefixes
    CALL getarg(1,dir_prefix) !should include subtraj directory
    CALL getarg(2,file_prefix)
    prefix = TRIM(dir_prefix)//'/'//TRIM(file_prefix)


    !open output file
    open(50,file=TRIM(prefix)//'_prepcluster.out',form='formatted',&
        status='REPLACE')

    !read the number of line segments from traj_part output
    num_segments = 0
    open(unit=30,file=TRIM(prefix)//'_subtraj.txt',iostat=ios)
    do
        read(30, *, iostat=ios)
        if ( ios /= 0 ) exit
        num_segments = num_segments + 1
    end do
    close(30)

    write(50,*) 'num_segments: ', num_segments

    allocate(segments(7,num_segments))
    num_distances = num_segments * (num_segments-1) / 2
    allocate(distance(num_distances))
    half_num_distances = num_distances / 2


    !read in the traj_part output
    open(unit=30,file=TRIM(prefix)//'_subtraj.txt',iostat=ios)
    do idbl = 1, num_segments
        read(unit=30, fmt="(6(f10.1,1x),1x,f7.0)") segments(:,idbl)
    end do
    close(30)


    !If distance_matrix file doesn't already exist, calculate and store in file.

    !This section assumes only one distance.bin file.
    inquire(file=TRIM(prefix)//'_distance.bin', exist=file_exists)
    if ( file_exists ) then
        write(50,*) TRIM(prefix)//'_distance.bin'
        open(unit=40, file=TRIM(prefix)//'_distance.bin', &
            action="read", form='unformatted', iostat=ios)
        if ( ios /= 0 ) stop "Error opening distance matrix file "
        read(40) (distance(idbl), idbl=1,num_distances)
        close(40)
    else
        stop ('Problem with distance matrix file')
        ! write(50,*) 'calling distance_matrix, will store in '
        ! write(50,*) TRIM(prefix)//'_distance.bin'
        ! CALL distance_matrix(segments, num_segments, distance, num_distances)
        ! open(unit=40, file=TRIM(prefix)//'_distance.bin', action="write",&
        !     form='unformatted')
        ! write(40) (distance(idbl), idbl=1,num_distances)
        ! close(40)
    end if



    write(50,*) maxval(distance)
    write(50,*) minval(distance)


    !Perform the simulated annealing optimization to find the epsilon value
    !associated with the minimum entropy (the optimum).  See the
    !wonderful documentation within the simulated_anneal module, including
    !the original source for this code.
    !  Set underflows to zero on IBM mainframes.
    !     CALL XUFLOW(0)

    ! Set input parameters. It's fine leaving these alone, but if you really
    ! want to muck with them, limit yourself to changing eps, rt, and nt and
    ! leave the rest alone per recommendation from Groffe et al. Again, read
    ! that documentation for suggestions.
    max = .false.  !we want to minimze entropy, not maximize it
    eps = 1.0D-4   !error tolerance for termination
    rt = .5        !temperature reduction factor
    iseed1 = 1     !random number seeds
    iseed2 = 2
    ns = 20        !number of cycles
    nt = 20  !Corona et al. suggested value MAX(100, 5*N), using Goffe et al.
            !recommendation instead
    maxevl = 100000 !max no. of function evaluations. IER=1 if exceeded
    iprint = 1      !set to higher value for more output
    DO i = 1, n
      lb(i) =  1.0D-25  !lower bound of reasonable epsilon values
      ub(i) =  maxval(distance)   !upper bound of reasonable epsilon values
      c(i) = 2.0  !step length adjustment factor.  Don't change.
    END DO

    !Pick a point at which to start epsilon.
    opt_epsilon = 300.
    x(1) = opt_epsilon

    !Set input values of input/output parameters.
    t = 5.0         !starting temperature
    vm(1:n) = 1.0   !step length vector. Will be reset inside fuction, so doesn't matter.


    WRITE(50,1000) n, max, t, rt, eps, ns, nt, neps, maxevl, iprint,iseed1,iseed2

    CALL prtvec(x, n, 'STARTING VALUES')
    CALL prtvec(vm, n, 'INITIAL STEP LENGTH')
    CALL prtvec(lb, n, 'LOWER BOUND')
    CALL prtvec(ub, n, 'UPPER BOUND')
    CALL prtvec(c, n, 'C VECTOR')
    WRITE(50, '(/, "  ****   END OF DRIVER ROUTINE OUTPUT   ****"/,  &
          &     "  ****   before CALL TO sa.             ****")')

    CALL sa(n, x, max, rt, eps, ns, nt, neps, maxevl, lb, ub, c, iprint, &
            iseed1, iseed2, t, vm, xopt, fopt, nacc, nfcnev, nobds, ier, &
            distance, num_segments, opt_MinLns)

    opt_epsilon = xopt(1)

    WRITE(50, '(/, "  ****   RESULTS AFTER SA   ****   ")')
    CALL prtvec(xopt, n, 'SOLUTION')
    CALL prtvec(vm, n, 'FINAL STEP LENGTH')
    WRITE(50,1001) fopt, nfcnev, nacc, nobds, t, ier

    WRITE(50,1002) opt_epsilon, opt_MinLns
    
    close(50)

    1000 FORMAT(/,' SIMULATED ANNEALING EXAMPLE',/,/,  &
                 ' NUMBER OF PARAMETERS: ',i3,'   MAXIMIZATION: ',l5, /, &
                 ' INITIAL TEMP: ', g9.2, '   RT: ',g9.2, '   EPS: ',g9.2, /, &
                 ' NS: ',i3, '   NT: ',i2, '   NEPS: ',i2, /,  &
                 ' MAXEVL: ',i10, '   IPRINT: ',i1, '   ISEED1: ',i4,  &
                 '   ISEED2: ',i4)
    1001 FORMAT(/,' OPTIMAL FUNCTION VALUE: ',g20.13  &
                /,' NUMBER OF FUNCTION EVALUATIONS:     ',i10,  &
                /,' NUMBER OF ACCEPTED EVALUATIONS:     ',i10,  &
                /,' NUMBER OF OUT OF BOUND EVALUATIONS: ',i10,  &
                /,' FINAL TEMP: ', g20.13,'  IER: ', i3)

    1002 FORMAT(/,' OPTIMAL EPSILON: ', g12.5 &
                /,' OPTIMAL MinLns, add 1~3: ', g12.5)



    deallocate(segments)
    deallocate(distance)

end program prep_cluster
