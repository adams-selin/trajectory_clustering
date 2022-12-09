program cluster
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DESCRIPTION:  Given a list of subtrajectory line segments output from
! traj_part, calculate a set of line segment clusters. Following the line
! segment clustering pseudocode provided in Lee et al., Fig. 12.
! First run through will create a distance matrix binary file.
! Line segment clusters are in the following format:
! segments(7,num_segments) = ( (x1,y1,z1, x2,y2,z2, orig_trajnum); num_segments)
!
! How to run: ./cluster traj_part_output_dir traj_part_output_filename_prefix epsilon MinLns
!
! OUTPUT: multiple files, one for each cluster, each containing all segments
! within that cluster.
!
! gfortran -g -fcheck=all -Wall -fbacktrace -frecursive -o cluster distance_functions.f90 cluster_module.f90 cluster.f90
!  gfortran -o cluster distance_functions.f90 cluster_module.f90 cluster.f90
! ifort -qopenmp -o cluster distance_functions.f90 cluster_module.f90 cluster.f90
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    use cluster_module

    implicit none

    real, allocatable :: segments(:,:), distance(:)
    real, allocatable :: segments_chars(:,:)  !storm characteristics of each segment
    integer(8) :: num_segments
    real :: epsilon
    integer ::  MinLns
    !character(len=100) :: char_epsilon, char_MinLns

    !distance matrix variables
    character(len=100) :: dir_prefix, file_prefix
    character(len=200) :: input_prefix, output_prefix, logfile
    character(len=4) :: char_MinLns
    character(len=8) :: char_epsilon
    integer :: ios, k
    integer(8) :: i, j
    !logical :: file_exists, file_exists_again

    !variables we produce
    integer :: num_clusters
    integer(8):: num_distances, half_num_distances
    real, allocatable :: segments_cluster(:,:) !all segments within one cluster
    real, allocatable :: segments_cluster_chars(:,:) !segments characteristics
                                                     !within one cluster
    integer :: segments_cluster_idx
    integer :: num_segIDs !number of segments less than epsilon away
    !arrays we have to define as having an arbitrarily large dimension
    integer(8), dimension(100000) :: Neps_segID !indices of segments array less
                                                !than epsilon away
    integer, dimension(20000) :: num_elements_cluster
    real, dimension(20000) :: cluster_distance !distance of each Neps_segID away
                                               !from main segment
    !allocatable arrays we will know the size of before we need them
    integer :: cluster_idx !cluster ID index
    integer, allocatable :: cluster_id(:) !cluster number assigned each segment
    real, allocatable :: min_cluster_distance(:) !min distance from main
        !segment in cluster; only assign cluster id if new distance is smaller
    integer(8), allocatable :: trajnums_cluster(:) !all orig_trajnums from the
                                                   !segments in this cluster
    integer :: MinTrajCard !min value for trajectory cardinality of the cluster
    integer :: unique
    integer(8) :: min_val, diff, min_trajnums_cluster, max_trajnums_cluster
    character(len=4) :: str_i

    !for counting noise
    integer(8) :: noise_count


    !read in filename prefix
    CALL getarg(1,dir_prefix)
    CALL getarg(2,file_prefix)
    CALL getarg(3,char_epsilon)
    CALL getarg(4,char_MinLns)
    !print *, char_epsilon, char_MinLns
    read (char_epsilon, *) epsilon
    read (char_MinLns, *) MinLns

    input_prefix = TRIM(dir_prefix)//'/subtraj/'//TRIM(file_prefix)
    output_prefix = TRIM(dir_prefix)//'/cluster/'//TRIM(file_prefix)//&
        '_'//TRIM(char_epsilon)//'_'//TRIM(char_MinLns)


    if ( MinLns .lt. 8 ) then
        MinTrajCard = MinLns
    else
        MinTrajCard = 8 !MinLns !this could change
    end if

    !create an output log file
    logfile = TRIM(output_prefix)//'_cluster.log'
    open(unit=110, file=logfile, status='REPLACE', form='formatted', &
        action='write')

    write(110,*) 'Starting cluster.'
    write(110,*) 'epsilon, MinLns, MinTrajCard: ', epsilon, MinLns, MinTrajCard

    !read the number of line segments from traj_part output
    num_segments = 0
    write(110,*) 'Reading in '//TRIM(input_prefix)//'_subtraj.txt'
    open(unit=30, file=TRIM(input_prefix)//'_subtraj.txt', iostat=ios)
    do
        read(30, *, iostat=ios)
        if ( ios /= 0 ) exit
        num_segments = num_segments + 1
    end do
    close(30)

    write(110,*) 'num_segments: ', num_segments

    allocate(segments(7,num_segments))
    allocate(segments_chars(20,num_segments))
    num_distances = num_segments * (num_segments-1) / 2
    allocate(distance(num_distances))
    half_num_distances = num_distances / 2


    !read in the traj_part output
    open(unit=30, file=TRIM(input_prefix)//'_subtraj.txt', iostat=ios)
    do i = 1, num_segments
        read(unit=30, fmt="(6(f10.1,1x),1x,f7.0)") segments(:,i)
    end do
    close(30)
    !read in the traj_part characteristic output
    open(unit=31, file=TRIM(input_prefix)//'_subtraj_chars.txt', iostat=ios)
    do i = 1, num_segments
        read(unit=31, fmt="(20(f12.4,1x))") segments_chars(:,i)
    end do
    close(31)


    !Read in the distance matrix calculated in prep_cluster. This section assumes
    !just one distance.bin file.
    write(110,*) 'reading distance matrix in from '
    write(110,*) TRIM(input_prefix)//'_distance.bin'
    open(unit=40,file=TRIM(input_prefix)//'_distance.bin',action="read",&
        form='unformatted', iostat=ios)
    if ( ios /= 0 ) stop "Error opening distance matrix file "
    read(40) (distance(i), i=1,num_distances)
    close(40)

    !This section assumes only two distance.bin files, with half the distance array in each.
    ! inquire(file=TRIM(input_prefix)//'_distance1.bin', exist=file_exists)
    ! inquire(file=TRIM(input_prefix)//'_distance2.bin', exist=file_exists_again)
    ! if ( file_exists .and. file_exists_again ) then
    !     write(50,*) TRIM(input_prefix)//'_distance1.bin'
    !     open(unit=40, file=TRIM(input_prefix)//'_distance1.bin', &
    !         action="read", form='unformatted', iostat=ios)
    !     if ( ios /= 0 ) stop "Error opening 1st distance matrix file "
    !     read(40) (distance(i), i=1,half_num_distances)
    !     close(40)
    !
    !     write(50,*) TRIM(input_prefix)//'_distance2.bin'
    !     open(unit=60, file=TRIM(input_prefix)//'_distance2.bin', &
    !         action="read", form='unformatted', iostat=ios)
    !     if ( ios /= 0 ) stop "Error opening 2nd distance matrix file "
    !     read(60) (distance(i), i=half_num_distances+1,num_distances)
    !     close(60)
    ! else
    !     stop "No distrance matrix files"
    ! end if


    write(110,*)'max distance: ', maxval(distance)
    close(110)
    open(unit=110, file=logfile, status='OLD', form='formatted', &
        action='write',position='append')

    !declare arrays we'll need for clustering
    allocate(cluster_id(num_segments))
    cluster_id(:) = 0 !initial cluster ID (i.e., unclassified). -1 if noise.
    cluster_idx = 1 !cluster ID index, starts at 1.
    !Need to keep track of minimum distance from main line segment within
    ! cluster, and only assigned cluster id if new distance is smaller
    allocate(min_cluster_distance(num_segments))
    min_cluster_distance(:) = 9999999.


    !loop through all segments
    num_elements_cluster(:) = 0
    do i = 1, num_segments
        if ( cluster_id(i) .eq. 0 ) then
            !find line segments within epsilon of this one (L)
            CALL calc_Neps(i, REAL(epsilon,kind=8), distance, num_distances, &
                num_segments, num_segIDs, Neps_segID, cluster_distance)

            !if we have at least MinLns segments, assign them to a cluster
            if ( num_segIDs .ge. MinLns ) then
                do j = 1, num_segIDs
                    !Has a segment already been assigned a cluster? If so,
                    !check how far away it was, and only change its assignment
                    !if its distance here is smaller.
                    if ( cluster_id(Neps_segID(j)) .gt. 0 ) then
                        if ( cluster_distance(j) .lt. &
                             min_cluster_distance(Neps_segID(j)) ) then
                            !remove this segment from the previous cluster
                            num_elements_cluster(cluster_id(Neps_segID(j))) = &
                              num_elements_cluster(cluster_id(Neps_segID(j)))-1
                            !add this segment to the current cluster
                            cluster_id(Neps_segID(j)) = cluster_idx
                            num_elements_cluster(cluster_idx) = &
                                num_elements_cluster(cluster_idx) + 1
                            min_cluster_distance(Neps_segID(j)) = &
                                cluster_distance(j)
                        else
                            cycle
                        endif
                    else
                        cluster_id(Neps_segID(j)) = cluster_idx
                        num_elements_cluster(cluster_idx) = &
                            num_elements_cluster(cluster_idx) + 1
                    end if
                end do

                !print *, 'cluster_idx, i: ', cluster_idx, i
                !print *, '  Neps_segID: ', Neps_segID(1:num_segIDs)
                !print *, '  num_elements_cluster: ', num_elements_cluster(cluster_idx)


                !find all density-connected sets to this segment
                CALL ExpandCluster(Neps_segID, num_segIDs, cluster_id, &
                  num_segments, cluster_idx, num_elements_cluster(cluster_idx),&
                  REAL(epsilon,kind=8), MinLns, distance, num_distances)
                cluster_idx = cluster_idx + 1
            else
                !mark as noise
                cluster_id(i) = -1
            end if !end if statement for at least MinLns segments in cluster
        end if !end if statement for unclassified line segment
    end do !end do loop through all segments
    !added 1 to cluster_idx at end of loop, so subtract for real value
    cluster_idx = cluster_idx - 1
    write(110,*) 'done classifying line segments'
    write(110,*) 'identified ', cluster_idx, ' clusters'
    print *, cluster_idx

    !determine how many segments have been classified as noise
    noise_count = 0
    do i = 1, num_segments
        if ( cluster_id(i) .lt. 0) then
            noise_count = noise_count + 1
        end if
    end do
    write(110,*) noise_count, ' segments are noise'
    print *, noise_count


    !If there are enough segements in a cluster, output to a file
    num_clusters = 0
    !loop through each cluster
    do i = 1, cluster_idx
        write(110,*) ''
        write(110,*) i, ' of ', cluster_idx, ' clusters'
        write(110,*) num_elements_cluster(i), ' segments in this cluster'
        allocate(segments_cluster(7,num_elements_cluster(i)))
        allocate(segments_cluster_chars(20,num_elements_cluster(i)))
        segments_cluster_idx = 1
        do j = 1, num_segments
            !pull out all segments in the current cluster
            if ( cluster_id(j) .eq. i ) then
                segments_cluster(:,segments_cluster_idx) = segments(:,j)
                segments_cluster_chars(:,segments_cluster_idx) = &
                    segments_chars(:,j)
                segments_cluster_idx = segments_cluster_idx + 1
            end if
        end do !end segment loop

        !how many different trajectories did these cluster segments come from?
        write(110,*) 'determining unique trajectories'
        if ( num_elements_cluster(i) .eq. 0 ) then
            write(110,*) 'cluster is empty, next'
            deallocate(segments_cluster)
            deallocate(segments_cluster_chars)
            deallocate(trajnums_cluster)
            cycle !next cluster
        end if
        allocate(trajnums_cluster(num_elements_cluster(i)))
        trajnums_cluster(:) = INT(segments_cluster(7,:), kind=8)
        min_trajnums_cluster = minval(trajnums_cluster)
        max_trajnums_cluster = maxval(trajnums_cluster)
        diff = max_trajnums_cluster - min_trajnums_cluster
        if ( diff .lt. MinTrajCard ) then
            write(110,*) 'fewer than ', diff, ' trajectories'
            deallocate(segments_cluster)
            deallocate(segments_cluster_chars)
            deallocate(trajnums_cluster)
            cycle  !go to next cluster
        end if
        !fortran "unique" code from
        !https://stackoverflow.com/questions/44198212/a-fortran-equivalent-to-unique
        unique = 0
        min_val = min_trajnums_cluster - 1
        do while ( min_val .lt. max_trajnums_cluster )
            unique = unique + 1
            min_val = minval(trajnums_cluster, &
                      mask=(trajnums_cluster .gt. min_val))
        end do
        write(110,*) 'found ', unique, ' unique trajectories'
        write(110,*) 'greater than ', MinTrajCard, '?'

        !if we have enough unique trajectories in this cluster, write out
        if ( unique .gt. MinTrajCard ) then
            !print *, 'writing out file ', i
            write (str_i, '(i0)') i
            write(110,*) 'writing to ', TRIM(output_prefix)//'_cluster_'//&
                TRIM(str_i)//'.txt'
            !write out the trajectories themselves
            open(100,file=TRIM(output_prefix)//'_cluster_'//TRIM(str_i)&
                //'.txt',status="replace", action="write")
            write(100,'(I0)') num_elements_cluster(i)
            do k = 1, num_elements_cluster(i)
              write(100,"(7(f10.1,1x))") (segments_cluster(j,k), j=1,7)
            end do

            !write out the trajectory characteristics
            open(101,file=TRIM(output_prefix)//'_cluster_chars_'//TRIM(str_i)&
                //'.txt',status="replace", action="write")
            write(101,'(I0)') num_elements_cluster(i)
            do k = 1, num_elements_cluster(i)
              write(101,"(20(f12.4,1x))") (segments_cluster_chars(j,k), j=1,20)
            end do
            num_clusters = num_clusters + 1
        end if

        !deallocate arrays for next cluster iteration
        deallocate(segments_cluster)
        deallocate(segments_cluster_chars)
        deallocate(trajnums_cluster)

    end do

close(110)
close(100)
close(101)

end program cluster
