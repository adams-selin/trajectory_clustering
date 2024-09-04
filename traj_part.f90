program traj_part
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DESCRIPTION:  For a set of trajectories, determine the set of characteristic
! points along each trajectory.  Output a file of all subtrajectories (determined
! by those characteristic points), along with the original trajectory source
! number. Based on 2d algorithm provided in Fig. 8 of Lee et al. TRACLUS algorithm.
! Modified by RAS to be 3D. MDL = minimum descriptive length.
!
!  gfortran -o traj_partition distance_functions.f90 MDL.f90 traj_part.f90
!   -I/opt/local/include  -L/opt/local/lib -lnetcdf -lnetcdff
! OR
!  gfortran -c -o distance_functions.o distance_functions.f90
!  gfortran -c -o MDL.o MDL.f90
!  gfortran -o traj_partition traj_part.f90 distance_functions.o MDL.o
!   -I/Users/rselin/miniconda3/include  -L/Users/rselin/miniconda3/lib -lnetcdf
! note - for some reason the above doesn't work anymore.  sigh...
!
! on cheyenne:
!   ifort -o traj_part.exe distance_functions.f90 MDL.f90 traj_part.f90 -I/glade/u/apps/ch/opt/netcdf/4.6.1/intel/17.0.1/include -L/glade/u/apps/ch/opt/netcdf/4.6.1/intel/17.0.1/lib -lnetcdf -lnetcdff
!
! on derecho:
!   ifort -o traj_part.exe distance_functions.f90 MDL.f90 traj_part.f90 -I/glade/u/apps/derecho/23.09/spack/opt/spack/netcdf/4.9.2/oneapi/2023.2.1/yzvj/include -L/glade/u/apps/derecho/23.09/spack/opt/spack/netcdf/4.9.2/oneapi/2023.2.1/yzvj/lib -lnetcdf -lnetcdff
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    use MDL
    use netcdf

    implicit none

    real, allocatable :: pts(:,:)
    real, allocatable :: cpts(:,:)
    real, allocatable :: pts_chars(:,:)
    real, allocatable :: cpts_chars(:,:)  !trajectory characteristics
        !that we will want to pass along through the clustering to provide
        !info at the end about things experienced by the representative traj.
        !Array structure: (num_chars,cplen) where first dimension is the following
        ! (sec,d,dense,ts,fw,vt,ri,rw,tc,w,u,v)
    integer, parameter :: num_chars = 12
    character(len=2) :: str_num_chars
    character(len=40) :: char_file_formatting

    integer :: plen !number of points in this trajectory
    integer :: cplen !number of characteristic points in this trajectory
    integer :: all_cplen !number of char. points in all trajectories
    integer :: startIndex, length, cpIndex, currIndex
    real :: delta !fudge factor to increase length of trajectory partitions by 20-30%
    real :: costPar, costNoPar

    !variables to read in the file
    character(len=100) :: file_prefix, dir_prefix
    character(len=200) :: filename
    character(10) :: thresh_str
    character(6) :: rst_str
    character(20) :: case_name
    integer :: ncid, varid, dimid, ISTAT
    logical :: THERE

    !variables actually in the data
    integer :: numbighail, numtimes
    character(len=40) :: dum !dimension name
    real, allocatable :: bigx(:,:), bigy(:,:), bigz(:,:)
    real, allocatable :: bigd(:,:), bigdense(:,:), bigts(:,:), bigfw(:,:), &
        bigvt(:,:), bigri(:,:), bigrw(:,:), bigtc(:,:), bigsec(:), bigw(:,:), &
        bigu(:,:), bigv(:,:)
    integer :: lastIndex !how many points are in each trajectory

    !variables to write out a file
    character(len=200) :: outfile

    !looping, internal variables
    integer :: i, j, k

    
    !read in command line arguments for the file name
    CALL getarg(1, thresh_str)  !threshold string (e.g., ge15lt19)
    CALL getarg(2, case_name)
    CALL getarg(3, rst_str)
    file_prefix= 'haildata_'//TRIM(case_name)//'_'//TRIM(thresh_str)//'_Wrel'
    dir_prefix = '/glade/derecho/scratch/radams/'//TRIM(case_name)//'/rst'//TRIM(rst_str)
    
    !read in data
    !filename = TRIM(dir_prefix)//'/'//TRIM(file_prefix)//'.nc'
    filename = TRIM(dir_prefix)//'/haildata_'//TRIM(thresh_str)//'_Wrel.nc'
    print *, filename
    INQUIRE( FILE=TRIM(filename), EXIST=THERE)
    IF ( THERE ) THEN
        call check( nf90_open(filename, NF90_NOWRITE,ncid) )
        call check( nf90_inq_dimid(ncid, 'traj', dimid))
        call check( nf90_inquire_dimension(ncid, dimid, dum, numbighail))
        call check( nf90_inq_dimid(ncid, 'time', dimid))
        call check( nf90_inquire_dimension(ncid, dimid, dum, numtimes))
        print *, numbighail, numtimes
        !allocate our variables
        allocate(bigx(numtimes, numbighail))
        allocate(bigy(numtimes, numbighail))
        allocate(bigz(numtimes, numbighail))
        allocate(bigd(numtimes, numbighail))
        allocate(bigdense(numtimes, numbighail))
        allocate(bigts(numtimes, numbighail))
        allocate(bigfw(numtimes, numbighail))
        allocate(bigvt(numtimes, numbighail))
        allocate(bigri(numtimes, numbighail))
        allocate(bigrw(numtimes, numbighail))
        allocate(bigtc(numtimes, numbighail))
        allocate(bigw(numtimes, numbighail))
        allocate(bigu(numtimes, numbighail))
        allocate(bigv(numtimes, numbighail))
        allocate(bigsec(numtimes))
        !now read them in
        call check( nf90_inq_varid(ncid, 'x', varid))
        call check( nf90_get_var(ncid, varid, bigx))
        call check( nf90_inq_varid(ncid, 'y', varid))
        call check( nf90_get_var(ncid, varid, bigy))
        call check( nf90_inq_varid(ncid, 'z', varid))
        call check( nf90_get_var(ncid, varid, bigz))
        call check( nf90_inq_varid(ncid, 'd', varid))
        call check( nf90_get_var(ncid, varid, bigd))
        call check( nf90_inq_varid(ncid, 'dense', varid))
        call check( nf90_get_var(ncid, varid, bigdense))
        call check( nf90_inq_varid(ncid, 'ts', varid))
        call check( nf90_get_var(ncid, varid, bigts))
        call check( nf90_inq_varid(ncid, 'fw', varid))
        call check( nf90_get_var(ncid, varid, bigfw))
        call check( nf90_inq_varid(ncid, 'tv', varid))
        call check( nf90_get_var(ncid, varid, bigvt))
        call check( nf90_inq_varid(ncid, 'qice', varid))
        call check( nf90_get_var(ncid, varid, bigri))
        call check( nf90_inq_varid(ncid, 'qliq', varid))
        call check( nf90_get_var(ncid, varid, bigrw))
        call check( nf90_inq_varid(ncid, 'tc', varid))
        call check( nf90_get_var(ncid, varid, bigtc))
        call check( nf90_inq_varid(ncid, 'w', varid))
        call check( nf90_get_var(ncid, varid, bigw))
        call check( nf90_inq_varid(ncid, 'u', varid))
        call check( nf90_get_var(ncid, varid, bigu))
        call check( nf90_inq_varid(ncid, 'v', varid))
        call check( nf90_get_var(ncid, varid, bigv))
        call check( nf90_inq_varid(ncid, 'time', varid))
        call check( nf90_get_var(ncid, varid, bigsec))
        ISTAT = NF90_CLOSE(ncid)
    ELSE
        stop "hail trajectory file doesn't exist" 
    ENDIF

    !get rid of the extra 10^9 factor that made it into all the Wrel.nc times
    bigsec = bigsec * 1E-9

    !Delta = small constant to encourage non-partitioning.  Helps with clustering later.
    delta = 1. !tested on first trajectory. This value reduced number of char. pts by ~25%,
               !following recommendation from Lee et al. top of page 599.
    delta = 17. ! Higher value to make distance_matrix manageable

    !Initialize all points counter
    all_cplen = 0

    !Open output files
    !change file_prefix here to make it smaller, if you want to
    !file_prefix = TRIM(case_name)//'_'//TRIM(rst_str)//'_Wrelative_srwinds'
    outfile = TRIM(dir_prefix)//'/subtraj/'//TRIM(file_prefix)//'_subtraj.txt'
    open(30,file=TRIM(outfile),form='formatted', status='REPLACE')
    outfile = TRIM(dir_prefix)//'/subtraj/'//TRIM(file_prefix)//'_subtraj_chars.txt'
    open(31,file=TRIM(outfile),form='formatted', status='REPLACE')

    !Set up a formatting string to output double the number of available characteristics
    write(str_num_chars,'(i2)') num_chars-1
    char_file_formatting = "(f6.0,"//TRIM(str_num_chars)//"(f12.4,1x),"//&
                           "f6.0,"//TRIM(str_num_chars)//"(f12.4,1x))"

    !Start the loop through each trajectory
    do i = 1, numbighail
        print *, ''
        print *, 'Trajectory: ', i

        !put all the non-missing parts of this trajectory into the pts array
        !plen = FINDLOC(bigx(:,i), -9999., dim=1) - 1
        !recode so not using FINDLOC, so can use more compilers
        plen = 0
        do j = 1, numtimes
            if (bigx(j,i) .eq. -9999.) then
                exit
            endif
            plen = j
        enddo
        !don't waste our time if it isn't at least 10 pts long
        if ( plen .ge. 10 ) then

            !assign everything to our pts or pts_chars array
            allocate(pts(3,plen))
            allocate(cpts(3,plen))
            allocate(pts_chars(num_chars,plen))
            allocate(cpts_chars(num_chars,plen))
            do j = 1, plen
                pts(1,j) = bigx(j,i)
                pts(2,j) = bigy(j,i)
                pts(3,j) = bigz(j,i)
                pts_chars(1,j) = bigsec(j)
                pts_chars(2,j) = bigd(j,i)
                pts_chars(3,j) = bigdense(j,i)
                pts_chars(4,j) = bigts(j,i)
                pts_chars(5,j) = bigfw(j,i)
                pts_chars(6,j) = bigvt(j,i)
                pts_chars(7,j) = bigri(j,i)
                pts_chars(8,j) = bigrw(j,i)
                pts_chars(9,j) = bigtc(j,i)
                pts_chars(10,j) = bigw(j,i)
                pts_chars(11,j) = bigu(j,i)
                pts_chars(12,j) = bigv(j,i)
            end do !end pts assignment loop
            print *, ' plen: ', plen

            !starting point is first characteristic point
            !Note - check first point - may have to skip because it doesn't have defined
            !updraft, etc.
            !cpts(:,1) = pts(:,2)
            !cpts_chars(:,1) = pts_chars(:,2)
            !startIndex = 2 !last characteristic point index in pts
            
            !Or, if you everything is defined at point 1:
            cpts(:,1) = pts(:,1)
            cpts_chars(:,1) = pts_chars(:,1)
            startIndex = 1 !last characteristic point index in pts
            length = 2 !number of trajectory points since startIndex
                       !start at two so you have at least two non-char pts to compare
            cpIndex = 2 !keep track of where we are in the char pts array

            do while ( (startIndex + length) .le. plen )
                currIndex = startIndex + length
                CALL MDLpar(length+1, pts(:,startIndex:currIndex), costPar)
                CALL MDLnopar(length+1, pts(:,startIndex:currIndex), costNoPar)
                costNoPar = costNoPar + delta

                !if partitioning at current point makes the MDL cost larger than not,
                ! partition at the previous point
                if ( costPar .gt. costNoPar ) then
                    cpts(:,cpIndex) = pts(:,currIndex-1)

                    !trajectory characteristics are averaged over this new partition
                    cpts_chars(:,cpIndex) = &
                        SUM(pts_chars(:,startIndex:currIndex-1), dim=2) / &
                        (currIndex-1 - startIndex + 1)
                    !print *, ' cpIndex, ave sec, interval: ', &
                    !    cpIndex, cpts_chars(1,cpIndex), startIndex, currIndex-1
                    !print *, '  pts_chars(1,startIndex:currIndex-1): ', &
                    !    pts_chars(1,startIndex:currIndex-1)


                    startIndex = currIndex - 1
                    length = 2
                    cpIndex = cpIndex + 1
                else
                    length = length + 1
                endif
            end do !end while loop through subtrajectory

            !include last point in characteristic points
            cpts(:,cpIndex) = pts(:,plen)
            cpts_chars(:,cpIndex) = pts_chars(:,plen)
            cplen = cpIndex
            all_cplen = all_cplen + cplen

            print *, '  num char. pts: ', cplen

            !write data out to file as a series of subtrajectory segments
            do j = 1, cplen-1
                write(unit=30, fmt="(6(f10.1,1x),1x,i7)") &
                    cpts(:,j),cpts(:,j+1),i
                !write(unit=31, fmt="(20(f12.4,1x))") &
                !    cpts_chars(:,j),cpts_chars(:,j+1)
                write(unit=31, fmt=TRIM(char_file_formatting)) &
                    cpts_chars(:,j),cpts_chars(:,j+1)
            end do

            if (allocated(pts)) deallocate(pts)
            if (allocated(cpts)) deallocate(cpts)
            if (allocated(pts_chars)) deallocate(pts_chars)
            if (allocated(cpts_chars)) deallocate(cpts_chars)
        end if !end check for long enough trajectory

    end do !end loop of all hail trajectories

    !close output files
    close(unit=30)
    close(unit=31)

    !deallocate all the netcdf data
    if (allocated(bigx)) deallocate(bigx)
    if (allocated(bigy)) deallocate(bigy)
    if (allocated(bigz)) deallocate(bigz)
    if (allocated(bigd)) deallocate(bigd)
    if (allocated(bigdense)) deallocate(bigdense)
    if (allocated(bigts)) deallocate(bigts)
    if (allocated(bigfw)) deallocate(bigfw)
    if (allocated(bigvt)) deallocate(bigvt)
    if (allocated(bigri)) deallocate(bigri)
    if (allocated(bigrw)) deallocate(bigrw)
    if (allocated(bigtc)) deallocate(bigtc)
    if (allocated(bigw)) deallocate(bigw)
    if (allocated(bigu)) deallocate(bigu)
    if (allocated(bigv)) deallocate(bigv)
    if (allocated(bigsec)) deallocate(bigsec)

    !print out total number of subtrajectory segments
    print *, '  num_segments: ', all_cplen

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! SUBROUTINE CHECK : Checks NETCDF error status, stops program if there is an error.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE check(status)

  integer, intent (in) :: status

  if(status /= nf90_noerr) then
     print *, trim(nf90_strerror(status))
     stop 2
  end if

  END SUBROUTINE check


end program
