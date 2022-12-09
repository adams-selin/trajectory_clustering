module cluster_module

    implicit none

CONTAINS

    subroutine calc_Neps(i_seg, epsilon, distance, num_distances, num_segments,&
        num_segIDs, Neps_segID, cluster_distance)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! DESCRIPTION:  Find the segment ID of all line segments within a
    ! distance away of epsilon of line segment L (original segment).
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        use distance_functions, only: indices_1d
        implicit none

        INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(14, 60)

        integer(8), intent(in) :: i_seg
        integer(8), intent(in) :: num_distances, num_segments
        real(dp), intent(in) :: epsilon
        real, intent(in) :: distance(num_distances)
        integer, intent(inout) :: num_segIDs !# segments within epsilon dist.
        integer(8), dimension(20000), intent(inout) :: Neps_segID !segment indices
        real, dimension(20000), intent(inout) :: cluster_distance
            !distance of each Neps_segID away from main segment

        !internal variables
        integer(8) :: i
        integer(8) :: ind1d

        !include line segment L itself in our neighborhood count
        num_segIDs = 1
        Neps_segID(:) = 0
        Neps_segID(1) = i_seg
        cluster_distance(1) = 0.
        do i = 1, num_segments
            if ( i .ne. i_seg ) then
                CALL indices_1d(i,i_seg,num_segments,ind1d)
                if ( (distance(ind1d) .le. epsilon) ) then
                    num_segIDs = num_segIDs + 1
                    Neps_segID(num_segIDs) = i
                    cluster_distance(num_segIDs) = distance(ind1d)
                end if
            end if
        end do

    end subroutine calc_Neps


    subroutine ExpandCluster(Neps_segID, num_segIDs, cluster_id, &
      num_segments, cluster_idx, num_elements_thiscluster, epsilon, &
      MinLns, distance, num_distances)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! DESCRIPTION:  Given the line segments within the epsilon neighborhood
    ! of L, find out if any of those are density-connected. If so add those
    ! to the cluster.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        use distance_functions, only: indices_1d
        implicit none

        integer, intent(in) :: num_segIDs, MinLns
        integer(8), intent(in) :: num_distances, num_segments
        integer(8), dimension(20000), intent(in) :: Neps_segID
        integer, intent(in) :: cluster_idx !cluster ID of line segment L
        real(8), intent(in) :: epsilon
        real, intent(in) :: distance(num_segments,num_segments)
        integer, intent(inout) :: num_elements_thiscluster
        integer, intent(inout) :: cluster_id(num_segments) !cluster id of each
            !segment, including those found within the subroutine

        !list of segment IDs for testing that could be within this cluster
        integer(8), dimension(100000) :: queue
        !number of elements in the queue
        integer :: num_queue
        !number of segments within epsilon neighborhood of segment m
        integer :: num_segm
        !segment IDs of the segments within epsilon neighborhood of segment m
        integer(8), dimension(20000) :: Neps_segm
        !other
        integer :: new_mem_q, i
        integer(8) :: m_id
        real, dimension(20000) :: cluster_distance
            !distance of each segment away from main segment in cluster.
            !Not needed here - routine already checks for classification.

        !queue starts with the line segments within epsilon of L
        num_queue = num_segIDs
        queue(1:num_segIDs) = Neps_segID(1:num_segIDs)

        !check all surrounding line segments
        do while ( num_queue .gt. 0 )
            m_id = queue(num_queue)
            !find line segments within epsilon of this one (L)
            CALL calc_Neps(m_id, epsilon, distance, &
                num_distances,num_segments,num_segm, Neps_segm,cluster_distance)
            !how many new members did we find to check?
            new_mem_q = 0
            if ( num_segm .ge. MinLns ) then
                !check all the new segments we found
                do i = 1, num_segm
                    !were they unclassified or noise?
                    if ( cluster_id(Neps_segm(i)) .le. 0 ) then
                        !if they were unclassified, add to queue for further expansion
                        !as they could be a core segment
                        if ( cluster_id(Neps_segm(i)) .eq. 0 ) then
                            new_mem_q = new_mem_q + 1
                            queue(num_queue+new_mem_q) = Neps_segm(i)
                        end if
                        !classify these points
                        cluster_id(Neps_segm(i)) = cluster_idx
                        num_elements_thiscluster = num_elements_thiscluster + 1
                    end if
                end do !end loop of segments within eps of m
            end if ! end if check for > MinLns

            !remove segment m from queue
            do i = num_queue, num_queue+new_mem_q-1
                queue(i) = queue(i+1)
            end do
            num_queue = num_queue + new_mem_q - 1

        end do !end surrounding line segments check

    end subroutine ExpandCluster

end module cluster_module
