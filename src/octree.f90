module octree_mod

	use iso_c_binding

  implicit none
! 
!   public
! 
!   public point_type
!   public octree_init
!   public octree_final
!   public octree_build
!   public octree_update
!   public octree_search
!   public print_tree
!   public print_node
!   public get_all_depths
!   public get_all_leaves
!   public find_leaves
!   public find_all_neighbors
!   public find_all_neighbors_recursive
!   public compute_all_vertices
!   public compute_vertices_from_bbox
!   public check_if_on_boundary
!   public check_if_they_share_vertex
!   public clean_up_junk_tree
!   public write_octree
!   public clean_node
!   public clean_tree
!   public check_if_lie_on_edge
!   public check_if_lie_on_plane
!   public copy_node
!   
!   public tree
!   public config
!   public node_type
!   public nodes_type
!   public tree_type
!   public node_list_type
! 

  type config_type
    integer max_num_point ! Maximum point number contained in leaf node.
    integer max_depth     ! Maximum level of branch and leaf nodes
    real(8) bbox(2, 3)
  end type config_type

  ! Points should be indexed by their id.
  type point_type
    integer id
    real(8) x(3)
  end type point_type

  ! There are two kinds of nodes:
  !   1. Branch node with children;
  !   2. Leaf node without child but containing points.
  type node_type
    integer depth, i, node_id, cell, leaf_id
    integer child_num
    real(8) bbox(2, 3)
    real(8) vertices(8,3) 
    real(8) :: centroid(3) = 0.0d0
    real(8) :: length = 0.0d0
    integer :: num_point = 0
    integer :: nneighbor = 0
    logical :: computed_moment = .false. 
    logical :: leaf = .false. 
    integer, allocatable :: point_ids(:)
    type(node_type), pointer :: parent
    type(node_type), pointer :: children(:)
    type(node_type), pointer :: neighbors(:)

  end type node_type

  type nodes_type
  	type(node_type), pointer :: node
  end type nodes_type  

  type (nodes_type), allocatable :: tmp_nodes(:)

  type node_list_type
	integer :: ncell = 0
	integer :: nonempty = 0
	type(nodes_type), allocatable :: nodes(:)
  end type node_list_type


  type tree_type
	 integer, pointer :: ids(:)
	 real(8), pointer :: points(:,:)
	 integer :: nleaves = 0
	 type(nodes_type), allocatable :: leaves(:)
     type(node_type), pointer :: root_node
     type(node_list_type), allocatable :: depth_list(:)
     integer :: nnodes
     integer :: max_pts
     integer :: max_depth
  end type tree_type

  type(config_type) config
  type(tree_type) tree
  
  logical, allocatable :: flag(:,:)
  
contains









  subroutine octree_init(max_num_point, max_depth, bbox)

    integer, intent(in), optional :: max_num_point
    integer, intent(in), optional :: max_depth
    real(8), intent(in), optional :: bbox(2, 3)

    
    if(present(max_num_point)) then
    	config%max_num_point = max_num_point
    else
	     config%max_num_point = 3
    end if

    if(present(max_depth)) then
	     config%max_depth = max_depth
    else
	     config%max_depth = 10
    end if
    
    if(.not. allocated(tree%depth_list)) then
    	allocate( tree%depth_list(0:config%max_depth) )
    end if
    
    
    
    if(present(bbox)) then
    	config%bbox = bbox
    else
	    config%bbox =  reshape([0.0d0, 1.0d0, 0.0d0, 1.0d0, 0.0d0, 1.0d0], [2, 3])
    end if
    
   if (.not. associated(tree%root_node)) allocate(tree%root_node)
    call reset_node(tree%root_node)
    tree%root_node%depth = 0
    tree%root_node%child_num = 1
    tree%root_node%bbox = config%bbox
	tree%nleaves = 0
	tree%nnodes = 0
	tree%max_pts = config%max_num_point
	tree%max_depth = 0


  end subroutine octree_init

  subroutine octree_final()

    call clean_node(tree%root_node)
    deallocate(tree%root_node)

  end subroutine octree_final

  recursive subroutine octree_build(ids, points, node_)

!--	Sohail Reddy
	integer, intent(in), target :: ids(:)
	real(8), intent(in), target :: points(:,:)
!--	
	
    type(node_type), intent(inout), target, optional :: node_

    type(node_type), pointer :: node
    integer i, j, counter
    integer num_contained_point, depth
    
!--	Sohail Reddy
	 integer :: tmp_index(size(points,1))
	 integer, allocatable :: contained_ids(:)
	 real(8), allocatable :: contained_points(:,:)
!--

    if (present(node_)) then
      node => node_
    else
      tree%points => points
      tree%ids => ids
      node => tree%root_node
    end if
    
   if(.not. allocated(flag)) then
   	allocate(flag(0:config%max_depth, size(points,1)))
   	flag = .false. 
   end if


	depth = node%depth
    tree%depth_list(depth)%ncell = tree%depth_list(depth)%ncell + 1
	tree%nnodes = tree%nnodes + 1
	node%node_id = tree%nnodes
	node%leaf_id = 0
	
	
    num_contained_point = 0
	tmp_index = 0

    do i = 1, size(points,1)
    
      if (points(i,1) < node%bbox(1, 1) .or. points(i,1) > node%bbox(2, 1) .or. &
          points(i,2) < node%bbox(1, 2) .or. points(i,2) > node%bbox(2, 2) .or. &
          points(i,3) < node%bbox(1, 3) .or. points(i,3) > node%bbox(2, 3)) cycle
 
!       if (points(i,1) >= node%bbox(1, 1) .and. points(i,1) < node%bbox(2, 1) .and. &
!          points(i,2) >= node%bbox(1, 2) .and. points(i,2) < node%bbox(2, 2) .and. &
!          points(i,3) >= node%bbox(1, 3) .and. points(i,3) < node%bbox(2, 3)) then

   
    	if(.not. flag(depth,ids(i))) then
    	  num_contained_point = num_contained_point + 1            
    	  tmp_index(num_contained_point) = i
	      flag(depth,ids(i)) = .true.    	  
    	end if
    	
!    	end if
    end do
	
	
	if (num_contained_point == 0) then
		node%leaf = .true.
     	tree%nleaves = tree%nleaves + 1
     	return
	end if
	
	if(allocated(node%point_ids)) then
		deallocate(node%point_ids)
		allocate(node%point_ids(num_contained_point))	
	else
		allocate(node%point_ids(num_contained_point))
	end if
	node%num_point = num_contained_point
	node%point_ids = ids(tmp_index(1:num_contained_point))
	
	! its a leaf
	if( node%num_point <= config%max_num_point .or. &
		node%depth >= config%max_depth) then
		node%leaf = .true.
     	tree%nleaves = tree%nleaves + 1		
	 	return
    end if
	



    allocate(contained_ids(num_contained_point))
    allocate(contained_points(num_contained_point,3))

	do i = 1, num_contained_point
		contained_ids(i) = ids(tmp_index(i))
		contained_points(i,:) = points(tmp_index(i),:)
	end do
	

    ! Subdivide node and run into the child nodes.
    call subdivide_node(node)
    do i = 1, size(node%children)
	  node%children(i)%child_num = i
      call octree_build(contained_ids,contained_points, node%children(i))
    end do

	return
	
  end subroutine octree_build

  subroutine octree_update(node_)

    type(node_type), intent(inout), target, optional :: node_

    type(node_type), pointer :: node

    if (present(node_)) then
      node => node_
    else
      node => tree%root_node
    end if

  end subroutine octree_update

  recursive subroutine octree_search(x, distance, cnode, node_)

    real(8), intent(in) :: x(3)
    real(8), intent(in) :: distance
	type(node_type), pointer, optional :: cnode
    type(node_type), intent(in), target, optional :: node_

    type(node_type), pointer :: node
    real(8) d2, dx(3)
    integer i

    if (present(node_)) then
      node => node_
    else
      node => tree%root_node
    end if

    if (associated(node%children)) then
      ! We are at branch node.
      do i = 1, size(node%children)
        if ((x(1) + distance) > node%children(i)%bbox(1,1) .and. &
            (x(1) - distance) < node%children(i)%bbox(2,1) .and. &
            (x(2) + distance) > node%children(i)%bbox(1,2) .and. &
            (x(2) - distance) < node%children(i)%bbox(2,2) .and. &
            (x(3) + distance) > node%children(i)%bbox(1,3) .and. &
            (x(3) - distance) < node%children(i)%bbox(2,3)) then

!        if ((x(1) + distance) > node%children(i)%bbox(1,1) .and. &
!            (x(1) - distance) < node%children(i)%bbox(2,1) .and. &
!            (x(2) + distance) > node%children(i)%bbox(1,2) .and. &
!            (x(2) - distance) < node%children(i)%bbox(2,2) ) then
          call octree_search(x, distance, cnode = cnode, node_ = node%children(i))
        end if
      end do
    else
      ! We are at leaf node.
 		cnode => node
    end if

  end subroutine octree_search





  subroutine reset_node(node)

    type(node_type), intent(inout) :: node

    node%num_point = 0
    if (.not. allocated(node%point_ids)) allocate(node%point_ids(config%max_num_point))
    nullify(node%parent)
    nullify(node%children)
    nullify(node%neighbors)


  end subroutine reset_node

  subroutine subdivide_node(node)

    type(node_type), intent(inout), target :: node

    integer i, j, k, l
    real(8) bbox(2, 3)

    allocate(node%children(4))
    l = 1
!    do k = 1, 2
      do j = 1, 2
        do i = 1, 2
          call reset_node(node%children(l))
          node%children(l)%depth = node%depth + 1
          node%children(l)%parent => node
          node%children(l)%bbox(1, 1) = node%bbox(1, 1) + (i - 1) * (node%bbox(2, 1) - node%bbox(1, 1)) * 0.5d0
          node%children(l)%bbox(2, 1) = node%bbox(2, 1) - (2 - i) * (node%bbox(2, 1) - node%bbox(1, 1)) * 0.5d0
          node%children(l)%bbox(1, 2) = node%bbox(1, 2) + (j - 1) * (node%bbox(2, 2) - node%bbox(1, 2)) * 0.5d0
          node%children(l)%bbox(2, 2) = node%bbox(2, 2) - (2 - j) * (node%bbox(2, 2) - node%bbox(1, 2)) * 0.5d0
          node%children(l)%bbox(1, 3) = node%bbox(1, 3) !+ (k - 1) * (node%bbox(2, 3) - node%bbox(1, 3)) * 0.5d0
          node%children(l)%bbox(2, 3) = node%bbox(2, 3) !- (2 - k) * (node%bbox(2, 3) - node%bbox(1, 3)) * 0.5d0
          node%children(l)%parent => node


		  node%children(l)%centroid = ((node%children(l)%bbox(2, :) + node%children(l)%bbox(1, :))/2.0d0) 
		  node%children(l)%length = node%children(l)%bbox(2, 1) - node%children(l)%bbox(1, 1)
		  
		  call compute_vertices_from_bbox(node%children(l)%bbox, node%children(l)%vertices)          
          l = l + 1
        end do
      end do
 !   end do



  end subroutine subdivide_node

  recursive subroutine clean_node(node)

    type(node_type), intent(inout) :: node

    integer i

    if (associated(node%children)) then
      do i = 1, size(node%children)
        call clean_node(node%children(i))
        deallocate(node%children(i)%point_ids)
      end do
      deallocate(node%children)
    end if
    
	if(associated(node%neighbors)) nullify(node%neighbors)

  end subroutine clean_node

	subroutine clean_tree
	
		implicit none
		if(.not. associated(tree%root_node)) return
		call clean_node(tree%root_node)
		if(allocated(tree%depth_list)) deallocate(tree%depth_list)
		if(allocated(tree%leaves)) deallocate(tree%leaves)

		if(associated(tree%ids)) nullify(tree%ids)
		if(associated(tree%points)) nullify(tree%points)

	end subroutine clean_tree

	subroutine clean_up_junk_tree
		implicit none
		
		if(allocated(flag)) deallocate(flag)
		if(allocated(tmp_nodes)) deallocate(tmp_nodes)

	end subroutine clean_up_junk_tree

  subroutine compute_vertices_from_bbox(bbox, vertices)

	implicit none
	real(8) :: bbox(2,3), vertices(8,3), vector(3)
	integer i, j, k, l
	
	l = 0
	do i = 0, 1
		do j = 0, 1
			do k = 0, 1
				vector = (/i,j,k/)
				l = l + 1
				vertices(l,:) = bbox(1,:) + (bbox(2,:) - bbox(1,:))*real(vector,8)
			end do
		end do
	end do
	
	return
  end subroutine compute_vertices_from_bbox
  
  recursive subroutine compute_all_vertices(node)
	implicit none
    type(node_type), intent(in) :: node
	integer i
	
	call compute_vertices_from_bbox(node%bbox, node%vertices)
    if (associated(node%children)) then
      do i = 1, size(node%children)
        call compute_all_vertices(node%children(i))
      end do
    end if

  end subroutine compute_all_vertices


  subroutine copy_node(node,copy)
  	implicit none
    type(node_type), intent(inout) :: node
    type(node_type) :: copy
    integer n
    
     
 	copy%num_point = node%num_point
    copy%depth = node%depth
    copy%node_id = node%node_id
    copy%child_num = node%child_num
    copy%bbox = node%bbox
    copy%centroid = node%centroid
    copy%nneighbor = node%nneighbor
	copy%leaf = node%leaf
	copy%length = node%length
	copy%cell = node%cell
  	
  end subroutine 

  subroutine associate_node(main_node,node)
  	implicit none
    type(node_type), intent(inout) :: main_node
    type(node_type) :: node
    integer n
    
    main_node = node

    if(allocated(node%point_ids)) then
    	if(.not. allocated(main_node%point_ids)) allocate(main_node%point_ids(node%num_point))
	    main_node%point_ids = node%point_ids
    end if
    
 	main_node%num_point = node%num_point
    main_node%depth = node%depth
    main_node%node_id = node%node_id
    main_node%child_num = node%child_num
    main_node%bbox = node%bbox
    main_node%vertices = node%vertices
    main_node%centroid = node%centroid
    main_node%nneighbor = node%nneighbor
	main_node%leaf = node%leaf
	main_node%length = node%length
	main_node%cell = node%cell
    
    
    if(associated(node%parent)) then
    	main_node%parent =>node%parent
    end if

    if(associated(node%children)) then
    	main_node%children =>node%children
    end if

    if(associated(node%neighbors)) then
		main_node%neighbors => node%neighbors
	end if
  	
  end subroutine 
	
   function check_if_they_share_vertex(listA, listB) result(flag)

	implicit none
	real(8) :: listA(8,3), listB(8,3)
	integer :: flag
	integer i, j, n

	flag = 0
	if(sum( abs( listA - listB ) ) == 0.0d0  ) then

		flag = -1
		return
	end if
	
	n = 8
	do i = 1, n
		do j = 1, n
			if( sum ( abs(listA(i,:) - listB(j,:))   ) == 0.0d0 ) then
				flag = 1
				return
			end if
		end do
	end do
 
	return
  end function check_if_they_share_vertex

function check_if_lie_on_edge(node1, node2) result(flag)

	implicit none
	logical :: flag
	type(node_type) :: node1, node2
	real(8) :: eps = 1.0d-10
	flag = .false.
	
	if( (abs(node1%centroid(1) - node2%centroid(1)) - (node1%length + node2%length)*0.5d0) <= eps .and. &
	&	(abs(node1%centroid(2) - node2%centroid(2)) - (node1%length + node2%length)*0.5d0) <= eps .and. &
	&	(abs(node1%centroid(3) - node2%centroid(3)) - (node1%length + node2%length)*0.5d0) <= eps ) flag = .true.

	return
end function check_if_lie_on_edge

function check_if_lie_on_plane(node1, origin, normal) result(flag)

	implicit none
	logical :: flag
	type(node_type) :: node1, node2
	real(8) :: origin(3), normal(3), projected(3)
	real(8) :: eps = 1.0d-10
	flag = .false.
	
	projected = normal*node1%centroid
	
	if(normal(1) /= 0.0d0) then
		if((abs(origin(1) - projected(1)) - node1%length*0.5d0) <= eps) flag = .true.	
	else if (normal(2) /= 0.0d0) then
		if((abs(origin(2) - projected(2)) - node1%length*0.5d0) <= eps) flag = .true.	
	else if (normal(3) /= 0.0d0) then	
		if((abs(origin(3) - projected(3)) - node1%length*0.5d0) <= eps) flag = .true.	
	else 
		print*, 'plane projection not defined'
		stop
	end if

	return
end function check_if_lie_on_plane


function check_if_on_boundary(point,bbox) result(ans)

	implicit none
	real(8) :: point(3)
	real(8) :: bbox(:,:)
	real(8) :: shift(3)!, ans(3)
	logical :: ans, flag
	
	ans = .false.
	flag = .false. 

	shift = abs(bbox(1,:))

!	bbox(2,:) = bbox(2,:) + shift	
!	point = point + shift
!	bbox(1,:) = bbox(1,:) + shift
	
	if(point(1) == bbox(1,1)) flag = .true. !point(1) = point(1) + 1.0d-14
	if(point(2) == bbox(1,2)) flag = .true. ! point(2) = point(2) + 1.0d-14
	if(point(3) == bbox(1,3)) flag = .true. ! point(3) = point(3) + 1.0d-14

	if(point(1) == bbox(2,1)) flag = .true. ! point(1) = point(1) - 1.0d-14
	if(point(2) == bbox(2,2)) flag = .true. ! point(2) = point(2) - 1.0d-14
	if(point(3) == bbox(2,3)) flag = .true. ! point(3) = point(3) - 1.0d-14

	ans = flag!point - shift

end function check_if_on_boundary


  subroutine compute_neighbors(node1, node2, ith)
  
  	implicit none
  	integer ith, i
    type(node_type), target :: node1, node2
  	
  	

	if(associated(node2%children)) then
		do i = 1, size(node2%children)	

			if( check_if_lie_on_edge(node1, node2%children(i)) .and.  &
			    node1%depth==node2%children(i)%depth .and.	&
!				node1%length==node2%children(i)%length .and.	&
				node1%node_id /= node2%children(i)%node_id) then


		
!			if( check_if_they_share_vertex(node1%vertices, node2%children(i)%vertices) == 1 .and. &
!			&   node1%depth==node2%children(i)%depth) then
				ith = ith + 1
!				call associate_node(node2%children(i),tmp_nodes(ith)%node)
				tmp_nodes(ith)%node => node2%children(i)
			end if 
		end do
	end if  	
  	
  	return

  end subroutine compute_neighbors
  
  
   recursive subroutine find_all_neighbors_recursive(node)
  
	implicit none
    type(node_type), intent(inout) :: node
  	integer i, j, ith, n
  	  	
	ith = 0
	allocate(tmp_nodes(26))

	if(associated(node%parent)) then
		call compute_neighbors(node, node%parent, ith)
		if(associated(node%parent%neighbors)) then
			n = node%parent%nneighbor
			if(n > 0) then
			do i = 1, n
				call compute_neighbors(node, node%parent%neighbors(i), ith)

!				if(.not. associated(node%children) .and. .not. associated(node%parent%neighbors(i)%children)) then
!				if(node%leaf .and. node%parent%neighbors(i)%leaf) then
!					if( check_if_they_share_vertex(node%vertices, node%parent%neighbors(i)%vertices) == 1 ) then
!						ith = ith + 1
!						tmp_nodes(ith)%node => node%parent%neighbors(i)
!					end if
!				end if
				
				!call find_neighbor(node, node%parent%neighbors(i), ith)
				
			end do
			end if

		end if

	end if

	if(ith > 26) then
		print*, '====================='
		print*, 'more neighbors than possible'
		print*, '====================='
	end if

	if(ith /= 0) then
	   	allocate(node%neighbors(ith))
		node%nneighbor = ith
		do i = 1, ith
			call associate_node(node%neighbors(i),tmp_nodes(i)%node)
	  	end do
	else
		nullify(node%neighbors)
	end if

  	deallocate(tmp_nodes)  

	if(associated(node%children)) then
		do i = 1, size(node%children)
		 call find_all_neighbors_recursive(node%children(i))
		end do
	end if	
  
  	return
  end subroutine find_all_neighbors_recursive
  

 subroutine find_all_neighbors
  
	implicit none
	type(node_type), pointer :: node, cnode
  	integer depth, icell, jcell  	  	
	integer tmpi(200)
	integer i, j, k , l	, n, nneighbors
	integer ith, id
	logical flag
	
	call get_all_leaves

	do depth =  1, tree%max_depth
		if(tree%depth_list(depth)%ncell /= 0 ) then
		do icell = 1, tree%depth_list(depth)%ncell


			cnode => tree%depth_list(depth)%nodes(icell)%node
			id = cnode%node_id
			allocate(tmp_nodes(1000))

			ith = 0
		
!			Direct Neighbours		
			do jcell = 1, tree%depth_list(depth)%ncell
				node => tree%depth_list(depth)%nodes(jcell)%node
				if( check_if_they_share_vertex(cnode%vertices, node%vertices) == 1 ) then
					ith = ith + 1
					tmp_nodes(ith)%node => node				
				end if				
				nullify(node)	
			end do


			if(cnode%leaf) then
			

			
!				If parent of leaf shares vertex			
				do jcell = 1, tree%nleaves
				node => tree%leaves(jcell)%node
				flag = .false.
				
!				if(cnode%node_id /= node%node_id .and. node%depth /= cnode%depth) then
				if(cnode%node_id /= node%node_id .and. abs(node%depth - cnode%depth) < 2) then	
							
				if( (check_if_they_share_vertex(cnode%vertices, node%parent%vertices) == 1) .or.	&
				&	(check_if_they_share_vertex(cnode%parent%vertices, node%vertices) == 1 ) ) then
				
				
					do i = 1, ith-1
						if(tmp_nodes(i)%node%node_id == node%node_id) then
							print*, 'cell already included in the neighbour list'
							go to 100
						end if
					end do
				
				
					ith = ith + 1
					tmp_nodes(ith)%node => node
					flag = .true.							
					
				end if			
				
							
				
				end if
100				continue				
				nullify(node)						
				end do



			end if


			if(ith /= 0) then
			   	allocate(cnode%neighbors(ith))
				cnode%nneighbor = ith
				do i = 1, ith
					call associate_node(cnode%neighbors(i),tmp_nodes(i)%node)
	  			end do
			else
				nullify(cnode%neighbors)
			end if

		  	deallocate(tmp_nodes)  


			nullify(cnode)

		end do
		end if
	end do	
		
  
  	return
  end subroutine find_all_neighbors
  

 


  recursive subroutine get_all_depths(node)

	implicit none
    type(node_type), target :: node
	integer i, depth, ncell
	
	
	depth = node%depth
	if(tree%max_depth < depth) tree%max_depth = depth
	if(.not. allocated(tree%depth_list(depth)%nodes)) then
		ncell = tree%depth_list(depth)%ncell
		if(ncell /= 0) then
			allocate(tree%depth_list(depth)%nodes(ncell))
			tree%depth_list(depth)%ncell = 0
		end if
	end if
	
	if(allocated(tree%depth_list(depth)%nodes)) then
		tree%depth_list(depth)%ncell = tree%depth_list(depth)%ncell + 1
		if(node%num_point>0) then
			tree%depth_list(depth)%nonempty = tree%depth_list(depth)%nonempty + 1
		end if
		ncell = tree%depth_list(depth)%ncell
		node%cell = ncell
		tree%depth_list(depth)%nodes(ncell)%node => node
	end if
	
	if(associated(node%children)) then
		do i = 1, size(node%children)
			call get_all_depths(node%children(i))
		end do
	end if

	return

  end subroutine get_all_depths

  subroutine get_all_leaves
  
  	implicit none
  	
  	integer :: i 
  	i = 1
  	if(.not. allocated(tree%leaves)) allocate(tree%leaves(tree%nleaves))
    call find_leaves(tree%root_node, tree%leaves, i)
  
  end subroutine get_all_leaves



  recursive subroutine find_leaves(node, leaves, ith)

    type(node_type), intent(in), target :: node
  	type(nodes_type) :: leaves(:)
  	integer ith
    integer i

    if (associated(node%children)) then
      do i = 1, size(node%children)
        call find_leaves(node%children(i), leaves, ith)
      end do
    else
!      if (node%num_point == 0) then
!       	 return
!      end if

  	leaves(ith)%node => node
  	leaves(ith)%node%leaf_id = ith
  	
  	ith = ith + 1     
    end if

  end subroutine find_leaves



  subroutine print_node(node)

    type(node_type), intent(in) :: node
    write(6, "('----------------------------------------------------------------')")
    write(6, "('Bounding box: ', 6F8.2)") node%bbox
    write(6, "('  Node id: ', I5)") node%node_id      
    write(6, "('Depth: ', I3)") node%depth
    write(6, "('Child num: ', I3)") node%child_num      
    write(6, "('Point number: ', I3)") node%num_point
    write(6, "('Leaf?: ', L1)") .not. associated(node%children)
    write(6, "('  Neighbors: ', I3)") node%nneighbor


  end subroutine print_node

  recursive subroutine print_tree(node)

    type(node_type), intent(in) :: node

    integer i

    if (associated(node%children)) then
      write(6, "('----------------------------------------------------------------')")
      write(6, "('Branch node: ')") 
      write(6, "('  Node id: ', I5)") node%node_id      
      write(6, "('  Child num: ', I3)") node%child_num      
      write(6, "('  Bounding box: ', 6F8.2)") node%bbox
      write(6, "('  Centroid: ', 3F8.2)") node%centroid
      write(6, "('  Depth: ', I3)") node%depth
      write(6, "('  Neighbors: ', I3)") node%nneighbor

      do i = 1, size(node%children)
        call print_tree(node%children(i))
      end do
    else
      if (node%num_point == 0) then
       	 return
      end if

      write(6, "('Leaf node: ')") 
      write(6, "('  Node id: ', I5)") node%node_id      
      write(6, "('  Child num: ', I3)") node%child_num      
      write(6, "('  Bounding box: ', 6F8.2)") node%bbox
      write(6, "('  Centroid: ', 3F8.2)") node%centroid
      write(6, "('  Depth: ', I3)") node%depth
      write(6, "('  Neighbors: ', I3)") node%nneighbor
      write(6, "('  Points:')", advance='no')
      write(6, *) (node%point_ids(i), i = 1, node%num_point)


    end if
    

  end subroutine print_tree





  subroutine write_octree
  
  	implicit none
  	
  	integer i, j, k, n, ncell, id
	type(node_type), pointer :: cnode
  	


  	
  	open(unit = 287, file = 'octree.dat', status = 'unknown')
	write(287,*) 'VARIABLES="X","Y","Z"'

  	open(unit = 288, file = 'octree_leaves.dat', status = 'unknown')
	write(288,*) 'VARIABLES="X","Y","Z"'
	

	do i = 0, tree%max_depth
			
	if(tree%depth_list(i)%ncell /= 0 ) then
		

		ncell = tree%depth_list(i)%ncell
							
			
  cells: do j = 1, ncell
				
			cnode => tree%depth_list(i)%nodes(j)%node
  			
  			if (cnode%num_point == 0) cycle cells
  			
			write(287,*) 'Zone'
			write(287,*) 'ZoneType=FEBrick'
			write(287,*) 'Nodes=8'
			write(287,*) 'Elements=1'
			write(287,*) 'datapacking=point'
  			
			do k = 1, 8
				write(287,*) cnode%vertices(k,:)
			end do
  			write(287,*) 1, 5, 7, 3, 2, 6, 8, 4
  			  			
  			if(cnode%leaf) then
 
	 			write(288,*) 'Zone'
				write(288,*) 'ZoneType=FEBrick'
				write(288,*) 'Nodes=8'
				write(288,*) 'Elements=1'
				write(288,*) 'datapacking=point'
  			
				do k = 1, 8
					write(288,*) cnode%vertices(k,:)
				end do
	  			write(288,*) 1, 5, 7, 3, 2, 6, 8, 4 			
  			
  			end if
  			
  		end do cells
  	end if
  	end do
  
  close(287)
  close(288)
  return
  end subroutine write_octree







end module octree_mod
