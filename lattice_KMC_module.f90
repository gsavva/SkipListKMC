module lattice_KMC
   use skipList_module

   implicit none
   integer                               :: Ns, c, empty, occupied
   integer,  allocatable, dimension(:)   :: reaction_effect
   integer,  allocatable, dimension(:,:) :: lattice, coords, neighbours
   real*8,   allocatable, dimension(:)   :: reaction_const
   real*8,   allocatable, dimension(:,:) :: propensities
   integer*8,allocatable, dimension(:)   :: h
   real*8                                :: ads_const  = 1., des_const= 1.
   real*8                                :: diff_const = 0.01

contains
   subroutine init(queue_struct,n,percentage, p_val)
      implicit none
      class(SkipList_Type) ::  queue_struct
      integer :: n!, i
      real*8  :: percentage
      real*8, optional :: p_val

      Ns = n*n ! total number of sites

      allocate(lattice(n,n))
      allocate(coords(Ns,2))
      allocate(neighbours(Ns,4))
      allocate(propensities(Ns,6))
      
      if (present(p_val)) call queue_struct % set_parameters(p_val)
      ! else, the default value for p will be used

      call queue_struct%initialize(6*Ns)
      call find_coords(n)
      call find_neighbours(n)
      call randomize_coverage(percentage)

      call insert_events_to_skiplist(queue_struct) ! passes necessary info to skiplist
!~       call queue_struct % print_list_values()
   end

   subroutine execute_reaction(queue_struct,reaction_occured,t_kmc)
      implicit none
      class(SkipList_Type) ::  queue_struct
      real*8     :: t_kmc
      integer    :: reaction_occured, site_or, site_des, r_type!,i,j

      !specify which site is affected and the reaction_type (diff,des,ads)
      site_or = int((reaction_occured-1)/6) + 1 ! ∈ [1, Ns]
      r_type  = mod( reaction_occured-1, 6) + 1 ! ∈ [1,  6]
!~       print*, "Reaction Occured, Site, Type, Label:", reaction_occured, site_or, r_type, label(site_or, r_type)
      ! find the destination site and update the state of the lattice
      ! update matrix of propensites that lattice_KMC "sees" AND do the checks
      select case(r_type)
         case(1:4)! diffusions
            site_des = neighbours(site_or,r_type) ! neighbour according to reaction_type
            lattice(coords(site_or, 1), coords(site_or, 2)) = 0 ! change state of ORIGIN site
            lattice(coords(site_des,1), coords(site_des,2)) = 1 ! change state of DESTINATION site
            !----------------ORIGIN SITE-----------------------------------------------------------
            propensities(site_or, 1:4) = 0.0_8     ! further diffusions are disabled
            propensities(site_or, 5  ) = 0.0_8     ! DEsorption is  disabled
            propensities(site_or, 6  ) = ads_const ! ADsorption is   enabled
            call enable_diffusions(site_or)    ! diffusion BACK is permitted by default but needs special treatment
            !----------------DESTINATION SITE------------------------------------------------------
            propensities(site_des, 5 ) = des_const ! DEsorption is  enabled
            propensities(site_des, 6 ) = 0.0_8     ! ADsorption is disabled
            call check_permitted_diffusions(site_des)
         case(5)! DEsorption
            site_des = site_or !! needed
            lattice(coords(site_or, 1), coords(site_or, 2)) = 0 ! change state of ORIGIN site
            propensities(site_or, 1:4     ) = 0.0_8     ! diffusions are disabled
            propensities(site_or, r_type  ) = 0.0_8     ! DEsorption is  disabled
            propensities(site_or, r_type+1) = ads_const ! ADsorption is   enabled
            call enable_diffusions(site_or)!diffusion of the neighbours to the empty site HAS TO BE enabled
         case(6)! ADsorption
            site_des = site_or !! needed
            lattice(coords(site_or, 1), coords(site_or, 2)) = 1
            call check_permitted_diffusions(site_or) ! OK
            propensities(site_or, r_type-1) = des_const ! DEsorption is  enabled
            propensities(site_or, r_type  ) = 0.0_8     ! ADsorption is disabled
      end select

      ! update queuing structure (tree or heap) according to NEW positions/propensities

      call update_skiplist(queue_struct, site_or, site_des, t_kmc)

   end

!=======================================================================
!----------------------- Internally Used Subroutines -------------------
!=======================================================================

   subroutine insert_events_to_skiplist(queue_struct)
      implicit none
      class(SkipList_Type) ::  queue_struct
      real*8, dimension(6) :: rand_nums
      real*8  :: t
      integer :: i,j
      
      do i=1, Ns
         call random_number(rand_nums)
         do j=1,6
            t = -LOG(rand_nums(j)) / propensities(i,j)
!~             if ( t > huge(0.0_8) ) cycle ! only non-INFINITY times are inserted in the Skip_List
            call queue_struct % insert(t, label(i,j) )
         enddo
      enddo
   end

   subroutine update_skiplist(queue_struct, site_or, site_des, t_kmc)
      implicit none
      class(SkipList_Type) ::  queue_struct
      real*8, dimension(-1:9,4) :: rand_nums
      real*8  :: t_kmc, t_new1, t_new2
      integer :: i,j, site_or, site_des, n_or, n_des
      integer, dimension(4) :: mir=(/2,1,4,3/)
      call random_number(rand_nums)

      if(site_or == site_des) then
         do i=1,4
            t_new1 = t_kmc - LOG(rand_nums(-1,i)) / propensities(site_or, i)
            call queue_struct % update(label(site_or ,i), t_new1)

            j=mir(i)
            n_or  = neighbours(site_or , i)
            t_new1 = t_kmc - LOG(rand_nums(2*j-1,j)) / propensities(n_or ,j)
            call queue_struct % update( label(n_or ,j), t_new1)

         enddo
         do i=5,6 ! update AD & DES propensities of origin and destination sites
            j = i-4 ! j=1,2
            t_new1 = t_kmc - LOG(rand_nums(9, 2*j-1)) / propensities(site_or, i) ! rands = (9,1) (9,3)
            call queue_struct % update(label(site_or ,i), t_new1)

         enddo
      else
         do i=1,4 ! update origin & destination DIFFUSION propensities
            t_new1 = t_kmc - LOG(rand_nums(-1,i)) / propensities(site_or, i)
            t_new2 = t_kmc - LOG(rand_nums( 0,i)) / propensities(site_des,i)
            ! update(current_value, label, updated_value)
            call queue_struct % update( label(site_or ,i), t_new1)
            call queue_struct % update( label(site_des,i), t_new2)

            j=mir(i)
   !~          do j=1,4 ! update NEIGHBOURS of origin & destination sites
            n_or  = neighbours(site_or , i)
            n_des = neighbours(site_des, i)                     
            t_new1 = t_kmc - LOG(rand_nums(2*j-1,j)) / propensities(n_or ,j)
            t_new2 = t_kmc - LOG(rand_nums(2*j  ,j)) / propensities(n_des,j) 
            call queue_struct % update( label(n_or ,j), t_new1)
            call queue_struct % update( label(n_des,j), t_new2)

   !~          enddo
         enddo
         do i=5,6 ! update AD & DES propensities of origin and destination sites
            j = i-4 ! j=1,2
            t_new1 = t_kmc - LOG(rand_nums(9, 2*j-1)) / propensities(site_or, i) ! rands = (9,1) (9,3)
            t_new2 = t_kmc - LOG(rand_nums(9, 2*j  )) / propensities(site_des,i) ! rands = (9,2) (9,4)
   
            call queue_struct % update( label(site_or ,i), t_new1)
            call queue_struct % update( label(site_des,i), t_new2)

         enddo
      endif
   end

!------------------Subroutines to handle lattice-private data structures

   subroutine check_permitted_diffusions(site_of_interest)
      implicit none
      integer :: i, site_of_interest, neighbour, row, col
      integer, dimension(4) :: mir=(/2,1,4,3/)

      do i=1,4 !check the 4 neighbours for permitted diffusions
         neighbour = neighbours(site_of_interest, i) ! destination site
         row = coords(neighbour,1)
         col = coords(neighbour,2)
!~          print*,neighbour,lattice(row,col),mir(i),propensities(site_of_interest,i),propensities(neighbour,mir(i))
!~          if(lattice( coords(neighbours(site_of_interest,i),1),coords(neighbours(site_of_interest,i),2) ) == 0) then
         if (lattice(row, col)==0) then
            propensities(site_of_interest,i) = diff_const
         else ! if lattice site is already occupied
            propensities(site_of_interest,i) = 0.0_8
            ! need to disable mirrored reaction, 1-2, 2-1, 3-4, 4-3
            propensities(neighbour,mir(i))=0.0_8
         endif
!~          print*,neighbour,lattice(row,col),mir(i),propensities(site_of_interest,i),propensities(neighbour,mir(i))
      enddo
   end

   subroutine enable_diffusions(site_of_interest)
      implicit none
      integer :: i, site_of_interest, neighbour, row, col
      integer, dimension(4) :: mir=(/2,1,4,3/)
      do i=1,4
         neighbour = neighbours(site_of_interest, i) ! destination site
         row = coords(neighbour,1)
         col = coords(neighbour,2)
         if (lattice(row, col)==1) then
            propensities(neighbour,mir(i)) = diff_const
         endif
      enddo
   end

   subroutine find_neighbours(n)
      implicit none
      integer :: n, i
      ! find all NORTH--------------------------------------------------
      do i=1, Ns
         if (mod(i,n)==1) then
            neighbours(i,1) = i + n - 1
         else
            neighbours(i,1) = i - 1
         endif
      enddo
      ! find all SOUTH--------------------------------------------------
      do i=1, Ns
         if (mod(i,n)==0) then
            neighbours(i,2) = i - n + 1
         else
            neighbours(i,2) = i + 1
         endif
      enddo
      ! find all EAST---------------------------------------------------
      do i=1, Ns-n
         neighbours(i,3) = i + n
      enddo
      do i=1, n
         neighbours(Ns-n+i,3) = i
      enddo
      ! find all WEST---------------------------------------------------
      do i=n+1, Ns
         neighbours(i,4) = i - n
      enddo
      do i=1, n
         neighbours(i,4) = Ns - n + i
      enddo
   end

   subroutine find_coords(n)
      implicit none
      integer :: i,n
      do i=1, Ns ! fill coords array
         coords(i,1) = mod( i-1,n)  + 1 ! row of i-lattice site
         coords(i,2) = int((i-1)/n) + 1 ! col of i-lattice site
      enddo      
   end
   
   integer function label(l_site, r_type)
      implicit none
      integer, intent(in) :: l_site, r_type
      label = (l_site-1)*6 + r_type
   end

   subroutine randomize_coverage(covrg)
      implicit none
      integer :: i,j, site, c
      real*8  :: randN, covrg
      integer, dimension(Ns) :: rand_positions
      
      ! initialize an empty lattice, only adsorption can occur
      lattice = 0
      propensities(:, 1:5) = 0.0_8     ! diffusions + DEsorption
      propensities(:, 6  ) = ads_const ! ADsorption
      if ( covrg > 0. ) then
         do i=1, Ns ! fill shuffled array
            call random_number(randN)
            j = 1 + FLOOR(i*randN)
            if(j /= i) rand_positions(i) = rand_positions(j)
            rand_positions(j) = i
         enddo
         c = NINT(covrg*Ns) ! number of randomly selected iccupied 
         print*,"Initially Occupied: ", c, " out of", Ns
         do i=1, c ! insert c sites in lattice
            site = rand_positions(i)
!~             print*,"Site Occupied: ", site
            lattice(coords(site,1), coords(site,2)) = 1
            call check_permitted_diffusions(site)
            propensities(site, 5) = des_const ! DEsorption is  enabled
            propensities(site, 6) = 0.0_8     ! ADsorption is disabled
         enddo
      endif
   end

end module lattice_KMC
