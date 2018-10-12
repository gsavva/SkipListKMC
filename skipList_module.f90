module skipList_module
   use node_module
   implicit none

   real(8), parameter :: d_QpInf = transfer(z'7FF0000000000000',1.0_8)

!~    type, extends(ProcessQueue_Type) :: SkipList_Type
   type :: SkipList_Type
      type(node), pointer :: head => null()
      type(node), pointer :: tail => null()
      type(container), allocatable :: updt(:) ! array of pointers, pointing to nodes
      type(container), allocatable :: mapArr(:)

      integer :: nsize
      integer :: maxlevels
      integer :: current_list_level
      real*8  :: p_val
   contains
      procedure :: initialize => skiplist_initialize
      procedure :: insert => skiplist_insert_element
      procedure :: remove => skiplist_remove_element
      procedure :: search => skiplist_search_element
      procedure :: update => skiplist_update_element
      procedure :: get_first => skiplist_get_highest_priority_element
      
      procedure :: printAll
      procedure :: print_list_values
      procedure :: deleteList
      procedure :: choose_random_height
   end type

contains

   subroutine skiplist_initialize(this, n0, prop_val)
      implicit none
      class(SkipList_Type) :: this
      real *8 :: prop_val
      integer :: n0, i
      this % p_val = prop_val
      this % nsize = 0
      this % current_list_level = 1
      this % maxlevels = int( log(n0*1.) / log(1./this%p_val) ) + 1
                    !node(level, value, ID)
      this % head => node(this % maxlevels,   0.0_8, -11) ! head should have a zero value
      this % tail => node(this % maxlevels, d_QpInf, -99) ! tail should have +infinity value
      allocate(this % updt(this % maxlevels))
      allocate(this % mapArr(n0))

      do i=1, this % maxlevels ! settle initial connectivity: head --> tail --> null
         this % head % forward_nodes(i)%p => this % tail  ! head is pointing to tail
         this % tail % forward_nodes(i)%p => null()       ! tail is pointing NOwhere
         this % head %backward_nodes(i)%p => null()
         this % tail %backward_nodes(i)%p => this % head
         this % updt(i)%p => null()
      enddo
   end

   subroutine skiplist_insert_element(this, element_to_insert, label)
      implicit none
      class(SkipList_Type) :: this
      type (node), pointer :: new_node, x
      real*8  :: element_to_insert
      integer :: height, i, label

      this % nsize = this % nsize + 1

      x => this % head
      do i = this % maxlevels, 1, -1
         do while(x % forward_nodes(i)%p % node_val < element_to_insert)
            x =>  x % forward_nodes(i)%p
         enddo
         this % updt(i)%p => x
      enddo

      height = this % choose_random_height()
      if (height > this % current_list_level) then
         height = this % current_list_level + 1  ! level is increased at max by 1 per new node 
         this % current_list_level = this % current_list_level + 1
      endif

      new_node => node(height, element_to_insert, label) ! create the new node
      do i=1, height ! iterate new_node's levels
         new_node %  forward_nodes(i) = this % updt(i)%p % forward_nodes(i) ! SHOULD NOT USE  " => "
         new_node % backward_nodes(i) = this % updt(i)

         this % updt(i)%p % forward_nodes(i)%p => new_node ! this assignment will also change the above if => is used instead of =
         new_node % forward_nodes(i)%p % backward_nodes(i)%p => new_node
      enddo
      this % mapArr(label)%p => new_node
   end

   subroutine skiplist_remove_element(this, label_to_remove)
      implicit none
      class(SkipList_Type) :: this
      type(node), pointer :: x
      integer  :: label_to_remove, i
      
      this % nsize = this % nsize - 1

      x => this % mapArr(label_to_remove)%p ! direct access to the element of interest
      do i = 1, size( x % forward_nodes )
         x % backward_nodes(i)%p %  forward_nodes(i)%p => x %  forward_nodes(i)%p
         x %  forward_nodes(i)%p % backward_nodes(i)%p => x % backward_nodes(i)%p
      enddo
      deallocate(x)
      this % mapArr(label_to_remove)%p => null()

   end

   subroutine skiplist_search_element(this, element_to_search)
      implicit none
      class(SkipList_Type) :: this
      type(node), pointer :: x
      real*8   :: element_to_search
      integer  :: i

      x => this % head
      do i = this % maxlevels, 1, -1
         do while(x % forward_nodes(i)%p % node_val < element_to_search)
            x =>  x % forward_nodes(i)%p
         enddo
      enddo
      x => x % forward_nodes(1)%p
      if ( (x%node_val - element_to_search)<1e-6) then
         write(*,*) "number: ", element_to_search, " found. ID: ", x%node_id
      else
         write(*,*) "number: ", element_to_search, " NOT found"
      endif
   end

   subroutine skiplist_update_element(this, label, updated_value)
      implicit none
      class(SkipList_Type) :: this
      real*8   :: updated_value
      integer  :: label
      
      if( updated_value > huge(0.0_8) .AND. .NOT. ASSOCIATED(this % mapArr(label)%p) ) return

      if( updated_value > huge(0.0_8) .AND.       ASSOCIATED(this % mapArr(label)%p) ) then
         call this % remove(label)
         return
      endif
      ! if updated_value < +infinity
      if ( .NOT. ASSOCIATED(this % mapArr(label)%p) ) then
         call this % insert(updated_value, label)
      else ! if mapArr(label)%p IS associated
         call this % remove(label)
         call this % insert(updated_value, label)
      endif

   end

   subroutine skiplist_get_highest_priority_element(this,event_id, event_t)
      implicit none
      class(skipList_Type) :: this
      integer :: event_id
      real*8  :: event_t

      event_id = this % head % forward_nodes(1)%p % node_id
      event_t  = this % head % forward_nodes(1)%p % node_val
   end

   subroutine printAll(this)
      implicit none
      class(SkipList_Type) :: this
      integer :: i, j
      type (node), pointer :: run

      run => this%head
      do i=1, this%nsize+2
         print*, i
         print*,"node value", run%node_val
         do j = size(run%forward_nodes), 1, -1
            if(associated(run % forward_nodes(j) % p)) then
               print*,"Node ",i, "level ",j," points to ", run % forward_nodes(j) % p % node_val, "FRONT"
            endif
            if(associated(run %backward_nodes(j) % p)) then
               print*,"Node ",i, "level ",j," points to ", run %backward_nodes(j) % p % node_val, "BACK"
            endif
         enddo
         print*, char(10), "Next Node:"
         run => run % forward_nodes(1) % p
      enddo
      print*, "Mapping Array data: (label, value)"
      do i=1, size(this%mapArr)
         if ( associated(this % mapArr(i)%p) ) then
            print*, this % mapArr(i)%p%node_id, this % mapArr(i)%p%node_val
         endif
      enddo
   end

   subroutine print_list_values(this)
      implicit none
      class(SkipList_Type) :: this
      type(node), pointer :: x
      
      x => this % head
      print*,"==============================LIST VALUES START====================================="
      print*, "Size of list: ", this % nsize
      print*, "               Value               ID       Level"
      do while ( associated(x%forward_nodes(1)%p) )
         print*, x % node_val, x % node_id , size(x%forward_nodes)
         x => x % forward_nodes(1)%p
      enddo
      print*, x % node_val, x % node_id, size(x%forward_nodes) ! prints the +inf, the last value on the list.
      print*,"==============================LIST VALUES END======================================="
   end

   subroutine deleteList(this)
      implicit none
      class(SkipList_Type) :: this
      type (node), pointer :: current, next

      current => this % head
      do while (associated(current))
         next => current % forward_nodes(1)%p
         deallocate(current)
         current => next
      enddo
      deallocate(this % updt)
      deallocate(this % mapArr)
   end

!--------------INTERNAL procedures -------------------------------------
   integer function choose_random_height(this)
      implicit none
      class(SkipList_Type) :: this
      real*8 :: temp

      choose_random_height = 1
      call random_number(temp)
      do while (temp < this % p_val .AND. choose_random_height < this % maxlevels)
         choose_random_height = choose_random_height + 1
         call random_number(temp)
      enddo
      ! may result to a faster alternative
!~       call random_number(temp)
!~       do while (temp < this%p_val)
!~          choose_random_height = choose_random_height + 1
!~          call random_number(temp)
!~       enddo
!~       choose_random_height = min(choose_random_height, this%maxlevels)
   end

end module
