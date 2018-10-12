module node_module
   implicit none

   type :: container
      type(node), pointer :: p => null() ! the pointer "p" can point only to "node" objects
   end type

   type :: node
      real*8   :: node_val
      integer  :: node_id
      type(container), allocatable :: forward_nodes(:) ! array of pointers to forward nodes      
      type(container), allocatable :: backward_nodes(:)      
   contains
      procedure :: getValue
      procedure :: getID
      procedure :: getLvl
   end type
   
   interface node
      procedure create_new_node
   end interface

contains

   real*8 function getValue(this)
      implicit none
      class(node) :: this
      getValue = this % node_val
   end function

   real*8 function getID(this)
      implicit none
      class(node) :: this
      getid = this % node_id
   end function

   real*8 function getLvl(this)
      implicit none
      class(node) :: this
      getLvl = size(this % forward_nodes)
   end function

   function create_new_node(new_level, new_value, new_id)
      implicit none
      type(node), pointer :: create_new_node 
      integer     :: new_level, new_id
      real*8      :: new_value

      allocate( create_new_node )
      allocate( create_new_node %  forward_nodes(new_level) )
      allocate( create_new_node % backward_nodes(new_level) )
      create_new_node % node_val = new_value
      create_new_node % node_id  = new_id
   end

end module
