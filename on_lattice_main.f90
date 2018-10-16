program on_lattice_MAIN
!~    use node_module
   use skipList_module
   use lattice_KMC

   implicit none
   class(SkipList_Type), allocatable ::  queue_struct
   character(len=20) :: write_file
   real*8    :: t_kmc, t1, t2, cov, ts, tf, p
   integer   :: i, j, iters, reaction_occured, Ldim, seed(33)=101
!~    call random_seed(PUT=seed)
!~    call cpu_time(ts) !-----------------global start time----------------

   write_file= 'lattice.txt'
!~    write(*,*) 'Give lattice Dimension, # of iterations:' !, coverage:'
   read*, Ldim, iters, p!, cov
!~    Ldim = 5
!~    iters = 1000
   cov = 0.0
   p = p/100. ! used to enable batch script run since it does not support decimals
   allocate(queue_struct)
   call init(queue_struct,Ldim,cov, p)

   t_kmc = 0
   reaction_occured=-1

!~    print*,"===============START========================================="
   open(unit=88,file=write_file)
   do j=1, Ldim
      write(88,*) lattice(j,:)
   enddo
!~    open(unit=55,file='statistics')
   call cpu_time(ts)
   do i=1, iters
!~       print*,"**************************************************** NEW ITERATION ABOUT TO BEGIN &
!~             & ****************************************************"
      call queue_struct % get_first(reaction_occured, t_kmc)

!~       call queue_struct % print_list_values()
!~       write(55,*) reaction_occured, t_kmc, queue_struct % nsize
!~       call queue_struct % print_list_values()
!~       if (mod(i,1000) == 0) then
!~           print*,i, reaction_occured, t_kmc
!~       endif

!~       call cpu_time(t1)
      call execute_reaction(queue_struct, reaction_occured, t_kmc)
!~       call cpu_time(t2)

!~       print*,"time to EXEC next reaction: ",t2-t1
!~       print*,"a0 AFTER reaction execution:",reaction_occured, queue_struct%tree_elements(sums_tree%head_node_indx)

!~        if (mod(i,5000)==0) then
!~          do j=1, Ldim
!~             write(88,*) lattice(j,:)
!~          enddo
!~        endif

   enddo
   print*, queue_struct % st(:)
   call cpu_time(tf) !-----------------global finish time---------------
   open(unit=77,file='stats',position='append')
   write(77,*) Ldim*Ldim, iters, tf-ts
   print*,"Total time: ",tf-ts, " seconds"

end program
