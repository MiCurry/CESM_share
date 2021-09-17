MODULE shr_mem_mod

  use shr_kind_mod, only : shr_kind_r8
  use shr_log_mod, only: s_logunit => shr_log_Unit
  use shr_sys_mod, only: shr_sys_abort

  implicit none
  private

  ! PUBLIC: Public interfaces

  public ::  shr_mem_getusage, &
       shr_mem_init
  public :: record_memusage

  ! PUBLIC: Public interfaces

  real(shr_kind_r8) :: mb_blk = 0.0_shr_kind_r8

  !===============================================================================
CONTAINS
  !===============================================================================

  subroutine shr_mem_init(prt, strbuf)

    implicit none

    !----- arguments -----

    logical, optional :: prt
    character(len=*), optional :: strbuf
    !----- local -----

    ! --- Memory stats ---
    integer :: msize                   ! memory size (high water)
    integer :: mrss0,mrss1,mrss2       ! temporary rss
    integer :: mshare,mtext,mdatastack
    logical :: lprt
    integer :: ierr

    integer :: GPTLget_memusage

    real(shr_kind_r8),allocatable :: mem_tmp(:)

    character(*),parameter :: subname = "(shr_mem_init)"
    !---------------------------------------------------

    lprt = .false.
    if (present(prt)) then
       lprt = prt
    endif

    ierr = GPTLget_memusage (msize, mrss0, mshare, mtext, mdatastack)
    if (ierr .ne. 0) call shr_sys_abort(trim(subname)//': GPTLget_memusage mrss0 failed')

    allocate(mem_tmp(1024*1024), stat=ierr)    ! 1 MWord, 8 MB
    if (ierr .ne. 0) call shr_sys_abort(trim(subname)//': allocate failed')

    mem_tmp = -1.0
    ierr = GPTLget_memusage (msize, mrss1, mshare, mtext, mdatastack)
    if (ierr .ne. 0) call shr_sys_abort(trim(subname)//': GPTLget_memusage mrss1 failed')

    deallocate(mem_tmp, stat=ierr)
    if (ierr .ne. 0) call shr_sys_abort(trim(subname)//': deallocate failed')

    ierr = GPTLget_memusage (msize, mrss2, mshare, mtext, mdatastack)
    if (ierr .ne. 0) call shr_sys_abort(trim(subname)//': GPTLget_memusage mrss2 failed')

    mb_blk = 0.0_shr_kind_r8
    if (mrss1 - mrss0 > 0) then
       mb_blk = (8.0_shr_kind_r8)/((mrss1-mrss0)*1.0_shr_kind_r8)
    endif

    if (lprt) then
       write(s_logunit,'(A,f16.2)') '8 MB memory   alloc in MB is ',(mrss1-mrss0)*mb_blk
       write(s_logunit,'(A,f16.2)') '8 MB memory dealloc in MB is ',(mrss1-mrss2)*mb_blk
       write(s_logunit,'(A,f16.2)') 'Memory block size conversion in bytes is ',mb_blk*1024.0_shr_kind_r8*1024.0_shr_kind_r8
    endif
    if (present(strbuf)) then
       write(strbuf,'(3(A,f16.2))') '8 MB memory   alloc in MB is ',(mrss1-mrss0)*mb_blk, &
            '\n8 MB memory dealloc in MB is ',(mrss1-mrss2)*mb_blk, &
            '\nMemory block size conversion in bytes is ',mb_blk*1024.0_shr_kind_r8*1024.0_shr_kind_r8
    endif


  end subroutine shr_mem_init

  !===============================================================================

  subroutine get_mpi_info(mpicom, masterprocid, masterproc, npes)

     use mpi, only: mpi_comm_rank, mpi_comm_size

     integer, intent(in) :: mpicom
     integer, intent(out) :: masterprocid
     logical, intent(out) :: masterproc
     integer, intent(out) :: npes
     integer :: iam
     integer :: ierr

     call mpi_comm_rank(mpicom, iam, ierr)
     call mpi_comm_size(mpicom, npes, ierr)

     if (iam == 0) then
         masterproc = .true.
     else
         masterproc = .false.
     endif

     masterprocid = 0

  end subroutine get_mpi_info

  ! Get an IO unit that is currently not opened for writing/reading
  subroutine get_new_unit(new_unit)

      implicit none

      integer, intent(out) :: new_unit
      integer, parameter :: max_units = 99
      integer :: i
      logical :: opened

      do i = 1, max_units
          inquire(i, opened=opened)
          if (opened) then
              new_unit = i
              return
          end if
      end do

  end subroutine get_new_unit

  ! Given a filename, create a file, fname, and gather memory statistics from shr_mem_getusage
  ! across all processors and have the master proc record all them as a CSV file.
  subroutine record_memusage(mpicom, fname)

      use mpi, only : MPI_Gather, MPI_REAL8

      implicit none

      integer, intent(in) :: mpicom
      character(len=*), intent(in) :: fname

      real(shr_kind_r8) :: mem_hw
      real(shr_kind_r8) :: mem
      real(shr_kind_r8), dimension(:), allocatable :: arry_mem_hw
      real(shr_kind_r8), dimension(:), allocatable :: arry_mem
      logical :: masterproc
      integer :: masterprocid
      integer :: npes
      integer :: i, ierr
      integer :: output_unit

      call get_mpi_info(mpicom, masterprocid, masterproc, npes)

      call shr_mem_getusage(mem_hw, mem)

      if (masterproc) then
          call get_new_unit(output_unit)
          open(file=fname, unit=output_unit)
          allocate(arry_mem_hw(npes))
          allocate(arry_mem(npes))
      end if

      if (masterproc) then ! Writer CSV headers
          ! Write CSV headers ..
          write(output_unit,*) "ncpu,", "mem_hw,","mem,"
      end if

      call MPI_Gather(mem_hw, 1, MPI_REAL8, arry_mem_hw, 1, MPI_REAL8, masterprocid, mpicom, ierr)
      call MPI_Gather(mem, 1, MPI_REAL8, arry_mem, 1, MPI_REAL8, masterprocid, mpicom, ierr)

      if (masterproc) then
          ! Write output..
          do i = 1, npes
              write(output_unit, '(I10,A,F15.1,A,F15.1,A)') i, ',', arry_mem_hw(i), ',', arry_mem(i), ','
              if (mod(i, 36) == 0) then
                  write(output_unit, *) ", , ,"
              end if
          end do

          close(unit=output_unit)
          deallocate(arry_mem_hw)
          deallocate(arry_mem)
      end if

  end subroutine record_memusage

  subroutine shr_mem_getusage(r_msize,r_mrss,prt)

    implicit none

    !----- arguments ---
    real(shr_kind_r8) :: r_msize,r_mrss
    logical, optional :: prt

    !----- local ---
    integer :: msize,mrss
    integer :: mshare,mtext,mdatastack
    integer :: ierr
    integer :: GPTLget_memusage, GPTLprint_memusage

    !---------------------------------------------------

    ierr = GPTLget_memusage (msize, mrss, mshare, mtext, mdatastack)
    r_msize = msize / 1024.0_shr_kind_r8
    r_mrss  = mrss  / 1024.0_shr_kind_r8

    if (present(prt)) then
       if (prt) then
          ierr = GPTLprint_memusage(' ')
       endif
    endif


  end subroutine shr_mem_getusage

  !===============================================================================

END MODULE shr_mem_mod
