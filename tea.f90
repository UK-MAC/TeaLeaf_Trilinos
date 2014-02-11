!Crown Copyright 2014 AWE.
!
! This file is part of TeaLeaf.
!
! TeaLeaf is free software: you can redistribute it and/or modify it under 
! the terms of the GNU General Public License as published by the 
! Free Software Foundation, either version 3 of the License, or (at your option) 
! any later version.
!
! TeaLeaf is distributed in the hope that it will be useful, but 
! WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
! FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more 
! details.
!
! You should have received a copy of the GNU General Public License along with 
! TeaLeaf. If not, see http://www.gnu.org/licenses/.

!>  @brief Driver for the heat conduction kernel
!>  @author David Beckingsale, Wayne Gaudin
!>  @details Invokes the user specified kernel for the heat conduction

MODULE tea_leaf_module

CONTAINS

SUBROUTINE tea_leaf()
 
  USE clover_module
  USE tea_leaf_kernel_module
  USE update_halo_module

  IMPLICIT NONE

!$ INTEGER :: OMP_GET_THREAD_NUM
  INTEGER :: c, n, j,k
  REAL(KIND=8) :: ry,rx, error

  INTEGER :: fields(NUM_FIELDS)

  REAL(KIND=8) :: kernel_time,timer

  DO c=1,number_of_chunks

    IF(chunks(c)%task.EQ.parallel%task) THEN

      fields=0
      fields(FIELD_ENERGY1) = 1
      fields(FIELD_DENSITY1) = 1
      CALL update_halo(fields,2)

      ! INIT
      IF(profiler_on) kernel_time=timer()
      IF(use_fortran_kernels) THEN
          CALL tea_leaf_kernel_init(chunks(c)%field%x_min, &
              chunks(c)%field%x_max,                       &
              chunks(c)%field%y_min,                       &
              chunks(c)%field%y_max,                       &
              chunks(c)%field%celldx,                      &
              chunks(c)%field%celldy,                      &
              chunks(c)%field%volume,                      &
              chunks(c)%field%density1,                    &
              chunks(c)%field%energy1,                     &
              chunks(c)%field%work_array1,                 &
              chunks(c)%field%u,                           &
              chunks(c)%field%work_array2,                 &
              chunks(c)%field%work_array3,                 &
              chunks(c)%field%work_array4,                 &
              chunks(c)%field%work_array5,                 &
              chunks(c)%field%work_array6,                 &
              chunks(c)%field%work_array7,                 &
              coefficient)
      ELSEIF(use_C_kernels) THEN
          CALL tea_leaf_kernel_init_c(chunks(c)%field%x_min, &
              chunks(c)%field%x_max,                       &
              chunks(c)%field%y_min,                       &
              chunks(c)%field%y_max,                       &
              chunks(c)%field%celldx,                      &
              chunks(c)%field%celldy,                      &
              chunks(c)%field%volume,                      &
              chunks(c)%field%density1,                    &
              chunks(c)%field%energy1,                     &
              chunks(c)%field%work_array1,                 &
              chunks(c)%field%u,                           &
              chunks(c)%field%work_array2,                 &
              chunks(c)%field%work_array3,                 &
              chunks(c)%field%work_array4,                 &
              chunks(c)%field%work_array5,                 &
              chunks(c)%field%work_array6,                 &
              chunks(c)%field%work_array7,                 &
              coefficient)
      ENDIF


      ! JACOBI

      rx = dt/(chunks(c)%field%celldx(chunks(c)%field%x_min)**2);
      ry = dt/(chunks(c)%field%celldy(chunks(c)%field%y_min)**2);

      IF(.NOT. use_trilinos_kernels) THEN

      DO n=1,max_iters

        IF(use_fortran_kernels) THEN
            CALL tea_leaf_kernel_solve(chunks(c)%field%x_min,&
                chunks(c)%field%x_max,                       &
                chunks(c)%field%y_min,                       &
                chunks(c)%field%y_max,                       &
                rx,                                          &
                ry,                                          &
                chunks(c)%field%work_array6,                 &
                chunks(c)%field%work_array7,                 &
                error,                                       &
                chunks(c)%field%work_array1,                 &
                chunks(c)%field%u,                           &
                chunks(c)%field%work_array2)
        ELSEIF(use_C_kernels) THEN
            CALL tea_leaf_kernel_solve_c(chunks(c)%field%x_min,&
                chunks(c)%field%x_max,                       &
                chunks(c)%field%y_min,                       &
                chunks(c)%field%y_max,                       &
                rx,                                          &
                ry,                                          &
                chunks(c)%field%work_array6,                 &
                chunks(c)%field%work_array7,                 &
                error,                                       &
                chunks(c)%field%work_array1,                 &
                chunks(c)%field%u,                           &
                chunks(c)%field%work_array2)
        ENDIF

        ! CALL update_halo
        fields=0
        fields(FIELD_U) = 1
        CALL update_halo(fields,2)

        CALL clover_max(error)

        IF (error .LT. eps) EXIT

      ENDDO
    ELSEIF(use_trilinos_kernels) THEN
        CALL trilinos_solve(                                   &
              grid%x_cells, &
              grid%y_cells, &
              chunks(c)%field%left,                            &
              chunks(c)%field%right,                           &
              chunks(c)%field%bottom,                          &
              chunks(c)%field%top,                             &
              1,                                               &
              grid%x_cells,                                    &
              1,                                               &
              grid%y_cells,                                    &
              rx,                                              &
              ry,                                              &
              chunks(c)%field%work_array6,                     &
              chunks(c)%field%work_array7,                     &
              chunks(c)%field%u)
      ENDIF

      IF (parallel%boss) THEN
!$      IF(OMP_GET_THREAD_NUM().EQ.0) THEN
          WRITE(g_out,"('Conduction error ',e14.7)") error
          WRITE(g_out,"('Iteration count ',i8)") n-1
          WRITE(0,"('Conduction error ',e14.7)") error
          WRITE(0,"('Iteration count ', i8)") n-1
!$      ENDIF
      ENDIF

      ! RESET
      IF(use_fortran_kernels) THEN
          CALL tea_leaf_kernel_finalise(chunks(c)%field%x_min, &
              chunks(c)%field%x_max,                           &
              chunks(c)%field%y_min,                           &
              chunks(c)%field%y_max,                           &
              chunks(c)%field%energy1,                         &
              chunks(c)%field%density1,                        &
              chunks(c)%field%u)
      ELSEIF(use_C_kernels) THEN
          CALL tea_leaf_kernel_finalise_c(chunks(c)%field%x_min, &
              chunks(c)%field%x_max,                           &
              chunks(c)%field%y_min,                           &
              chunks(c)%field%y_max,                           &
              chunks(c)%field%energy1,                         &
              chunks(c)%field%density1,                        &
              chunks(c)%field%u)
      ENDIF

      fields=0
      fields(FIELD_ENERGY1) = 1
      CALL update_halo(fields,1)

    ENDIF

  ENDDO
  IF(profiler_on) profiler%PdV=profiler%tea+(timer()-kernel_time)

END SUBROUTINE tea_leaf

END MODULE tea_leaf_module
