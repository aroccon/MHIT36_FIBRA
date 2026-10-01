module mod_linspace
    use iso_fortran_env, only: dp => real64
    implicit none

contains

    !>  
    !   Return evenly spaced numbers over a specified interval.
    !
    !   Returns `num` evenly spaced samples, calculated over the interval `[start, stop]`. 
    !
    !   Ported from the numpy routine.
    !
    !   Author: Ivan Pribec
    !
    function linspace(start,end,num,endpoint,step) result(samples)
        
        ! PARAMETERS
        real(dp), intent(in) :: start 
            !! The starting value of the sequence.
        real(dp), intent(in) :: end
            !! The end value of the sequence, unless `endpoint` is set to `.false.`. 
            !! In that case, the sequence consists of all but the last of `num + 1` 
            !! evenly spaced samples, so that `end` is excluded. Note that the 
            !! step size changes when `endpoint` is `.false.`.
        integer, intent(in), optional :: num
            !! Number of samples to generate. Default value is 50.
        logical, intent(in), optional :: endpoint
            !! If `.true.`, `end` is the last sample. Otherwise, it is not included. Default is `.true.`.
        real(dp), intent(out), optional :: step
            !! If present, `step` is the size of spacing between samples.

        ! RETURNS
        real(dp), allocatable :: samples(:)
            !! There are `num` equally spaced samples in the closed interval `[start, stop]` or 
            !! the half-open interval `[start, stop)` (depending on whether `endpoint` is `.true.` or `.false.`).

        integer :: num_, i
        logical :: endpoint_
        real(dp) :: step_

        num_ = 50
        if (present(num)) num_ = num

        endpoint_ = .true.
        if (present(endpoint)) endpoint_ = endpoint

        ! find step size
        if (endpoint_) then
            step_ = (end - start)/real(num_-1,dp)
        else
            step_ = (end - start)/real(num_,dp)
        end if

        if (present(step)) step = step_

        allocate(samples(num_))
        do i = 1, num_
            samples(i) = start + (i-1)*step_
        end do
    end function linspace
  

    !---     Quadrature rule     ----!
    subroutine quad_rule(noGPs, x, w)
    IMPLICIT NONE
    
    INTEGER, PARAMETER      :: rk = kind ( 1.0D+00 )
    INTEGER, INTENT(IN)     :: noGPs
    REAL(16), INTENT(INOUT) :: x(noGPs), w(noGPs)

    if (noGPs .EQ. 1) then
       x(1) = 0.0D+00
    else if (noGPs .EQ. 2) then
       x(1) = - 1.0D+00 / sqrt ( 3.0D+00 )
       x(2) =   1.0D+00 / sqrt ( 3.0D+00 )
    else if (noGPs .EQ. 3) then
       x(1) = - 1.0D+00 / sqrt ( 2.0D+00 )
       x(2) =   0.0D+00
       x(3) =   1.0D+00 / sqrt ( 2.0D+00 )
    else if (noGPS .EQ. 4) then
        x(1) = -sqrt ( ( 1.0D+00 + 2.0D+00 / sqrt( 5.0D+00 ) ) / 3.0D+00 )
        x(2) = -sqrt ( ( 1.0D+00 - 2.0D+00 / sqrt( 5.0D+00 ) ) / 3.0D+00 )
        x(3) = sqrt ( ( 1.0D+00 - 2.0D+00 / sqrt( 5.0D+00 ) ) / 3.0D+00 )
        x(4) = sqrt ( ( 1.0D+00 + 2.0D+00 / sqrt( 5.0D+00 ) ) / 3.0D+00 )
    else if (noGPS .EQ. 5) then
        x(1) = - sqrt ( ( 5.0D+00 + sqrt ( 11.0D+00 ) ) / 12.0D+00 )
        x(2) = - sqrt ( ( 5.0D+00 - sqrt ( 11.0D+00 ) ) / 12.0D+00 )
        x(3) =   0.0D+00
        x(4) =   sqrt ( ( 5.0D+00 - sqrt ( 11.0D+00 ) ) / 12.0D+00 )
        x(5) =   sqrt ( ( 5.0D+00 + sqrt ( 11.0D+00 ) ) / 12.0D+00 )
    else if (noGPS .EQ. 6) then
        x(1) = - 0.866246818107820591383598D+00
        x(2) = - 0.422518653761111529118546D+00
        x(3) = - 0.266635401516704720331534D+00
        x(4) =   0.266635401516704720331534D+00
        x(5) =   0.422518653761111529118546D+00
        x(6) =   0.866246818107820591383598D+00
    else if (noGPS .EQ. 7) then
        x(1) = - 0.883861700758049035704224D+00
        x(2) = - 0.529656775285156811385048D+00
        x(3) = - 0.323911810519907637519673D+00
        x(4) =   0.0D+00
        x(5) =   0.323911810519907637519673D+00
        x(6) =   0.529656775285156811385048D+00
        x(7) =   0.883861700758049035704224D+00
    else if (noGPS .EQ. 9) then
         x(1) = - 0.911589307728434473664949D+00
         x(2) = - 0.601018655380238071428128D+00
         x(3) = - 0.528761783057879993260181D+00
         x(4) = - 0.167906184214803943068031D+00
         x(5) =   0.0D+00
         x(6) =   0.167906184214803943068031D+00
         x(7) =   0.528761783057879993260181D+00
         x(8) =   0.601018655380238071428128D+00
         x(9) =   0.911589307728434473664949D+00
    else
       print*, "Quadrature rule is not defined for given no. of GPs."
       STOP
    endif

    w(1:noGPs) = 2.0D+00 / real ( noGPs, kind = rk )
    
    END SUBROUTINE quad_rule



!----- SORTING ------------------
SUBROUTINE sort_pick(n,arr)
IMPLICIT NONE
INTEGER, INTENT(IN) :: N
INTEGER, DIMENSION(n), INTENT(INOUT) :: arr
INTEGER :: i,j
INTEGER :: a
do j=2,n 
	a=arr(j)
	do i=j-1,1,-1 !Look for the place to insert it.
		if (arr(i) <= a) exit
		arr(i+1)=arr(i)
	end do
	arr(i+1)=a !Insert it.
end do
END SUBROUTINE sort_pick



end module mod_linspace

