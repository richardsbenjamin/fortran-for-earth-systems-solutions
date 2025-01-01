
! m_Ny attribute has been added to Config class to account for
! rectangular shapes. The nececessary adjustments in the other
! parts of the code have been made.

module NumericKinds
    implicit none

    integer, parameter :: &
        R_SP = selected_real_kind( 6, 37 ), &
        R_DP = selected_real_kind( 15, 307 ), &
        R_QP = selected_real_kind( 33, 4931 )

    integer, parameter :: RK = R_DP

    integer, parameter :: &
        I1B = selected_real_kind(2), &
        I2B = selected_real_kind(4), &
        I3B = selected_real_kind(9), &
        I4B = selected_real_kind(18)

    integer, parameter :: IK = I3B

    character(len=*), parameter :: R_SP_FMT = "f0.6", &
        R_DP_FMT = "f0.15", R_QP_FMT = "f0.33"
    character(len=*), parameter :: RK_FMT = R_DP_FMT 
end module NumericKinds

module Config_class
    use NumericKinds
    implicit none
    private

    type, public :: Config
        real(RK) :: m_diffusivity = 1.15E-6_RK, &
            m_tempA = 100._RK, &
            m_tempB = 75._RK, &
            m_tempC = 50._RK, &
            m_tempD = 25._RK, &
            m_side_length = 30._RK
        integer(IK) :: m_Nx = 200, m_Ny = 200
    end type Config

    interface Config 
        module procedure create_config
    end interface Config

contains
    type(Config) function create_config( cfg_filepath )
        character(len=*), intent(in) :: cfg_filepath
        integer :: cfg_fileID

        open( newunit=cfg_fileID, file=trim(cfg_filepath), &
            status='old', action='read')

        read(cfg_fileID,*) create_config%m_tempA
        read(cfg_fileID,*) create_config%m_tempB
        read(cfg_fileID,*) create_config%m_tempC
        read(cfg_fileID,*) create_config%m_tempD
        read(cfg_fileID,*) create_config%m_Nx
        read(cfg_fileID,*) create_config%m_Ny
        read(cfg_fileID,*) create_config%m_diffusivity
        read(cfg_fileID,*) create_config%m_side_length
        close(cfg_fileID)
    end function create_config
    
end module Config_class

module Solver_class
    use NumericKinds
    use Config_class
    implicit none
    private

    type, public :: Solver 
        private
        type(Config) :: m_config
        real(RK) :: m_Nt, &
            m_dx, m_dy, m_dt, m_A, m_B
        real(RK), allocatable, dimension(:,:) :: m_U, m_V 
        integer(IK) :: m_num_iters_max, m_curr_iter = 0
    contains
        private 
        procedure, public :: init
        procedure, public :: run
        procedure, public :: write_ascii
        procedure, public :: get_temp
        procedure :: advance_u
        procedure :: advance_v
        final :: clean_up
    end type Solver

contains
    subroutine init( this, cfg_filepath, sim_time )
        class(Solver), intent(inout) :: this 
        character(len=*), intent(in) :: cfg_filepath 
        real(RK), intent(in) :: sim_time

        integer(IK) :: nX, nY, i, j
        real(RK) :: lambda

        this%m_config = Config( cfg_filepath )
        nX = this%m_config%m_Nx
        nY = this%m_config%m_Ny

        this%m_Nt = nX**2
        this%m_num_iters_max = int( sim_time*this%m_Nt)
        this%m_dx = 1. / nX
        this%m_dy = 1. / nY
        this%m_dt = 1. / this%m_Nt 
        lambda = (2.*this%m_dt) / (this%m_dx**2)
        this%m_A = (1.-lambda) / (1.+lambda)
        this%m_B = lambda / (2.*(1.+lambda))

        allocate( this%m_U(0:nX, 0:nY), this%m_V(0:nX, 0:nY) )

        this%m_U = 1.
        this%m_U(:, nY) = [ ( 1./3.*(i/real(nX, RK))+2./3., i=0,nX) ]
        this%m_U(0, :)  = [ ( 1./3.*(j/real(nY, RK))+1./3., j=0,nY) ]
        this%m_U(:, 0)  = [ (-1./3.*(i/real(nX, RK))+1./3., i=0,Nx) ]
        this%m_U(nX, :) = [ (             j / real(nY, RK), j=0,nX) ]
        this%m_V = this%m_U 
    end subroutine init 

    real(RK) function get_temp( this, i, j )
        class(Solver), intent(in) :: this 
        integer(IK), intent(in) :: i, j
        get_temp = 0.5*( this%m_U(i,j) + this%m_V(i,j))
    end function get_temp

    subroutine run( this )
        class(Solver), intent(inout) :: this
        integer(IK) :: k

        do k=1, this%m_num_iters_max
            if( mod(k-1, (this%m_num_iters_max-1)/10) == 0) then
                write(*, '(i5,a)') nint((k*100.0)/this%m_num_iters_max), "%"
            end if 
            call this%advance_u()
            call this%advance_v()
            this%m_curr_iter = this%m_curr_iter + 1
        end do
    end subroutine run

    subroutine advance_u( this )
        class(Solver), intent(inout) :: this 
        integer(IK) :: i, j
        do j=1, this%m_config%m_Ny-1
            do i=1, this%m_config%m_Nx-1
                this%m_U(i,j) = this%m_A*this%m_U(i,j) + this%m_B*( &
                    this%m_U(i-1,j) + this%m_U(i+1,j) + this%m_U(i,j-1)+this%m_U(i,j+1) )
            end do 
        end do 
    end subroutine advance_u

    subroutine advance_v( this )
        class(Solver), intent(inout) :: this 
        integer(IK) :: i, j
        do j=this%m_config%m_Ny-1, 1, -1
            do i=this%m_config%m_Nx-1, 1, -1
                this%m_V(i,j) = this%m_A*this%m_V(i,j) + this%m_B*( &
                    this%m_V(i-1,j) + this%m_V(i+1,j) + this%m_V(i,j-1) + this%m_V(i,j+1) )
            end do 
        end do 
    end subroutine advance_v

    subroutine write_ascii( this, output_path )
        class(Solver), intent(in) :: this 
        character(len=*), intent(in) :: output_path
        integer(IK) :: x, y, outfileID

        open( newunit=outfileID, &
                file=trim(output_path), &
                status='replace', action='write')

        write(outfileID, '(a)') &
            "# output file for program heat_diffusion_v1.f90"
        write(outfileID, '(a,2x,a)') '"s"', "# time unit"
        write(outfileID, '(f0.8,2x,a)') &
            (this%m_curr_iter*this%m_config%m_side_length**2) / &
            (this%m_config%m_diffusivity*this%m_Nt), &
            "# current time"
        write(outfileID, '(a,2x,a)') '"m"', "# X unit"
        write(outfileID, '(i0,2x,a)') this%m_config%m_Nx, "# Nx"
        write(outfileID, '(a,2x,a)') '"m"', "# Y unit"
        write(outfileID, '(i0,2x,a)') this%m_config%m_Ny, "# Ny"
        write(outfileID, '(a,2x,a)') '"degree~C"', "# temperature unit"

        do x=0, this%m_config%m_Nx 
            write(outfileID, '(f0.8,2x)', advance='no') this%m_dx*this%m_config%m_side_length*x 
        end do 
        write(outfileID, '(a)') "# Xvals"
        do y=0, this%m_config%m_Ny 
            write(outfileID, '(f0.8,2x)', advance='no') this%m_dy*this%m_config%m_side_length*y
        end do 
        write(outfileID, '(a)') "# Yvals"

        write(outfileID, '(a)') "# from next line to end: simulated temperature"
        do y=0, this%m_config%m_Ny
            do x=0, this%m_config%m_Nx 
                write(outfileID, '(f0.8,2x)', advance='no') &
                    this%m_config%m_tempD+this%get_temp(x,y)*(this%m_config%m_tempA-this%m_config%m_tempD)
            end do
            write(outfileID, *)
        end do 
        close(outfileID)
    end subroutine write_ascii

    subroutine clean_up( this )
        type(Solver), intent(inout) :: this 
        deallocate( this%m_U, this%m_V )
    end subroutine clean_up

end module Solver_class

program exercise13
    use NumericKinds
    use Solver_class
    implicit none 

    type(Solver) :: square 
    real(RK) :: sim_time = 0.1

    character(len=200) :: config_file = "config_file_formatted.in", &
        output_file_name = "simulation_final_temp_field.dat"
    
    call square%init( config_file, sim_time )
    call square%run()

    call square%write_ascii( output_file_name )

end program exercise13
