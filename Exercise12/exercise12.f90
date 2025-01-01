

module Vec2d_class
    implicit none
    private

    type, public :: Vec2d
        private
        real :: m_u = 0., m_v = 0.
    contains
        private
        procedure, public :: init => init_vec2d
        procedure, public :: get_u => get_u_vec2d
        procedure, public :: set_u => set_u_vec2d
        procedure, public :: get_v => get_v_vec2d
        procedure, public :: set_v => set_v_vec2d
        procedure, public :: get_magnitude => get_magnitude_vec2d
        procedure, public :: print => print_vec2d
    end type Vec2d 

    interface Vec2d
        module procedure create_vec2d
    end interface Vec2d

    contains
        type(Vec2d) function create_vec2d( u, v )
            real, intent(in) :: u, v 
            create_vec2d%m_u = u 
            create_vec2d%m_v = v
        end function create_vec2d

        subroutine init_vec2d( this, u, v )
            class(Vec2d), intent(inout) :: this 
            real, intent(in) :: u, v 
            this%m_u = u 
            this%m_v = v 
        end subroutine init_vec2d

        real function get_u_vec2d( this )
            class(Vec2d), intent(in) :: this 
            get_u_vec2d = this%m_u
        end function get_u_vec2d

        subroutine set_u_vec2d( this, u )
            class(Vec2d), intent(inout) :: this
            real, intent(in) :: u 
            this%m_u = u 
        end subroutine set_u_vec2d

        real function get_v_vec2d( this )
            class(Vec2d), intent(inout) :: this 
            get_v_vec2d = this%m_v 
        end function get_v_vec2d

        subroutine set_v_vec2d( this, v )
            class(Vec2d), intent(inout) :: this 
            real, intent(in) :: v 
            this%m_v = v 
        end subroutine set_v_vec2d

        real function get_magnitude_vec2d( this ) result(mag)
            class(Vec2d), intent(in) :: this 
            mag = sqrt( this%m_u ** 2 + this%m_v**2 )
        end function get_magnitude_vec2d

        subroutine print_vec2d( this, name ) 
            class(Vec2d), intent(in) :: this
            character(len=*), intent(in) :: name

            write(*, '(a, ": U = ", f0.3, ", V = ", f0.3, ", Magnitude = ", f0.3)') &
                trim(name), this%m_u, this%m_v, this%get_magnitude()
        end subroutine print_vec2d

end module Vec2d_class

program ex12
    use Vec2d_class
    implicit none

    type(Vec2d) :: A, B

    A = Vec2d(u=1.1, v=9.4)
    call B%init(u=1.1, v=9.4)

    call A%print('A')
    call B%print('B')

    print '(A)', "Now we are going to change u of A and v of B and then reprint."
    call A%set_u(10.4)
    call B%set_v(20.4)

    call A%print('A')
    call B%print('B')

end program ex12



