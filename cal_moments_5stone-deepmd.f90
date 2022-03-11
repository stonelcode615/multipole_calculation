
module physics_constants
    implicit none
    real,    parameter :: e_charge      = 1.6021766208e-19        ! [C]
    real,    parameter :: Debye         = 0.208194                ! [e*A]
    real,    parameter :: pi            = 3.1415926
    real,    parameter :: B2A           = 0.529177249
    real,    parameter :: fg            = 2.9979307               ! [V*A^2/D]
end module physics_constants

module md_info
    use physics_constants
    implicit none

    integer            :: num_frame     = 0                       ! for calculation, one frame for each rerun
    real,    parameter :: L             = 13.023575/B2A           ! Box size [Angs, 1 Bohr=0.52917 Angs]
    
    character(100)     :: mole_name     = 'rerun'
    integer, parameter :: num_mole      = 64
    integer, parameter :: num_nuke      = 3                       ! nuclei number
    integer, parameter :: num_mole_cube = 4                       ! cube file number per molecule
    integer, parameter :: num_grid      = 80

    integer, parameter :: num_cube      = num_mole_cube*num_mole  ! cube file number
    real,    parameter :: delta         = L/num_grid              ! [Angstrome]

    real,    parameter :: field_point(3)= (/L*0.5,L*0.5,L*0.5/)
end module md_info

module component_info
    use md_info
    implicit none
    type mole_info                                            ! molecular info
        integer :: counter
        integer :: ndx_matched_cube(num_mole_cube)             ! index of cube files 
        
        real    :: location(3,num_nuke)                       ! molecule nuclei location xyz 
        real    :: mole_center(3)                             ! molecule center location xyz
        real    :: cube_vector(3,num_cube)                    ! cube files contributes to dipole vector (3)
        real    :: nuke_vector(3,num_nuke)
        real    :: cube_tensor(3,3,num_cube)                  ! cube files contributes to quadrupole tensor(3,3)
        real    :: nuke_tensor(3,3,num_nuke)                  ! molecule nuclei contributes to quadrupole tensor(3,3)

        real, dimension(3)   :: e_dipole,     N_dipole,     dipole
        real, dimension(3,3) :: e_quadrupole, N_quadrupole, quadrupole 

        real    :: mass

    end type mole_info
    type(mole_info) :: mole_list(num_mole)

    real, dimension(num_nuke) :: nuke_charge= (/6,1,1/)
    real, dimension(num_nuke) :: nuke_mass= (/15.9949,1.0079,1.0079/)
end module component_info

module tools
    use md_info
    use component_info
    implicit none

    contains
        subroutine pbc_mole(nuke_locations)
            ! place molecule (water) in primary pbc box, 
            ! com (oxygen) is placed in the primary box, while other nukes are assigned accordingly,
            ! no matter with the if they (ohter nukes) are beyond the primary box
            real, intent(inout) :: nuke_locations(3,num_nuke)
    
            integer :: i,j,k
            integer :: ndx_shift, ndx_nuke
            real    :: reference(3), remainning(3), delta(3)
            real    :: r0,r1,r2,increment
            
            ! step one, place reference com into the box
            reference = nuke_locations(:,1) ! oxygen
            delta = (/0.0, 0.0, 0.0/)
            do i=1, 3
                r0 = reference(i)
                do ndx_shift = -5, 5
                    r1 = r0+ndx_shift*L
                    if ( 0 < r1 .AND. r1 < L ) then
                        delta(i) = ndx_shift*L
                        exit
                    end if
                end do
            end do
            nuke_locations(:,1) = nuke_locations(:,1)+delta
            ! end step one

            ! step two, place other nuke acoordingly, no matter if its position is out of box
            do ndx_nuke = 2, num_nuke
                nuke_locations(:,ndx_nuke) = nuke_locations(:,ndx_nuke)+delta
            end do
            !reference = nuke_locations(:,1)
            !do j=1,num_nuke-1
            !    remainning = nuke_locations(:,j)
            !    do i = 1, 3
            !        r0 = remainning(i) - reference(i)
            !        do k = -5, 5
            !            increment = k*L
            !            if (abs(r0+increment)<1.6/B2A) then
            !                remainning(i) = remainning(i) + increment
            !                exit
            !            end if
            !        end do ! k
            !    end do     ! i
            !    nuke_locations(:,j) = remainning  ! update
            !end do         ! j other nuke
            
        end subroutine

        subroutine stone_inverse(a,c,n)
            integer,intent(in)  :: n
            real, intent(inout) :: a(3,3),c(3,3)
            real                :: detinv, det
            det      = (a(1,1)*a(2,2)*a(3,3)+a(1,2)*a(2,3)*a(3,1)+a(1,3)*a(2,1)*a(3,2) &
                       -a(1,1)*a(2,3)*a(3,2)-a(1,2)*a(2,1)*a(3,3)-a(1,3)*a(2,2)*a(3,1))
            detinv   = 1/(a(1,1)*a(2,2)*a(3,3)+a(1,2)*a(2,3)*a(3,1)+a(1,3)*a(2,1)*a(3,2) &
                         -a(1,1)*a(2,3)*a(3,2)-a(1,2)*a(2,1)*a(3,3)-a(1,3)*a(2,2)*a(3,1))
            c(1,1) = +detinv * (a(2,2)*a(3,3) - a(2,3)*a(3,2))
            c(2,1) = -detinv * (a(2,1)*a(3,3) - a(2,3)*a(3,1))
            c(3,1) = +detinv * (a(2,1)*a(3,2) - a(2,2)*a(3,1))
            c(1,2) = -detinv * (a(1,2)*a(3,3) - a(1,3)*a(3,2))
            c(2,2) = +detinv * (a(1,1)*a(3,3) - a(1,3)*a(3,1))
            c(3,2) = -detinv * (a(1,1)*a(3,2) - a(1,2)*a(3,1))
            c(1,3) = +detinv * (a(1,2)*a(2,3) - a(1,3)*a(2,2))
            c(2,3) = -detinv * (a(1,1)*a(2,3) - a(1,3)*a(2,1))
            c(3,3) = +detinv * (a(1,1)*a(2,2) - a(1,2)*a(2,1))
        end subroutine stone_inverse

end module tools

module cube_ana
    use physics_constants
    use md_info
    use component_info
    use tools
    implicit none

    integer :: ndx, ndx_mole,  ndx_cube, ndx_nuke
    integer :: i,j,k, dep, col, row   !(depth)

    real    :: wannier(num_grid,num_grid,num_grid)
    real    :: rho_e(num_grid,num_grid,num_grid)
    
    contains

        subroutine clean()
            do ndx_mole = 1, num_mole
                mole_list(ndx_mole)%counter = 0
            end do

        end subroutine

        subroutine moments(ndx_frame)
            integer, intent(in) :: ndx_frame

            real :: rxrx, rxry, rxrz, ryry, ryrz, rzrz
            real :: dipole(3), dipole_loc(3), quadrupole(3,3), quadrupole_loc(3,3)
            real :: ex0(3), ey0(3), ez0(3)
            real :: ex1(3), ey1(3), ez1(3), vec1(3)
            real :: magn  ! magatinude of a vector or norm
            real :: vec_oh0(3), vec_oh1(3), vec_oh2(3), crossed(3)
            real :: Q_rot(3,3),Q_rot_0(3,3), Q_rot_Inv(3,3), Q_rot_Inv_T(3,3),temp_matrix(3,3)
            integer :: ndim=3

            ex0 = (/1, 0, 0/)
            ey0 = (/0, 1, 0/)
            ez0 = (/0, 0, 1/)

            do ndx_mole = 1, num_mole
                !print*, '    '
                !print*, '-------------------------------------------------------'
                !print*, 'molecule  index:', ndx_mole
                !print*, 'cube file index:', mole_list(ndx_mole)%ndx_matched_cube(1:4)

                ! --------------------------------------calculate dipole ------------------------------
                mole_list(ndx_mole)%e_dipole     = 0
                mole_list(ndx_mole)%N_dipole     = 0
                
                do ndx = 1, num_cube
                    mole_list(ndx_mole)%e_dipole = mole_list(ndx_mole)%e_dipole +                 &
                                                   mole_list(ndx_mole)%cube_vector(:,ndx)*(-2.0)
                end do

                do ndx = 1, num_nuke
                    mole_list(ndx_mole)%N_dipole = mole_list(ndx_mole)%N_dipole + &
                                                   mole_list(ndx_mole)%nuke_vector(:,ndx)*nuke_charge(ndx)
                end do

                mole_list(ndx_mole)%dipole       = mole_list(ndx_mole)%e_dipole + &
                                                   mole_list(ndx_mole)%N_dipole

                dipole                      = mole_list(ndx_mole)%dipole/Debye*B2A
                mole_list(ndx_mole)%dipole  = dipole                                               ! [Angstrom units]

                ! -----------------------------------end calculate dipole ------------------------------

                ! -----------------------------------calculate quadrupole ------------------------------
                mole_list(ndx_mole)%e_quadrupole = 0
                mole_list(ndx_mole)%N_quadrupole = 0

                do ndx = 1, num_cube
                    mole_list(ndx_mole)%e_quadrupole = mole_list(ndx_mole)%e_quadrupole +        &
                                                       mole_list(ndx_mole)%cube_tensor(:,:,ndx)* &
                                                       0.5*(-2.0)
                end do

                do ndx = 1, num_nuke
                    mole_list(ndx_mole)%N_quadrupole = mole_list(ndx_mole)%N_quadrupole +        &
                                                       mole_list(ndx_mole)%nuke_tensor(:,:,ndx)* &
                                                       0.5*nuke_charge(ndx)
                end do

                mole_list(ndx_mole)%quadrupole       = mole_list(ndx_mole)%e_quadrupole +        &
                                                       mole_list(ndx_mole)%N_quadrupole

                quadrupole                           = mole_list(ndx_mole)%quadrupole/Debye*B2A**2 ! [Angstrom units ]
                mole_list(ndx_mole)%quadrupole       = quadrupole
                ! -----------------------------------end calculate quadrupole ------------------------------


                ! ----------------------- rotate to the local frame --------------------------
                vec_oh1 = mole_list(ndx_mole)%location(:,2)-mole_list(ndx_mole)%location(:,1)
                vec_oh2 = mole_list(ndx_mole)%location(:,3)-mole_list(ndx_mole)%location(:,1)
                !print*, 'vec_oh1', vec_oh1*B2A, norm2(vec_oh1)
                !print*, 'vec_oh2', vec_oh2*B2A, norm2(vec_oh2)

                if ( norm2(vec_oh1) >= norm2(vec_oh2)+0.000001 ) then
                    vec_oh0 = vec_oh1
                    vec_oh1 = vec_oh2
                    vec_oh2 = vec_oh0
                end if
                ex1         = vec_oh1/norm2(vec_oh1)

                crossed(1)  = vec_oh1(2)*vec_oh2(3) - vec_oh1(3)*vec_oh2(2)
                crossed(2)  = vec_oh1(3)*vec_oh2(1) - vec_oh1(1)*vec_oh2(3)
                crossed(3)  = vec_oh1(1)*vec_oh2(2) - vec_oh1(2)*vec_oh2(1)
                ez1         = crossed/norm2(crossed)

                crossed(1)  = ez1(2)*ex1(3)-ez1(3)*ex1(2)
                crossed(2)  = ez1(3)*ex1(1)-ez1(1)*ex1(3)
                crossed(3)  = ez1(1)*ex1(2)-ez1(2)*ex1(1)
                ey1         = crossed/(NORM2(crossed))                
                
                Q_rot(1,:)  = ex1
                Q_rot(2,:)  = ey1
                Q_rot(3,:)  = ez1
                Q_rot_0     = Q_rot
                !print*, 'Q_rot_0'
                !print*, Q_rot_0(1,1), Q_rot_0(1,2), Q_rot_0(1,3)
                !print*, Q_rot_0(2,1), Q_rot_0(2,2), Q_rot_0(2,3)
                !print*, Q_rot_0(3,1), Q_rot_0(3,2), Q_rot_0(3,3)

                call stone_inverse(Q_rot_0,Q_rot_Inv,ndim)
                Q_rot_Inv_T = transpose(Q_rot_Inv)
                !print*, 'Q_rot_Inv'
                !print*, Q_rot_Inv(1,1),Q_rot_Inv(1,2),Q_rot_Inv(1,3)
                !print*, Q_rot_Inv(2,1),Q_rot_Inv(2,2),Q_rot_Inv(2,3)
                !print*, Q_rot_Inv(3,1),Q_rot_Inv(3,2),Q_rot_Inv(3,3)

                dipole_loc(1)  = dot_product(dipole,Q_rot_Inv(:,1))
                dipole_loc(2)  = dot_product(dipole,Q_rot_Inv(:,2))
                dipole_loc(3)  = dot_product(dipole,Q_rot_Inv(:,3))

                quadrupole     = matmul(Q_rot_Inv_T,quadrupole)
                quadrupole_loc = matmul(quadrupole,Q_rot_Inv)

                ! convert
                mole_list(ndx_mole)%dipole     = dipole_loc
                mole_list(ndx_mole)%quadrupole = quadrupole_loc

                ! -----------------------end  rotate to the local frame --------------------------

                ! ----------------PRINT-----------------
                !print*, 'molecule  index:', ndx_mole
                !print*, 'cube file index:', mole_list(ndx_mole)%ndx_matched_cube(1:4)

                print '(6I5,11f10.5)', ndx_frame, ndx_mole, mole_list(ndx_mole)%ndx_matched_cube(1:4), & 
                    mole_list(ndx_mole)%dipole,   norm2(mole_list(ndx_mole)%dipole),                   &
                    mole_list(ndx_mole)%quadrupole(1,1), mole_list(ndx_mole)%quadrupole(1,2), mole_list(ndx_mole)%quadrupole(1,3), & 
                    mole_list(ndx_mole)%quadrupole(2,2), mole_list(ndx_mole)%quadrupole(2,3), mole_list(ndx_mole)%quadrupole(3,3), &
                    mole_list(ndx_mole)%quadrupole(1,1)+mole_list(ndx_mole)%quadrupole(2,2)+mole_list(ndx_mole)%quadrupole(3,3)

            end do       ! loop over water indices

        end subroutine moments


        subroutine match(ndx_cube,got_issue)

            integer, intent(in) ::  ndx_cube
            integer, intent(inout) :: got_issue
            real     :: rx_grid,ry_grid,rz_grid
            real     :: rx,ry,rz
            complex  :: integration_d(3)
            real     :: rxrx,rxry,rxrz, ryry,ryrz,rzrz
            complex  :: integration_q1(3,3)
            complex  :: integration_q2(3,3)

            real     :: r2, r0, increment

            complex  :: phase(3)

            real     :: cube_center(3), mole_shifted_center(3), mole_center(3)
            integer  :: counter

            real     :: cube_vector(3)
            real     :: cube_tensor(3,3)
            real     :: nuke_vector(3)
            real     :: nuke_tensor(3,3)

            real     :: omega
            real     :: sum_rho
            
            integer  :: accumulate_num_mole
            real     :: com_wfc_dist, com_wfc_dist_record
            integer  :: ndx_mole_record
            real     :: diff, diff_vec(3)

            rho_e = wannier**2
            sum_rho = 0.0

            integration_d  = 0
            integration_q1 = 0
            integration_q2 = 0
            do row =1, num_grid
                do col =1, num_grid
                    do dep =1, num_grid
                        rz_grid  = (dep-1)*delta  ! [A] over the box volume
                        ry_grid  = (col-1)*delta
                        rx_grid  = (row-1)*delta
                        phase(1) = cmplx(0, 2*pi/L*rx_grid)
                        phase(2) = cmplx(0, 2*pi/L*ry_grid)
                        phase(3) = cmplx(0, 2*pi/L*rz_grid)

                        sum_rho  = sum_rho +rho_e(dep,col,row)*delta**3

                        integration_d(1)   = integration_d(1)     + rho_e(dep,col,row)*delta**3*exp(phase(1))
                        integration_d(2)   = integration_d(2)     + rho_e(dep,col,row)*delta**3*exp(phase(2))
                        integration_d(3)   = integration_d(3)     + rho_e(dep,col,row)*delta**3*exp(phase(3))

                        integration_q1(1,1) = integration_q1(1,1) + rho_e(dep,col,row)*delta**3*exp(phase(1))*exp(-phase(1))
                        integration_q2(1,1) = integration_q2(1,1) + rho_e(dep,col,row)*delta**3*exp(phase(1))*exp(phase(1))
                        integration_q1(1,2) = integration_q1(1,2) + rho_e(dep,col,row)*delta**3*exp(phase(1))*exp(-phase(2))
                        integration_q2(1,2) = integration_q2(1,2) + rho_e(dep,col,row)*delta**3*exp(phase(1))*exp(phase(2))
                        integration_q1(1,3) = integration_q1(1,3) + rho_e(dep,col,row)*delta**3*exp(phase(1))*exp(-phase(3))
                        integration_q2(1,3) = integration_q2(1,3) + rho_e(dep,col,row)*delta**3*exp(phase(1))*exp(phase(3))
                        integration_q1(2,2) = integration_q1(2,2) + rho_e(dep,col,row)*delta**3*exp(phase(2))*exp(-phase(2))
                        integration_q2(2,2) = integration_q2(2,2) + rho_e(dep,col,row)*delta**3*exp(phase(2))*exp(phase(2))
                        integration_q1(2,3) = integration_q1(2,3) + rho_e(dep,col,row)*delta**3*exp(phase(2))*exp(-phase(3))
                        integration_q2(2,3) = integration_q2(2,3) + rho_e(dep,col,row)*delta**3*exp(phase(2))*exp(phase(3))
                        integration_q1(3,3) = integration_q1(3,3) + rho_e(dep,col,row)*delta**3*exp(phase(3))*exp(-phase(3))
                        integration_q2(3,3) = integration_q2(3,3) + rho_e(dep,col,row)*delta**3*exp(phase(3))*exp(phase(3))
                        
                    end do
                end do
            end do 
            
            do i = 1, 3
                ! Wannier center coodinates, x,y,z [A]
                r0  = L/(2*pi)*imag(log(integration_d(i)))        ! [A, Angstrome]
                cube_center(i) = r0
            end do

            counter             = 0
            accumulate_num_mole = 0
            ndx_mole_record     = 0
            com_wfc_dist_record = L

            do ndx_mole =1, num_mole                                 ! loop over molecule index

                accumulate_num_mole = accumulate_num_mole + 1

                if ( mole_list(ndx_mole)%counter == 4) then
                    cycle
                end if 

                mole_center = mole_list(ndx_mole)%mole_center
                diff_vec    = mole_center-cube_center
                where (diff_vec > 0.5*L)
                    mole_shifted_center = mole_center-L
                else where (-0.5*L <= diff_vec .AND. diff_vec <= 0.5*L)
                    mole_shifted_center = mole_center
                else where
                    mole_shifted_center = mole_center+L
                end where

                if (any(mole_shifted_center-cube_center > 0.8/B2A)) then
                    cycle
                end if

                com_wfc_dist   = NORM2(mole_shifted_center-cube_center)
                if (com_wfc_dist_record >= com_wfc_dist) then
                    ndx_mole_record     = ndx_mole
                    com_wfc_dist_record = com_wfc_dist
                end if
            end do  ! end loop over mole ndx

            !print*, 'cube ', ndx_cube, 'center: ', cube_center*B2A, 'mole ', ndx_mole_record, &
            !    'dist ', com_wfc_dist_record*B2A, 'counter', counter
            if ( com_wfc_dist_record*B2A <= 0.8) then
                counter = mole_list(ndx_mole_record)%counter
                counter = counter+1
                mole_list(ndx_mole_record)%counter                   = counter
                mole_list(ndx_mole_record)%ndx_matched_cube(counter) = ndx_cube
            else
                print*, 'Issue No.1: Cube Index ', ndx_cube, ' com_wfc_dist > 0.8 Angstrom'
                got_issue = 1
                return
            end if

            ! this loop is to place the cube center close to the com, no matter with pbc box
            mole_center = mole_list(ndx_mole_record)%mole_center
            do i = 1, 3
                r0 = cube_center(i)-mole_center(i)
                do k = -1, 1
                    increment = k*L
                    if ( abs(r0+increment)<0.5*L ) then
                        cube_center(i) = cube_center(i)+increment
                        exit
                    end if
                end do
            end do            

            rx   = cube_center(1)
            ry   = cube_center(2)
            rz   = cube_center(3)
            cube_center = (/rx, ry, rz/)     ! [Angs]

            rxrx = L**2/(16*pi**2)*(  log(abs(integration_q1(1,1))**2)    &
                                      -log(abs(integration_q2(1,1))**2) ) &
                    +rx*rx
            rxry = L**2/(16*pi**2)*(  log(abs(integration_q1(1,2))**2)    &
                                      -log(abs(integration_q2(1,2))**2) ) &
                    +rx*ry
            rxrz = L**2/(16*pi**2)*(  log(abs(integration_q1(1,3))**2)    &
                                      -log(abs(integration_q2(1,3))**2) ) &
                    +rx*rz
            ryry = L**2/(16*pi**2)*(  log(abs(integration_q1(2,2))**2)    &
                                      -log(abs(integration_q2(2,2))**2) ) &
                    +ry*ry
            ryrz = L**2/(16*pi**2)*(  log(abs(integration_q1(2,3))**2)    &
                                      -log(abs(integration_q2(2,3))**2) ) &
                    +ry*rz
            rzrz = L**2/(16*pi**2)*(  log(abs(integration_q1(3,3))**2)    &
                                      -log(abs(integration_q2(3,3))**2) ) &
                    +rz*rz

            omega = 2/(2*pi)**2*L**2*(1-abs(integration_d(1))+&
                                      1-abs(integration_d(2))+&
                                      1-abs(integration_d(3)))

            rxrx  = rxrx - mole_center(1)*rx - mole_center(1)*rx + mole_center(1)*mole_center(1)
            rxry  = rxry - mole_center(1)*ry - mole_center(2)*rx + mole_center(1)*mole_center(2)
            rxrz  = rxrz - mole_center(1)*rz - mole_center(3)*rx + mole_center(1)*mole_center(3)
            ryry  = ryry - mole_center(2)*ry - mole_center(2)*ry + mole_center(2)*mole_center(2)
            ryrz  = ryrz - mole_center(2)*rz - mole_center(3)*ry + mole_center(2)*mole_center(3)
            rzrz  = rzrz - mole_center(3)*rz - mole_center(3)*rz + mole_center(3)*mole_center(3)

            r2    = rxrx+ryry+rzrz

            cube_vector      = cube_center-mole_center
            cube_tensor(:,1) = (/3*rxrx-r2, 3*rxry,    3*rxrz/)
            cube_tensor(:,2) = (/3*rxry,    3*ryry-r2, 3*ryrz/)
            cube_tensor(:,3) = (/3*rxrz,    3*ryrz,    3*rzrz-r2/)

            mole_list(ndx_mole_record)%cube_vector(:,  counter) = cube_vector
            mole_list(ndx_mole_record)%cube_tensor(:,:,counter) = cube_tensor

        end subroutine match

        subroutine molecule_read(cube_name)
            character(100), intent(in) :: cube_name
            character(len=100) :: buffer_str
            integer :: buffer_int
            real    :: buffer_real

            open(unit=10, file=cube_name,status='old', action='read')

            !read head
            do ndx = 1, 6
                read(10,'(A)') buffer_str
            end do
            !read water xyz
            do ndx_mole = 1, num_mole
                do col = 1,num_nuke
                    read(10,*) buffer_int, buffer_real, mole_list(ndx_mole)%location(:,col)
                                             ! column wise read, row loop first
                    mole_list(ndx_mole)%location(:,col) = mole_list(ndx_mole)%location(:,col)
                end do
            end do

            close(10)

        end subroutine molecule_read

        subroutine molecule_comp()
            real :: molecule_mass
            real :: mole_center(3)
            real :: rx_n, ry_n, rz_n, r2
            real :: nuke_vector(3)
            real :: nuke_tensor(3,3)

            do ndx_mole = 1, num_mole

                if (any(mole_list(ndx_mole)%location(:,1)<0) .OR. &
                    any(mole_list(ndx_mole)%location(:,1)>L)) then
                    call pbc_mole(mole_list(ndx_mole)%location)  ! call pbc mole to place com (O) atom in the box and then place ohters accordingly
                end if

                mole_list(ndx_mole)%mole_center = mole_list(ndx_mole)%location(:,1)   ! stone regards oxygen as ce
                !print*, 'mole_ndx:', ndx_mole, mole_list(ndx_mole)%mole_center*B2A  !  2020

                do ndx_nuke = 1, num_nuke
                    rx_n = mole_list(ndx_mole)%location(1,ndx_nuke)-mole_list(ndx_mole)%mole_center(1)
                    ry_n = mole_list(ndx_mole)%location(2,ndx_nuke)-mole_list(ndx_mole)%mole_center(2)
                    rz_n = mole_list(ndx_mole)%location(3,ndx_nuke)-mole_list(ndx_mole)%mole_center(3)

                    nuke_vector      = (/rx_n,ry_n,rz_n/)
                    r2               = rx_n*rx_n+ry_n*ry_n+rz_n*rz_n
                    nuke_tensor(:,1) = (/3*rx_n*rx_n-r2, 3*rx_n*ry_n,    3*rx_n*rz_n/)
                    nuke_tensor(:,2) = (/3*ry_n*rx_n,    3*ry_n*ry_n-r2, 3*ry_n*rz_n/)
                    nuke_tensor(:,3) = (/3*rz_n*rx_n,    3*rz_n*ry_n,    3*rz_n*rz_n-r2/)
                    mole_list(ndx_mole)%nuke_vector(:,ndx_nuke)   = nuke_vector
                    mole_list(ndx_mole)%nuke_tensor(:,:,ndx_nuke) = nuke_tensor
                end do

            end do   ! end of loop over molecule index

        end subroutine molecule_comp

        subroutine cube_read(cube_name)
            character(100),intent(in) :: cube_name
            character(len=100) :: buffer_str
            integer :: buffer_int
            real    :: buffer_real

            open(unit=11, file=cube_name,status='old')

            !read head
            !read water xyz
            do ndx = 1, 6+num_mole*num_nuke
                read(11,'(A)') buffer_str
            end do

            !read wannier cube
            do row = 1, num_grid
                do col = 1, num_grid
                    do dep = 1, ceiling(num_grid/float(6))
                        if (dep == ceiling(num_grid/float(6))) then
                            if (mod(num_grid,6)==0) then
                                read(11,*) wannier(dep*6-5:dep*6,col,row)
                            else
                                read(11,*) wannier(num_grid-mod(num_grid,6)+1:num_grid,col,row)
                            end if
                        else
                            read(11,*) wannier(dep*6-5:dep*6,col,row)
                        end if
                    end do 
                end do
            end do

            close(11)

        end subroutine cube_read

end module cube_ana

program main
    use physics_constants
    use component_info
    use cube_ana
    implicit none

    character(100) :: cube_file, cube_ndx, frame_ndx, frame_ndx0, frame_ndx1
    integer        :: ndx_frame, ndx_frame0, ndx_frame1
    integer        :: got_issue

    call GET_COMMAND_ARGUMENT(1,frame_ndx0)
    call GET_COMMAND_ARGUMENT(2,frame_ndx1)
    read(frame_ndx0,*) ndx_frame0
    read(frame_ndx1,*) ndx_frame1

    do ndx_frame = ndx_frame0, ndx_frame1
        write(frame_ndx, "(I10)") ndx_frame

        got_issue=0

        do ndx_mole = 1, num_mole
            mole_list(ndx_mole)%counter = 0
        end do

        do ndx_cube = 1, num_cube
            write(cube_ndx, "(I10)") ndx_cube
            cube_file="./cubes/"//trim(adjustl(mole_name))//"-locHOMO_w"//&
                                  trim(adjustl(cube_ndx))//"_s1-1_"//&
                                  trim(adjustl(frame_ndx))//".cube"

            if (ndx_cube == 1) then
                call molecule_read(cube_file)
                call molecule_comp()
            end if
            call cube_read(cube_file)

            call match(ndx_cube,got_issue)
            if (got_issue == 1) then
                exit
            end if   

        end do

        if (got_issue == 1) then
            print*, 'STOP, frame', ndx_frame, 'report the above issue to be resolved'
            print*, ' '
            cycle
        end if
        call moments(ndx_frame)
        call clean()

    end do

end program main
