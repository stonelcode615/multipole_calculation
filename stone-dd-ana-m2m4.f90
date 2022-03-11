
module physics_constants
    implicit none
    real,    parameter :: e_charge       = 1.6021766208e-19        ! [C]
    real,    parameter :: Debye          = 0.208194                ! [e*A]
    real,    parameter :: pi             = 3.1415926
    real,    parameter :: B2A            = 0.529177249
    real,    parameter :: fg             = 2.9979307               ! [V*A^2/D]
end module physics_constants

module md_info
    use physics_constants
    implicit none

    integer, parameter :: num_frame     = 60000
    real,    parameter :: L             = 13.023575 !25.013800    ! Box size [Angs, 1 Bohr=0.52917 Angs]
    character(100)     :: mole_name     = 'water'
    integer, parameter :: num_mole      = 64
    integer, parameter :: num_nuke      = 3                       ! nuclei number
    
    real,    parameter :: cave_location(3)= (/L*0.5,L*0.5,L*0.5/)

    integer            :: ndx_frame0, ndx_frame1

    real               :: O_location(3,num_mole)
    real               :: H_location(3,num_mole*2)
    integer            :: num_seperate_H
    integer            :: ndx_H_assigned(2*num_mole) 
    real               :: dist_O_H_assigned(num_mole*2)

end module md_info

module component_info
    use md_info
    implicit none
    type mole_info                                            ! molecular info
        
        real    :: location(3,num_nuke)                       ! molecule nuclei location xyz 
        real    :: dipole(3)
        real    :: quadrupole(3,3)

        real    :: mass
        real    :: dist_to_cave
        real    :: T_vector(3)
        real    :: T_tensor(3,3)

        real    :: phi2 ! dipole     potential
        real    :: phi4 ! quadrupole potential

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
        subroutine build_water(ndx_mole)
            integer, intent(in)  :: ndx_mole
            integer              :: ndx_H, ndx_H_1st, ndx_H_2nd
            real                 :: OH_vector(3), H_location_shifted(3)
            real                 :: OH_dist, OH_dist_1st, OH_dist_2nd
            real                 :: H1_location(3), H2_location(3)

            OH_dist_1st  = 1000.00
            OH_dist_2nd  = 1000.00
            ndx_H_1st = num_mole*2+1
            ndx_H_2nd = num_mole*2+2

            do ndx_H=1 , num_mole*2
                OH_vector   = H_location(:,ndx_H) - O_location(:,ndx_mole)
                where (OH_vector > 0.5*L)
                    H_location_shifted = H_location(:,ndx_H) - L
                else where (-0.5*L <= OH_vector .AND. OH_vector <= 0.5*L)
                    H_location_shifted = H_location(:,ndx_H)
                else where
                    H_location_shifted = H_location(:,ndx_H) + L
                end where

                OH_vector   = H_location_shifted - O_location(:,ndx_mole)
                OH_dist     = NORM2(OH_vector)
                if ( OH_dist < OH_dist_1st ) then 
                    OH_dist_2nd  = OH_dist_1st
                    OH_dist_1st  = OH_dist
                    ndx_H_2nd    = ndx_H_1st
                    ndx_H_1st    = ndx_H
                    H2_location  = H1_location
                    H1_location  = H_location_shifted
                else if ( OH_dist < OH_dist_2nd ) then
                     OH_dist_2nd = OH_dist
                     ndx_H_2nd   = ndx_H
                     H2_location = H_location_shifted
                else 
                    cycle
                end if
            end do  ! end loop over H atom index
            !print*, 'ndx_mole = ', ndx_mole, 'ndx_H_1st = ', ndx_H_1st, 'ndx_H_2nd = ', ndx_H_2nd

            !if ( ANY(ndx_H_assigned .eq. ndx_H_1st) .or. ANY(ndx_H_assigned .eq. ndx_H_2nd)) then 
            !    print*, ' '
            !    print*, 'The water building is crashed '
            !    print*, 'ndx  O: ', ndx_mole, O_location(:,ndx_mole)
            !    print*, 'ndx H1: ', 2*ndx_mole-1, H_location(:, 2*ndx_mole-1), 'distance ',&
            !    NORM2(H_location(:,2*ndx_mole-1)-O_location(:,ndx_mole))
            !    print*, 'ndx H2: ', 2*ndx_mole, H_location(:, 2*ndx_mole), 'distance ',&
            !    NORM2(H_location(:,2*ndx_mole)-O_location(:,ndx_mole))
            !    print*, 'H1  H2:  distance ', NORM2(H_location(:,2*ndx_mole)-H_location(:,2*ndx_mole-1))
            !    print*, 'ndx_H_1st ', ndx_H_1st, H_location(:,ndx_H_1st), 'distance ',&
            !    NORM2(H_location(:,ndx_H_1st)-O_location(:,ndx_mole))    
            !    print*, 'ndx_H_2nd ', ndx_H_2nd, H_location(:,ndx_H_2nd), 'distance ',&
            !    NORM2(H_location(:,ndx_H_2nd)-O_location(:,ndx_mole))    
            !    print*, '1st 2nd: distance ', NORM2(H_location(:,ndx_H_1st)-H_location(:,ndx_H_2nd))

            !end if

            ndx_H_assigned(2*ndx_mole-1)      = ndx_H_1st
            ndx_H_assigned(2*ndx_mole)        = ndx_H_2nd
            dist_O_H_assigned(2*ndx_mole-1)   = OH_dist_1st
            dist_O_H_assigned(2*ndx_mole)     = OH_dist_2nd

            mole_list(ndx_mole)%location(:,1) = O_location(:,ndx_mole)
            mole_list(ndx_mole)%location(:,2) = H1_location
            mole_list(ndx_mole)%location(:,3) = H2_location


        end subroutine

        subroutine ana_seperate_H(ndx_H)
            integer, intent(in)  :: ndx_H
            real                 :: OH_vector(3), O_location_shifted(3)
            integer              :: ndx_mole, ndx_H_1st
            real                 :: OH_dist,  OH_dist_1st

            OH_dist_1st  = 1000.00
            ndx_H_1st = num_mole
            do ndx_mole = 1 , num_mole
                OH_vector   = H_location(:,ndx_H) - O_location(:,ndx_mole)
                where (OH_vector > 0.5*L)
                    O_location_shifted = O_location(:,ndx_mole) - L
                else where (-0.5*L <= OH_vector .AND. OH_vector <= 0.5*L)
                    O_location_shifted = O_location(:,ndx_mole)
                else where
                    O_location_shifted = O_location(:,ndx_mole) + L
                end where

                OH_vector   = H_location(:,ndx_H) - O_location_shifted
                OH_dist     = NORM2(OH_vector)
                if ( OH_dist < OH_dist_1st) then
                    OH_dist_1st  = OH_dist
                    ndx_H_1st    = ndx_mole
                end if
            end do
            print*, 'closest O is ', ndx_H_1st, O_location(:,ndx_H_1st),    'dist',  OH_dist_1st
            print*, 'the O 1st H assigned ', ndx_H_assigned(2*ndx_H_1st-1), 'dist ', dist_O_H_assigned(2*ndx_H_1st-1)
            print*, 'the O 2nd H assigned ', ndx_H_assigned(2*ndx_H_1st),   'dist ', dist_O_H_assigned(2*ndx_H_1st)
 
        end subroutine
            
        subroutine pbc_mole(nuke_locations)
            ! place molecule (water) in primary pbc box, 
            ! com (oxygen) is placed in the primary box, while other nukes are assigned accordingly,
            ! no matter with the if they (ohter nukes) are beyond the primary box
            real, intent(inout) :: nuke_locations(3,num_nuke)
    
            integer :: i
            integer :: ndx_shift, ndx_nuke
            real    :: reference(3),  delta(3)
            real    :: r0, r1
            
            ! step one,1place reference com into the box
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
            
        end subroutine

        subroutine stone_inverse(a,c)
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

        subroutine Rotate_from_loc_to_lab(ndx_mole)
            integer,intent(in)    :: ndx_mole
            real :: ex0(3), ey0(3), ez0(3)
            real :: ex1(3), ey1(3), ez1(3)
            real :: vec_oh1(3), vec_oh2(3)
            real :: crossed(3)

            real :: dipole(3), dipole_loc(3)
            real :: quadrupole(3,3), quadrupole_loc(3,3)
            real :: Q_rot(3,3), Q_rot_Inv(3,3), Q_rot_T(3,3)

            ex0 = (/1, 0, 0/)
            ey0 = (/0, 1, 0/)
            ez0 = (/0, 0, 1/)

            vec_oh1 = mole_list(ndx_mole)%location(:,2)-mole_list(ndx_mole)%location(:,1)
            vec_oh2 = mole_list(ndx_mole)%location(:,3)-mole_list(ndx_mole)%location(:,1)
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

            dipole_loc      = mole_list(ndx_mole)%dipole
            quadrupole_loc  = mole_list(ndx_mole)%quadrupole

            dipole(1)       = dot_product(dipole_loc,Q_rot(:,1))
            dipole(2)       = dot_product(dipole_loc,Q_rot(:,2))
            dipole(3)       = dot_product(dipole_loc,Q_rot(:,3))

            Q_rot_T         = transpose(Q_rot)
            quadrupole_loc  = matmul(Q_rot_T,quadrupole_loc)
            quadrupole      = matmul(quadrupole_loc,Q_rot)

            mole_list(ndx_mole)%dipole     = dipole
            mole_list(ndx_mole)%quadrupole = quadrupole
                
        end subroutine Rotate_from_loc_to_lab

end module tools

module compute
    use physics_constants
    use md_info
    use component_info
    use tools
    implicit none

    contains

        subroutine comp_Ttensor(ndx_mole)
            integer, intent(in) :: ndx_mole
            real                :: dist_R
            real                :: dist_vector(3)
            real                :: dist_tensor(3,3)
            real,    dimension(3,3) :: Kronecker 
            Kronecker(1,:) = (/ 1.0,0.0,0.0 /)
            Kronecker(2,:) = (/ 0.0,1.0,0.0 /)
            Kronecker(3,:) = (/ 0.0,0.0,1.0 /)
            
            dist_vector      = cave_location-mole_list(ndx_mole)%location(:,1)
            dist_tensor(:,1) = dist_vector*dist_vector(1)
            dist_tensor(:,2) = dist_vector*dist_vector(2)
            dist_tensor(:,3) = dist_vector*dist_vector(3)

            dist_R                           = NORM2(dist_vector)
            mole_list(ndx_mole)%dist_to_cave = dist_R
            mole_list(ndx_mole)%T_vector     = -1*(dist_vector)/dist_R**3 
            mole_list(ndx_mole)%T_tensor     = (3*dist_tensor-Kronecker*dist_R**2)/dist_R**5

        end subroutine comp_Ttensor

end module compute


program main
    use physics_constants
    use md_info
    use component_info
    use tools
    use compute
    implicit none

    integer :: ndx_frame, ndx_mole, ndx_record
    integer :: buf_int(50)
    real    :: buf_real(50)
    character(len=100) :: buf_str0, buf_str(50)

    integer, parameter :: num_bin         = 30

    integer :: bin_index     
    real    :: bin_size                   = 1.0
    integer :: bin_count(num_bin)         = 0
    real    :: bin_sum(num_bin)           = 0.0
    real    :: bin_ave(num_bin,num_frame) = 0.0

    real    :: bin_sum_px(num_bin)            = 0.0
    real    :: bin_ave_px(num_bin,num_frame)  = 0.0
    real    :: bin_sum_py(num_bin)            = 0.0
    real    :: bin_ave_py(num_bin,num_frame)  = 0.0
    real    :: bin_sum_pz(num_bin)            = 0.0
    real    :: bin_ave_pz(num_bin,num_frame)  = 0.0
    real    :: bin_sum_qxx(num_bin)           = 0.0
    real    :: bin_ave_qxx(num_bin,num_frame) = 0.0
    real    :: bin_sum_qyy(num_bin)           = 0.0
    real    :: bin_ave_qyy(num_bin,num_frame) = 0.0
    real    :: bin_sum_qzz(num_bin)           = 0.0
    real    :: bin_ave_qzz(num_bin,num_frame) = 0.0
    real    :: bin_sum_qxy(num_bin)           = 0.0
    real    :: bin_ave_qxy(num_bin,num_frame) = 0.0
    real    :: bin_sum_qxz(num_bin)           = 0.0
    real    :: bin_ave_qxz(num_bin,num_frame) = 0.0
    real    :: bin_sum_qyz(num_bin)           = 0.0
    real    :: bin_ave_qyz(num_bin,num_frame) = 0.0

    real    :: phi2_sum(num_frame)        = 0.0
    real    :: phi4_sum(num_frame)        = 0.0
    real    :: phi2_ave, phi4_ave

    real    :: bin_phi2_sum(num_bin)      = 0.0
    real    :: bin_phi4_sum(num_bin)      = 0.0
    real    :: bin_phi2_ave(num_bin,num_frame) = 0.0
    real    :: bin_phi4_ave(num_bin,num_frame) = 0.0

    real    :: bin_aof(num_bin)      = 0.0 
    real    :: bin_aof_px(num_bin)   = 0.0 
    real    :: bin_aof_py(num_bin)   = 0.0 
    real    :: bin_aof_pz(num_bin)   = 0.0 
    real    :: bin_aof_qxx(num_bin)  = 0.0 
    real    :: bin_aof_qyy(num_bin)  = 0.0 
    real    :: bin_aof_qzz(num_bin)  = 0.0 
    real    :: bin_aof_qxy(num_bin)  = 0.0 
    real    :: bin_aof_qxz(num_bin)  = 0.0 
    real    :: bin_aof_qyz(num_bin)  = 0.0 

    character(100) :: frame_ndx0, frame_ndx1

    call GET_COMMAND_ARGUMENT(1,frame_ndx0)
    call GET_COMMAND_ARGUMENT(2,frame_ndx1)
    read(frame_ndx0,*) ndx_frame0
    read(frame_ndx1,*) ndx_frame1
    
    open(unit=10,file='water.xyz', status='old',action='read')
    open(unit=11,file='dipole_loc',    status='old',action='read')
    open(unit=12,file='quadrupole_loc',status='old',action='read')
    open(unit=21,file='m2m4', status='replace',action='write')

    do ndx_frame = 1, ndx_frame0-1 
        do ndx_mole = 1, num_mole*3+3
            read(10,'(A)') buf_str0
        end do
    end do

    phi2_sum   = 0.0
    phi4_sum   = 0.0
    do ndx_frame = ndx_frame0, ndx_frame1
        ndx_record = ndx_frame
        
        read(10,*) buf_str(1)
        read(10,*) buf_str(1)
        do ndx_mole = 1, num_mole
            read(10,*) buf_str0, O_location(:,ndx_mole)
            read(10,*) buf_str0, H_location(:,2*ndx_mole-1)
            read(10,*) buf_str0, H_location(:,2*ndx_mole)
            read(11,*) buf_int(1:2), mole_list(ndx_mole)%dipole
            read(12,*) buf_int(1:2), buf_real(1:9)
            mole_list(ndx_mole)%quadrupole(1,1) = buf_real(1)
            mole_list(ndx_mole)%quadrupole(1,2) = buf_real(2)
            mole_list(ndx_mole)%quadrupole(1,3) = buf_real(3)
            mole_list(ndx_mole)%quadrupole(2,1) = buf_real(4)
            mole_list(ndx_mole)%quadrupole(2,2) = buf_real(5)
            mole_list(ndx_mole)%quadrupole(2,3) = buf_real(6)
            mole_list(ndx_mole)%quadrupole(3,1) = buf_real(7)
            mole_list(ndx_mole)%quadrupole(3,2) = buf_real(8)
            mole_list(ndx_mole)%quadrupole(3,3) = buf_real(9)
        end do
        read(10,*) buf_str(1)
    
        bin_count      = 0
        bin_sum        = 0.0
        ndx_H_assigned = 0
        bin_phi2_sum = 0.0
        bin_phi4_sum = 0.0
        do ndx_mole = 1, num_mole
            call build_water(ndx_mole)
            !call Rotate_from_loc_to_lab(ndx_mole)
            call comp_Ttensor(ndx_mole)

            mole_list(ndx_mole)%phi2 = -1*fg*&
                dot_product(mole_list(ndx_mole)%T_vector,mole_list(ndx_mole)%dipole)
            phi2_sum(ndx_record) = phi2_sum(ndx_record) + mole_list(ndx_mole)%phi2
            !print*, 'frame', ndx_frame, 'mole', ndx_mole, mole_list(ndx_mole)%dipole, norm2(mole_list(ndx_mole)%dipole)

            mole_list(ndx_mole)%phi4 =                            fg*&
                dot_product(mole_list(ndx_mole)%T_tensor(:,1),mole_list(ndx_mole)%quadrupole(:,1))
            mole_list(ndx_mole)%phi4 = mole_list(ndx_mole)%phi4 + fg*&
                dot_product(mole_list(ndx_mole)%T_tensor(:,2),mole_list(ndx_mole)%quadrupole(:,2))
            mole_list(ndx_mole)%phi4 = mole_list(ndx_mole)%phi4 + fg*&
                dot_product(mole_list(ndx_mole)%T_tensor(:,3),mole_list(ndx_mole)%quadrupole(:,3))
            phi4_sum(ndx_record) = phi4_sum(ndx_record) + mole_list(ndx_mole)%phi4
            !print*, 'frame', ndx_frame, 'mole', ndx_mole 
            !print*, mole_list(ndx_mole)%quadrupole(:,1)
            !print*, mole_list(ndx_mole)%quadrupole(:,2)
            !print*, mole_list(ndx_mole)%quadrupole(:,3)

            bin_index = floor(mole_list(ndx_mole)%dist_to_cave/bin_size)
            bin_count(bin_index)  = bin_count(bin_index) + 1
            bin_sum(bin_index)    = bin_sum(bin_index) + norm2(mole_list(ndx_mole)%dipole)
            bin_sum_px(bin_index) = bin_sum_px(bin_index)  + mole_list(ndx_mole)%dipole(1)
            bin_sum_py(bin_index) = bin_sum_py(bin_index)  + mole_list(ndx_mole)%dipole(2)
            bin_sum_pz(bin_index) = bin_sum_pz(bin_index)  + mole_list(ndx_mole)%dipole(3)
            bin_sum_qxx(bin_index)= bin_sum_qxx(bin_index) + mole_list(ndx_mole)%quadrupole(1,1)
            bin_sum_qyy(bin_index)= bin_sum_qyy(bin_index) + mole_list(ndx_mole)%quadrupole(2,2)
            bin_sum_qzz(bin_index)= bin_sum_qzz(bin_index) + mole_list(ndx_mole)%quadrupole(3,3)
            bin_sum_qxy(bin_index)= bin_sum_qxy(bin_index) + mole_list(ndx_mole)%quadrupole(1,2)
            bin_sum_qxz(bin_index)= bin_sum_qxz(bin_index) + mole_list(ndx_mole)%quadrupole(1,3)
            bin_sum_qyz(bin_index)= bin_sum_qyz(bin_index) + mole_list(ndx_mole)%quadrupole(2,3)

            bin_phi2_sum(bin_index) = bin_phi2_sum(bin_index) + mole_list(ndx_mole)%phi2
            bin_phi4_sum(bin_index) = bin_phi4_sum(bin_index) + mole_list(ndx_mole)%phi4
        !    print*, 'ndx_frame = ', ndx_frame, ndx_mole, bin_index, bin_count(bin_index)
        end do
        !print*, 'frame ', ndx_frame, 'phi2 ', phi2_sum(ndx_record), 'phi4 ', phi4_sum(ndx_record)

        do bin_index = 1, num_bin
            if ( bin_count(bin_index) .ne. 0 )  then
                bin_ave(bin_index,ndx_record)     = bin_sum(bin_index)/bin_count(bin_index)
                bin_ave_px(bin_index,ndx_record)  = bin_sum_px(bin_index)/bin_count(bin_index)
                bin_ave_py(bin_index,ndx_record)  = bin_sum_py(bin_index)/bin_count(bin_index)
                bin_ave_pz(bin_index,ndx_record)  = bin_sum_pz(bin_index)/bin_count(bin_index)
                bin_ave_qxx(bin_index,ndx_record) = bin_sum_qxx(bin_index)/bin_count(bin_index)
                bin_ave_qyy(bin_index,ndx_record) = bin_sum_qyy(bin_index)/bin_count(bin_index)
                bin_ave_qzz(bin_index,ndx_record) = bin_sum_qzz(bin_index)/bin_count(bin_index)
                bin_ave_qxy(bin_index,ndx_record) = bin_sum_qxy(bin_index)/bin_count(bin_index)
                bin_ave_qxz(bin_index,ndx_record) = bin_sum_qxz(bin_index)/bin_count(bin_index)
                bin_ave_qyz(bin_index,ndx_record) = bin_sum_qyz(bin_index)/bin_count(bin_index)

                bin_phi2_ave(bin_index,ndx_record) = bin_phi2_sum(bin_index)
                bin_phi4_ave(bin_index,ndx_record) = bin_phi4_sum(bin_index)
            end if
        end do

        print*, bin_ave(5,ndx_record)
        print*, bin_ave_px(5,ndx_record), bin_ave_py(5,ndx_record), bin_ave_pz(5,ndx_record)
        print*, bin_ave_qxx(5,ndx_record),bin_ave_qyy(5,ndx_record),bin_ave_qzz(5,ndx_record)
        print*, bin_ave_qxy(5,ndx_record),bin_ave_qxz(5,ndx_record),bin_ave_qyz(5,ndx_record)
        print*, '      '
        bin_sum      = 0.0 
        bin_sum_px   = 0.0 
        bin_sum_py   = 0.0 
        bin_sum_pz   = 0.0 
        bin_sum_qxx  = 0.0 
        bin_sum_qyy  = 0.0 
        bin_sum_qzz  = 0.0 
        bin_sum_qxy  = 0.0 
        bin_sum_qxz  = 0.0 
        bin_sum_qyz  = 0.0 
    end do  ! end loop over frames

    bin_aof      = sum(bin_ave,dim=2)/(ndx_record)
    bin_aof_px   = sum(bin_ave_px,dim=2)/(ndx_record)
    bin_aof_py   = sum(bin_ave_py,dim=2)/(ndx_record)
    bin_aof_pz   = sum(bin_ave_pz,dim=2)/(ndx_record)
    bin_aof_qxx  = sum(bin_ave_qxx,dim=2)/(ndx_record)
    bin_aof_qyy  = sum(bin_ave_qyy,dim=2)/(ndx_record)
    bin_aof_qzz  = sum(bin_ave_qzz,dim=2)/(ndx_record)
    bin_aof_qxy  = sum(bin_ave_qxy,dim=2)/(ndx_record)
    bin_aof_qxz  = sum(bin_ave_qxz,dim=2)/(ndx_record)
    bin_aof_qyz  = sum(bin_ave_qyz,dim=2)/(ndx_record)

    !bin_phi2_sum = sum(bin_phi2_ave,dim=2)/num_frame
    !bin_phi4_sum = sum(bin_phi4_ave,dim=2)/num_frame
    do bin_index = 1, num_bin
        write(21,'(f10.4,1x,f10.4,1x,&
            f10.4,1x,f10.4,1x,f10.4,&
            f10.4,1x,f10.4,1x,f10.4,&
            f10.4,1x,f10.4,1x,f10.4)') &
            bin_index*bin_size, bin_aof(bin_index), &
            bin_aof_px(bin_index) ,bin_aof_py(bin_index), bin_aof_pz(bin_index), &
            bin_aof_qxx(bin_index),bin_aof_qyy(bin_index),bin_aof_qzz(bin_index),&
            bin_aof_qxy(bin_index),bin_aof_qxz(bin_index),bin_aof_qyz(bin_index)
    end do

    !phi2_ave = sum(phi2_sum)/num_frame
    !phi4_ave = sum(phi4_sum)/num_frame
    !print*, 'phi2 ', sum(phi2_sum), phi2_ave
    !print*, 'phi4 ', sum(phi4_sum), phi4_ave

    close(10)
    close(11)
    close(21)

end program main
