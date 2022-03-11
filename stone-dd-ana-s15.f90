
module physics_constants
    implicit none
    real,    parameter :: e_charge       = 1.6021766208e-19        ! [C]
    real,    parameter :: Debye          = 0.208194                ! [e*A]
    real,    parameter :: pi             = 3.1415926
    real,    parameter :: B2A            = 0.529177249
    real,    parameter :: fg             = 2.9979307               ! [V*A^2/D]
    real,    parameter :: rho_w          = 0.033328493             ! [1/A^3]
    real,    parameter :: kT             = 0.02838                 ! [eV] at 330K
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

    real               :: O_locations(3,num_mole)
    real               :: H_locations(3,num_mole*2)
    integer            :: ndx_OH(num_mole*2)
    integer            :: num_seperate_H
    real               :: dist_OH(num_mole*2)

end module md_info

module component_info
    use md_info
    implicit none
    type mole_info                                            ! molecular info
        
        real    :: location(3,num_nuke)                       ! molecule nuclei location xyz 
        real    :: dipole(3)
        real    :: quadrupole(3,3)
        real    :: dipole_loc(3)
        real    :: quadrupole_loc(3,3)
        real    :: TrQ
        real    :: Pr, Qrr

        real    :: mass
        real    :: dist_to_cave
        real    :: T_vector(3)
        real    :: T_tensor(3,3)

        real    :: phi2 ! dipole     potential
        real    :: phi4 ! quadrupole potential

        real    :: R_vector(3)
        real    :: hoh_normal(3)
        real    :: ey1(3)

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
    
            integer :: i
            integer :: ndx_shift, ndx_nuke
            real    :: reference(3),  delta(3)
            real    :: r0, r1
            real    :: OH_vector(3)
            
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

            ! step three, place other nuke close to the reference to keep a whole molecule
            do ndx_nuke = 2, num_nuke
                OH_vector = nuke_locations(:, ndx_nuke) - nuke_locations(:,1)
                where (OH_vector > 0.5*L)
                    nuke_locations(:,ndx_nuke) = nuke_locations(:,ndx_nuke) - L
                else where (-0.5*L <= OH_vector .AND. OH_vector <= 0.5*L)
                    nuke_locations(:,ndx_nuke) = nuke_locations(:,ndx_nuke)
                else where 
                    nuke_locations(:,ndx_nuke) = nuke_locations(:,ndx_nuke) + L
                end where
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

        subroutine Rotate_from_lab_to_loc(ndx_mole)
            integer,intent(in)    :: ndx_mole
            real :: ex0(3), ey0(3), ez0(3)
            real :: ex1(3), ey1(3), ez1(3)
            real :: vec_oh1(3), vec_oh2(3), vec_oh0(3)
            real :: crossed(3)

            real :: dipole(3), dipole_loc(3)
            real :: quadrupole(3,3), quadrupole_0(3,3), quadrupole_loc(3,3)
            real :: Q_rot(3,3), Q_rot_T(3,3),   Q_rot_Inv(3,3), Q_rot_Inv_T(3,3)
            real :: Q_rot_0(3,3)

            ex0 = (/1, 0, 0/)
            ey0 = (/0, 1, 0/)
            ez0 = (/0, 0, 1/)

            dipole      = mole_list(ndx_mole)%dipole
            quadrupole  = mole_list(ndx_mole)%quadrupole

            vec_oh1 = mole_list(ndx_mole)%location(:,2)-mole_list(ndx_mole)%location(:,1)
            vec_oh2 = mole_list(ndx_mole)%location(:,3)-mole_list(ndx_mole)%location(:,1)

            !if ( norm2(vec_oh1) >= norm2(vec_oh2)+0.000001 ) then
            !    vec_oh0 = vec_oh1
            !    vec_oh1 = vec_oh2
            !    vec_oh2 = vec_oh0
            !end if
            !
            !ex1         = vec_oh1/norm2(vec_oh1)
            !crossed(1)  = vec_oh1(2)*vec_oh2(3) - vec_oh1(3)*vec_oh2(2)
            !crossed(2)  = vec_oh1(3)*vec_oh2(1) - vec_oh1(1)*vec_oh2(3)
            !crossed(3)  = vec_oh1(1)*vec_oh2(2) - vec_oh1(2)*vec_oh2(1)
            !ez1         = crossed/norm2(crossed)
            !crossed(1)  = ez1(2)*ex1(3)-ez1(3)*ex1(2)
            !crossed(2)  = ez1(3)*ex1(1)-ez1(1)*ex1(3)
            !crossed(3)  = ez1(1)*ex1(2)-ez1(2)*ex1(1)
            !ey1         = crossed/(NORM2(crossed))

            vec_oh0     = 0.5*(vec_oh1+vec_oh2)
            ex1         = vec_oh0/norm2(vec_oh0)

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

            call stone_inverse(Q_rot_0,Q_rot_Inv)
            Q_rot_Inv_T = transpose(Q_rot_Inv)

            dipole_loc(1)   = dot_product(dipole,Q_rot_Inv(:,1))
            dipole_loc(2)   = dot_product(dipole,Q_rot_Inv(:,2))
            dipole_loc(3)   = dot_product(dipole,Q_rot_Inv(:,3))

            quadrupole_0    = matmul(Q_rot_Inv_T,quadrupole)
            quadrupole_loc  = matmul(quadrupole_0,Q_rot_Inv)

            mole_list(ndx_mole)%dipole_loc     = dipole_loc
            mole_list(ndx_mole)%quadrupole_loc = quadrupole_loc
            mole_list(ndx_mole)%hoh_normal     = ez1
            mole_list(ndx_mole)%ey1            = ey1  

        end subroutine Rotate_from_lab_to_loc

        subroutine transform_cartesian_2_spherical(ndx_mole)

            integer, intent(in) :: ndx_mole
            real :: x, y, z, r, theta, psi
            real :: dist_vec(3)
            real :: cos_theta,sin_theta
            real :: cos_psi,  sin_psi
            real :: c1, s1, c2, s2
            real :: tensor(3,3), vector(3)
            real :: Qrr, Qtt, Qpp, Qrt, Qrp, Qtp
            real :: Pr, Pt, Pp

            dist_vec = mole_list(ndx_mole)%location(:,1)-cave_location ! [Angstrom]
            vector   = mole_list(ndx_mole)%dipole
            tensor   = mole_list(ndx_mole)%quadrupole
            x = dist_vec(1)
            y = dist_vec(2)
            z = dist_vec(3)
            r = norm2(dist_vec)
            cos_theta = z/r
            sin_theta = (1-cos_theta**2)**0.5
            if (sin_theta .lt. 1e-008) then
                cos_psi   = 0.0
                sin_psi   = 0.0
            else
                cos_psi   = x/(r*sin_theta)
                sin_psi   = y/(r*sin_theta)
            end if
            c1 = cos_theta
            s1 = sin_theta
            c2 = cos_psi
            s2 = sin_psi
            Qrr = tensor(1,1)*s1**2*c2**2+tensor(2,2)*s1**2*s2**2+tensor(3,3)*c1**2+&
                 (tensor(1,2)*s1*s2*c2   +tensor(1,3)*c1*c2      +tensor(2,3)*2*c1*s2)*2*s1
            Pr  = vector(1)*s1*c2+vector(2)*s1*s2+vector(3)*c1
            !Pr  = dot_product(vector,dist_vec/r) ! equal to the above line
            !print*, ndx_mole, r, c1, s1, c2, s2, cos_theta, sin_theta, Pr, Qrr
            mole_list(ndx_mole)%Pr  = Pr
            mole_list(ndx_mole)%Qrr = Qrr

        end subroutine transform_cartesian_2_spherical

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

            mole_list(ndx_mole)%R_vector = dist_vector/dist_R

        end subroutine comp_Ttensor

end module compute


program main
    use physics_constants
    use md_info
    use component_info
    use tools
    use compute
    implicit none

    integer :: ndx_frame, ndx_mole, ndx_record, num_record
    integer :: buf_int(50)
    real    :: buf_real(50)
    character(len=100) :: buf_str0, buf_str(50)

    integer, parameter :: num_bin    = 400
    real    :: rmax                  = 6.50
    real    :: bin_size              = 0.10
    real    :: rmin,dr,r
    integer :: ndx_rmin, ndx_rmax, ndx_bin
    real    :: dV

    real,    dimension(:),   allocatable :: bin_Pr, bin_Qrr, bin_TrQ
    real,    dimension(:),   allocatable :: bulk_Qrr

    real,    dimension(:,:), allocatable :: delta_phi_D
    real,    dimension(:,:), allocatable :: delta_phi_Q2
    real,    dimension(:,:), allocatable :: delta_phi_Q1
    real,    dimension(:,:), allocatable :: delta_phi

    real,    dimension(:),   allocatable :: delta_phi_QB
    real,    dimension(:),   allocatable :: delta_phi_B
    real,    dimension(:),   allocatable :: delta_phi_DQ2
    real,    dimension(:),   allocatable :: delta_phi_fluc2
    real,    dimension(:),   allocatable :: delta_phi_fluc3

    real,    dimension(:,:), allocatable :: angle_dip
    real,    dimension(:,:), allocatable :: angle_hoh
    real,    dimension(:,:), allocatable :: angle_ey1
    real,    dimension(:,:), allocatable :: project_dip

    integer, dimension(:,:), allocatable :: bin_count


    real    :: beta = 1/kT
    real    :: cumulant1, cumulant2, cumulant3
    real    :: raw_moment1, raw_moment2, raw_moment3

    real    :: angle
    real    :: TrQ_rmax
    real    :: rho_rmax

    character(100) :: frame_ndx0, frame_ndx1

    call GET_COMMAND_ARGUMENT(1,frame_ndx0)
    call GET_COMMAND_ARGUMENT(2,frame_ndx1)
    read(frame_ndx0,*) ndx_frame0
    read(frame_ndx1,*) ndx_frame1

    open(unit=10,file='water.xyz',      status='old',action='read')
    open(unit=11,file='dipole',         status='old',action='read')
    open(unit=12,file='quadrupole',     status='old',action='read')

    open(unit=20,file='phi_dist',       status='unknown',action='write')
    open(unit=21,file='delta_phi',      status='unknown',action='write')
    open(unit=22,file='dipole_loc',     status='unknown',action='write')
    open(unit=23,file='quadrupole_loc', status='unknown',action='write')
    open(unit=24,file='orientation',    status='unknown',action='write')
    
    num_record = ndx_frame1-ndx_frame0+1
   ! initialization
    allocate(bin_count(num_bin,num_record))

    allocate(bin_Pr(num_bin))
    allocate(bin_Qrr(num_bin))
    allocate(bin_TrQ(num_bin))

    allocate(bulk_Qrr(num_record))

    allocate(delta_phi_D(num_bin,num_record))
    allocate(delta_phi_Q2(num_bin,num_record))
    allocate(delta_phi_Q1(num_bin,num_record))
    allocate(delta_phi(num_bin,num_record))
    allocate(delta_phi_QB(num_record))
    allocate(delta_phi_B(num_record))
    allocate(delta_phi_DQ2(num_record))

    allocate(delta_phi_fluc2(num_record))
    allocate(delta_phi_fluc3(num_record))

    allocate(angle_dip(num_bin,num_record))
    allocate(angle_hoh(num_bin,num_record))
    allocate(angle_ey1(num_bin,num_record))
    allocate(project_dip(num_bin,num_record))
    
    bin_count       = 0.0
    
    bin_Pr          = 0.0
    bin_Qrr         = 0.0
    bin_TrQ         = 0.0

    bulk_Qrr        = 0.0

    delta_phi_D     = 0.0
    delta_phi_Q2    = 0.0
    delta_phi_Q1    = 0.0
    delta_phi       = 0.0

    delta_phi_QB    = 0.0
    delta_phi_B     = 0.0
    delta_phi_DQ2   = 0.0
    delta_phi_fluc2 = 0.0
    delta_phi_fluc3 = 0.0

    angle_dip       = 0.0
    angle_hoh       = 0.0
    angle_ey1       = 0.0

    dr              = bin_size
    rmin            = bin_size
    dV              = 4*pi*rmax**2*dr

    TrQ_rmax        = 0.0
    rho_rmax        = 0.0

    ndx_rmin        = floor(rmin/bin_size)
    ndx_rmax        = floor(rmax/bin_size)

   ! skip frames to the frame ndx_frame0
    do ndx_record = 1, ndx_frame0-1 
        print*, 'skip lines', ndx_record, ndx_frame0-1
        read(10,*) buf_str(1)
        read(10,*) buf_str(1)
        do ndx_mole = 1, num_mole
            read(10,*) buf_str0, mole_list(ndx_mole)%location(:,1)
            read(10,*) buf_str0, mole_list(ndx_mole)%location(:,2)
            read(10,*) buf_str0, mole_list(ndx_mole)%location(:,3)
            read(11,*) buf_int(1:2), mole_list(ndx_mole)%dipole
            read(12,*) buf_int(1:2), buf_real(1:9)
        end do
        read(10,*) buf_str(1)
    end do

    do ndx_record = ndx_frame0, ndx_frame1
        ndx_frame = ndx_record - ndx_frame0 + 1
        bin_Pr    = 0.0
        bin_Qrr   = 0.0
        bin_TrQ   = 0.0

        read(10,*) buf_str(1)
        read(10,*) buf_str(1)
        do ndx_mole = 1, num_mole
            read(10,*) buf_str0, mole_list(ndx_mole)%location(:,1)
            read(10,*) buf_str0, mole_list(ndx_mole)%location(:,2)
            read(10,*) buf_str0, mole_list(ndx_mole)%location(:,3)
            read(11,*) buf_int(1:2), mole_list(ndx_mole)%dipole
            read(12,*) buf_int(1:2), buf_real(1:6)
            mole_list(ndx_mole)%quadrupole(1,1) = buf_real(1)
            mole_list(ndx_mole)%quadrupole(2,1) = buf_real(2)
            mole_list(ndx_mole)%quadrupole(3,1) = buf_real(3)
            mole_list(ndx_mole)%quadrupole(1,2) = buf_real(2)
            mole_list(ndx_mole)%quadrupole(2,2) = buf_real(4)
            mole_list(ndx_mole)%quadrupole(3,2) = buf_real(5)
            mole_list(ndx_mole)%quadrupole(1,3) = buf_real(3)
            mole_list(ndx_mole)%quadrupole(2,3) = buf_real(5)
            mole_list(ndx_mole)%quadrupole(3,3) = buf_real(6)
            mole_list(ndx_mole)%TrQ = mole_list(ndx_mole)%quadrupole(1,1) +&
                                      mole_list(ndx_mole)%quadrupole(2,2) +&
                                      mole_list(ndx_mole)%quadrupole(3,3)
        end do
        read(10,*) buf_str0

        do ndx_mole = 1, num_mole
            call pbc_mole(mole_list(ndx_mole)%location)
            call Rotate_from_lab_to_loc(ndx_mole)
            call comp_Ttensor(ndx_mole)
            call transform_cartesian_2_spherical(ndx_mole)

            write(22, '(I8,I5,4f10.5)') ndx_frame, ndx_mole, &
            mole_list(ndx_mole)%dipole_loc, norm2(mole_list(ndx_mole)%dipole_loc)
            write(23, '(I8,I5,9f10.5)') ndx_frame, ndx_mole, &
    mole_list(ndx_mole)%quadrupole_loc(1,1), mole_list(ndx_mole)%quadrupole_loc(2,1), mole_list(ndx_mole)%quadrupole_loc(3,1), &
    mole_list(ndx_mole)%quadrupole_loc(1,2), mole_list(ndx_mole)%quadrupole_loc(2,2), mole_list(ndx_mole)%quadrupole_loc(3,2), &
    mole_list(ndx_mole)%quadrupole_loc(1,3), mole_list(ndx_mole)%quadrupole_loc(2,3), mole_list(ndx_mole)%quadrupole_loc(3,3)

            ndx_bin = floor(mole_list(ndx_mole)%dist_to_cave/bin_size)
            if ( ndx_bin .le. num_bin ) then
                bin_count(ndx_bin,ndx_frame) = bin_count(ndx_bin,ndx_frame)  + 1

                bin_Pr(ndx_bin)  = bin_Pr(ndx_bin) + mole_list(ndx_mole)%Pr
                bin_Qrr(ndx_bin) = bin_Qrr(ndx_bin)+ mole_list(ndx_mole)%Qrr
                bin_TrQ(ndx_bin) = bin_TrQ(ndx_bin)+ mole_list(ndx_mole)%TrQ

                angle = ACOS(dot_product(mole_list(ndx_mole)%R_vector, &
                    mole_list(ndx_mole)%dipole)/norm2(mole_list(ndx_mole)%dipole))/pi*180
                angle_dip(ndx_bin,ndx_frame) = angle_dip(ndx_bin,ndx_frame) + angle
                angle = ACOS(dot_product(mole_list(ndx_mole)%R_vector, &
                    mole_list(ndx_mole)%hoh_normal))/pi*180
                angle_hoh(ndx_bin,ndx_frame) = angle_hoh(ndx_bin,ndx_frame) + angle
                angle = ACOS(dot_product(mole_list(ndx_mole)%R_vector, &
                    mole_list(ndx_mole)%ey1))/pi*180
                angle_ey1(ndx_bin,ndx_frame) = angle_ey1(ndx_bin,ndx_frame) + angle
                project_dip(ndx_bin,ndx_frame) =  project_dip(ndx_bin,ndx_frame)+ mole_list(ndx_mole)%Pr
            end if

            if (ndx_bin .ge. ndx_rmax) then
                bulk_Qrr(ndx_frame) = bulk_Qrr(ndx_frame) + mole_list(ndx_mole)%Qrr
            end if

        end do  ! end loop over molecule index

        if (bin_count(ndx_rmax,ndx_frame) .ne. 0) then
            TrQ_rmax = TrQ_rmax + bin_TrQ(ndx_rmax)/bin_count(ndx_rmax,ndx_frame)
            !print*, 'TrQ_rmax = ', TrQ_rmax
        end if

        do ndx_bin = ndx_rmin, ndx_rmax
            r                               = bin_size*ndx_bin
            delta_phi_D(ndx_bin,ndx_frame)  = fg*bin_Pr(ndx_bin)/r**2
            delta_phi_Q2(ndx_bin,ndx_frame) = fg*(bin_TrQ(ndx_bin)-3*bin_Qrr(ndx_bin))/r**3
            delta_phi_Q1(ndx_bin,ndx_frame) = -1*fg*bin_Qrr(ndx_bin)*4*pi/(4*pi*r**2*bin_size)
            delta_phi(ndx_bin,ndx_frame)    = delta_phi_D(ndx_bin,ndx_frame)  + &
                                              delta_phi_Q2(ndx_bin,ndx_frame) + &
                                              delta_phi_Q1(ndx_bin,ndx_frame)

            if (bin_count(ndx_bin,ndx_frame) .ne. 0) then
                angle_dip(ndx_bin,ndx_frame)   = angle_dip(ndx_bin,ndx_frame)/bin_count(ndx_bin,ndx_frame)
                angle_hoh(ndx_bin,ndx_frame)   = angle_hoh(ndx_bin,ndx_frame)/bin_count(ndx_bin,ndx_frame)
                angle_ey1(ndx_bin,ndx_frame)   = angle_ey1(ndx_bin,ndx_frame)/bin_count(ndx_bin,ndx_frame)
                project_dip(ndx_bin,ndx_frame) = project_dip(ndx_bin,ndx_frame)/bin_count(ndx_bin,ndx_frame)

            else
                angle_dip(ndx_bin,ndx_frame)   = -10.0
                angle_hoh(ndx_bin,ndx_frame)   = -10.0
                angle_ey1(ndx_bin,ndx_frame)   = -10.0
                project_dip(ndx_bin,ndx_frame) = -100.0
            end if
        end do
        delta_phi_QB(ndx_frame) = -fg*bulk_Qrr(ndx_frame)*(4*pi)/(L**3-4*pi/3*((ndx_rmax-1)*bin_size)**3)
        delta_phi_B(ndx_frame)  = sum(delta_phi_D(ndx_rmin:ndx_rmax,ndx_frame))  &
                                + sum(delta_phi_Q2(ndx_rmin:ndx_rmax,ndx_frame)) &
                                + delta_phi_QB(ndx_frame)
        delta_phi_DQ2(ndx_frame)= sum(delta_phi_D(ndx_rmin:ndx_rmax,ndx_frame))  &
                                + sum(delta_phi_Q2(ndx_rmin:ndx_rmax,ndx_frame))

        write(24,'(I10, 2x, 8(f10.4,2x))') ndx_frame, angle_dip(41,ndx_frame), angle_dip(65,ndx_frame),&
                                                      angle_hoh(41,ndx_frame), angle_hoh(65,ndx_frame),&
                                                      angle_ey1(41,ndx_frame), angle_ey1(65,ndx_frame),&
                                                      project_dip(41,ndx_frame), project_dip(65,ndx_frame)

        print*, 'frame ', ndx_record
        print*, '      '


    end do  ! end loop over frames

    !print*,  '4.1, orientation angle to dip: ', sum(angle_dip(41,:))/num_record
    !print*,  '4.1, orientation angle to hoh: ', sum(angle_hoh(41,:))/num_record
    !9print*,  '6.5, orientation angle to dip: ', sum(angle_dip(65,:))/num_record
    !print*,  '6.5, orientation angle to hoh: ', sum(angle_hoh(65,:))/num_record

    do ndx_bin = ndx_rmin, ndx_rmax
        write(21,'(f5.2, 2x, 4(f10.4,2x))') ndx_bin*bin_size, &
            sum(delta_phi_D(ndx_rmin:ndx_bin,:))/num_record,           &
            sum(delta_phi_Q2(ndx_rmin:ndx_bin,:))/num_record,          &
            sum(delta_phi_Q1(ndx_bin,:))/num_record,          &
            sum(delta_phi_D(ndx_rmin:ndx_bin,:))/num_record + sum(delta_phi_Q2(ndx_rmin:ndx_bin,:))/num_record + &
            sum(delta_phi_Q1(ndx_bin,:))/num_record
    end do
    rho_rmax = float(sum(bin_count(ndx_rmax,:)))/float(num_record)/dV*1000

    print('(A,f10.2)'), 'bin_size  = ',    bin_size
    print('(A,f10.2)'), 'rmax      = ',    rmax
    print('(A,f10.2)'), 'delta phi  D',    0.01*NINT(sum(delta_phi_D(ndx_rmin:ndx_rmax,:))/num_record*100)
    print('(A,f10.2)'), 'delta phi Q2',    0.01*NINT(sum(delta_phi_Q2(ndx_rmin:ndx_rmax,:))/num_record*100)
    print('(A,f10.2)'), 'delta phi Q1',    0.01*NINT(sum(delta_phi_QB)/num_record*100)
    print('(A,f10.2)'), 'delta phi  =',    0.01*NINT(sum(delta_phi_B)/num_record*100)
    print*, '      '
    print'(A,f10.2)',   'rho(rmax) = ',    rho_rmax
    print'(A,f10.2)',   'TrQ(rmax) = ',    0.01*NINT(TrQ_rmax/num_record*100)
    print('(A,f10.2)'), 'bethe pot = ',    0.01*NINT(-fg*4*pi/3*TrQ_rmax/num_record*rho_rmax*100/1000)

    raw_moment1 = sum(delta_phi_B)/num_record
    do ndx_frame = 1, num_record
        delta_phi_fluc2(ndx_frame) = (delta_phi_B(ndx_frame)-raw_moment1)**2
        delta_phi_fluc3(ndx_frame) = (delta_phi_B(ndx_frame)-raw_moment1)**3
        write(20,'(I10,2x,f10.4,2x,f10.4)') ndx_frame, delta_phi_B(ndx_frame),delta_phi_DQ2(ndx_frame)
    end do
   
    raw_moment2 = sum(delta_phi_fluc2)/num_record
    raw_moment3 = sum(delta_phi_fluc3)/num_record

    cumulant1   = raw_moment1
    cumulant2   = raw_moment2
    cumulant3   = raw_moment3

    print*, 'cumulant2 = ', cumulant2
    print*, 'cumulant3 = ', cumulant3

    print*, '   '
    print'(A,f10.2)', 'cumulant1 = ', 0.01*NINT(cumulant1*100)
    print'(A,f10.2)', 'cumulant3 = ', 0.01*NINT(1./6.0*beta**2*cumulant3*100)
    
    close(10)
    close(11)
    close(12)

    close(20)
    close(21)   
    close(22)
    close(23)
    close(24)

end program main
