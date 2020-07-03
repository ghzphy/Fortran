module para
    implicit none
    real*8, parameter :: omega = 0.3d0
    real*8, parameter :: eta = 0.01d0
    real*8, parameter :: pi = 3.14159265358979d0
    real*8, parameter :: V0 = 0.05d0
    real*8, parameter :: c = 500d0 ! spin-orbit coupling parameter
    real*8, parameter :: qx_min = -0.314d0
    real*8, parameter :: qx_max = -qx_min
    real*8, parameter :: qy_min = -0.314d0
    real*8, parameter :: qy_max = -qy_min
    real*8, parameter :: kx_min = -0.628d0
    real*8, parameter :: kx_max = -kx_min
    real*8, parameter :: ky_min = -0.628d0
    real*8, parameter :: ky_max = -ky_min
    integer*4, parameter :: n_qx = 100
    integer*4, parameter :: n_qy = n_qx
    integer*4, parameter :: n_kx = 200
    integer*4, parameter :: n_ky = n_kx
    integer*4, parameter :: Nq = n_qx*n_qy
    integer*4, parameter :: Nk = n_kx*n_ky
    real*8, parameter :: delta_qx = (qx_max - qx_min)/n_qx
    real*8, parameter :: delta_qy = (qy_max - qy_min)/n_qy
    real*8, parameter :: delta_kx = (kx_max - kx_min)/n_kx
    real*8, parameter :: delta_ky = (ky_max - ky_min)/n_ky
    complex*16, parameter, dimension(2,2) :: sigma_x = transpose(reshape((/0.,1.,1.,0./), shape(sigma_x)))
    complex*16, parameter, dimension(2,2) :: sigma_y = transpose(reshape((/(0.,0.),(0.,-1.),(0.,1.),(0.,0.)/), shape(sigma_y)))
    complex*16, parameter, dimension(2,2) :: sigma_z = transpose(reshape((/1.,0.,0.,-1./), shape(sigma_z)))
    complex*16, parameter, dimension(2,2) :: sigma_0 = transpose(reshape((/1.,0.,0.,1./), shape(sigma_0)))

    contains
    !
    ! The function : H(k) = \sum_{k} d(k)*sigma_k
    !
    function dk(kx,ky)
        implicit none
        real*8 :: kx, ky
        real*8, dimension(0:3) :: dk
        real*8, parameter :: v = 2.55d0
        real*8, parameter :: lambda = 255d0
        dk(0) = 0d0
        dk(1) = -v*ky
        dk(2) = v*kx
        dk(3) = lambda*(-3*kx**2*ky+ky**3)
    end function
    !
    ! Get the inverse matrix.
    !
    function matinv(m)
        complex*16, dimension(2,2) :: m(2,2)
        complex*16, dimension(2,2) :: matinv(2,2)
        complex*16 :: det
        det = 1./(m(1,1)*m(2,2) - m(1,2)*m(2,1))
        matinv(1,1) = det*m(2,2)
        matinv(2,2) = det*m(1,1)
        matinv(1,2) = -det*m(1,2)
        matinv(2,1) = -det*m(2,1)
    end function
end module

program main
    use para
    implicit none
    !use mpi
    include 'mpif.h'
    
    integer*4 num_cpu, cpuid, ierr
    real*8 wtime_diff,wtime_end,wtime_start
    
    !
    ! Initialize MPI.
    !
    call MPI_INIT(ierr)
    !
    ! Get this process's ID
    !
    call MPI_COMM_RANK(MPI_COMM_WORLD,cpuid,ierr)
    !
    ! Find out how many processes are available.
    !
    call MPI_COMM_SIZE(MPI_COMM_WORLD,num_cpu,ierr)
    
    if (cpuid == 0 ) then
        write(*,'(a)') '================     QPI_MPI version     ================'
        write(*,'(a,i2,a)') '=> You are using ' ,num_cpu, ' CPUs.'
        call timestamp()
        wtime_start = MPI_Wtime()
    end if
    
    call T_matrix_QPI(cpuid,num_cpu)
    !call JDOS_QPI(cpuid,num_cpu)
   ! call Spin_Orbit_QPI(cpuid,num_cpu)
 
    if (cpuid == 0) then
        wtime_end = MPI_Wtime()
        wtime_diff = wtime_end - wtime_start
        
        write(*,'(a)') ' '
        write(*,'(a,g14.6,a)') '   Elapsed wall clock seconds = ',wtime_diff,' s.'
        write(*,'(a)') ' '
        write(*,'(a)') '=> The QPI calculation was done successfully!'
        call timestamp()
        write(*,'(a)') '========================================================='
    end if
    !
    ! Terminate the MPI.
    !
    call MPI_FINALIZE(ierr)
end program


subroutine T_matrix_QPI(cpuid,num_cpu)
    use para
    implicit none
    include 'mpif.h'
    
    integer*4 num_cpu, cpuid, ierr
    complex*16, allocatable, dimension(:,:,:,:) :: G0k_matrix,G0kq_matrix
    complex*16, dimension(2,2) :: V_matrix,int_G0k,T_matrix,G_matrix
    real*8, dimension(0:3) :: dk_tmp,dk_tmp1
    integer*4 :: ikx,iky,iqx,iqy,ikpx,ikpy,ikp,i,j,iq
    real*8 :: kx,ky,kpx,kpy,qx_min_shape,qx_max_shape,qy_min_shape,qy_max_shape
    real*8, allocatable :: rho(:),rho_part(:)
    integer*4, allocatable :: q12(:,:)
    real*8 wtime_diff,wtime_end,wtime_start
    
    allocate(G0k_matrix(2,2,1:n_kx,1:n_ky))
    allocate(G0kq_matrix(2,2,1:n_kx+n_qx,1:n_ky+n_qy))
    allocate(q12(2,1:Nq))
    allocate(rho(1:Nq),rho_part(1:Nq))
    
    if (cpuid == 0 ) then
        write(*,'(a)') '=> The integral of T_matrix QPI start:'
        do ikx = 1,n_kx
            do iky = 1,n_ky
                kx = (ikx - 1.0/2)*delta_kx + kx_min
                ky = (iky - 1.0/2)*delta_ky + ky_min
                dk_tmp = dk(kx,ky)
                G0k_matrix(:,:,ikx,iky) = matinv(complex(omega,eta)*sigma_0 - &
                & (dk_tmp(0)*sigma_0 + dk_tmp(1)*sigma_x + dk_tmp(2)*sigma_y + dk_tmp(3)*sigma_z))
            end do
        end do
    
        int_G0k(1,1) = sum(G0k_matrix(1,1,:,:))
        int_G0k(1,2) = sum(G0k_matrix(1,2,:,:))
        int_G0k(2,1) = sum(G0k_matrix(2,1,:,:))
        int_G0k(2,2) = sum(G0k_matrix(2,2,:,:))
        int_G0k = int_G0k*delta_kx*delta_ky
        
        V_matrix = 0d0
        V_matrix(1,1) = V0
        V_matrix(2,2) = V0
        !
        ! For scalar delta potential, T_matrix can be solved exactly.
        ! T = V/(I - V*sum G0)
        !
        T_matrix = matmul(matinv(sigma_0 - matmul(V_matrix,int_G0k)),V_matrix)
    end if
    !write(*,'(a,i8,a)') 'Process ', cpuid,'  is active.'

    call MPI_Bcast(T_matrix,size(T_matrix),MPI_double_complex,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(G0k_matrix,size(G0k_matrix),MPI_double_complex,0,MPI_COMM_WORLD,ierr)
    
    iq = 0
    do i = 1,n_qx
        do j = 1,n_qy
            iq = iq + 1
            q12(1,iq) = i
            q12(2,iq) = j
        end do
    end do
    
    rho_part = 0d0
    rho = 0d0
    
    do iq = cpuid + 1,Nq, num_cpu
        iqx = q12(1,iq)
        iqy = q12(2,iq)
        do ikx = 1,n_kx
            do iky = 1,n_ky
                ikpx = ikx + iqx
                ikpy = iky + iqy
                kpx = (ikx - 1.0/2)*delta_kx + kx_min + (iqx - 1.0/2)*delta_qx + qx_min
                kpy = (iky - 1.0/2)*delta_ky + ky_min + (iqy - 1.0/2)*delta_qy + qy_min
                dk_tmp1 = dk(kpx,kpy)
                G0kq_matrix(:,:,ikpx,ikpy) = matinv(complex(omega,eta)*sigma_0- &
                & (dk_tmp1(0)*sigma_0+dk_tmp1(1)*sigma_x+dk_tmp1(2)*sigma_y+dk_tmp1(3)*sigma_z))
                G_matrix = matmul(matmul(G0k_matrix(:,:,ikx,iky),T_matrix),G0kq_matrix(:,:,ikpx,ikpy))
                rho_part(iq)= rho_part(iq) + imagpart(G_matrix(1,1)+G_matrix(2,2))
            end do
        end do
        rho_part(iq) = rho_part(iq)*delta_kx*delta_ky
    end do

    call MPI_REDUCE(rho_part,rho,size(rho_part),MPI_real8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    
    if (cpuid == 0) then
        write(*,'(a)') '   the integral was done.'
        write(*,'(a)') "=> Now it's ready to write qpi_Tm.dat and qpi_Tm.gnu>."
        open(10,file = 'qpi_Tm.dat')
        do iq = 1, Nq
            if (iq==1) then
                qx_min_shape = qx_min+(q12(1,iq)-1.0/2)*delta_qx
                qy_min_shape = qy_min+(q12(2,iq)-1.0/2)*delta_qy
            end if
            if (iq==Nq) then
                qx_max_shape = qx_min+(q12(1,iq)-1.0/2)*delta_qx
                qy_max_shape = qy_min+(q12(2,iq)-1.0/2)*delta_qy
            end if
            write(10,*) qx_min+(q12(1,iq)-1.0/2)*delta_qx, qy_min+(q12(2,iq)-1.0/2)*delta_qy, rho(iq)
            if (mod(iq, n_qy)==0) write(10, *)' '
        end do
        close(10)

        open(11,file = 'qpi_Tm.gnu')
        write(11,'(a)') 'set encoding iso_8859_1'
        write(11,'(a)') 'set terminal png truecolor enhanced font ",50" size 1920, 1680'
        write(11,'(a)') "set output 'qpi_Tm.png'"
        write(11,'(a)') 'unset ztics'
        write(11,'(a)') 'unset key'
        write(11,'(a)') 'set pm3d'
        write(11,'(a)') 'set border lw 6'
        write(11,'(a)') 'set size ratio -1'
        write(11,'(a)') 'set view map'
        write(11,'(a)') 'set xtics'
        write(11,'(a)') 'set ytics'
        write(11,'(a)') 'set xlabel "q_x"'
        write(11,'(a)') 'set ylabel "q_y"'
        write(11,'(a)') 'set ylabel offset 1,0'
        write(11,'(a)') 'set colorbox'
        write(11,'(a)') 'unset cbtics'
        write(11,'(a,f8.6,a,f8.6,a)') 'set xrange [', qx_min_shape, ':', qx_max_shape, ']'
        write(11,'(a,f8.6,a,f8.6,a)') 'set yrange [', qy_min_shape, ':', qy_max_shape, ']'
        write(11,'(a)') 'set pm3d interpolate 2,2'
        write(11,'(a)') "splot 'qpi_Tm.dat' u 1:2:3 w pm3d"
        close(11)
    end if

    !call mpi_barrier(MPI_COMM_WOLRD, ierr)

    deallocate(G0k_matrix)
    deallocate(G0kq_matrix)
    deallocate(q12)
    deallocate(rho,rho_part)

    return
end subroutine

subroutine JDOS_QPI(cpuid,num_cpu)
    !
    ! The subroutine include two calculations : JDODS and JSDOS.
    ! Modified by Ghz. 2019.9.30
    ! Joint Density of states: JDOS
    ! Joint Spin Density of states: JSDOS
    ! \rho(q,w) = \sum_{i=0,x,y,z} \sum_{k} \rho_{i}(k,w)*\rho_{i}(k+q,w)
    ! \rho_{i}(k,w) = -Im[Tr{\sigma_{i}*G0(k,w)}],
    ! \rho_{i}(k+q,w) = -Im[Tr{\sigma_{i}*G0(k+q,w)}]
    ! \rho_{0}(k,w) : rhok0, \rho_{1}(k,w) : rhok1
    ! \rho_{2}(k,w) : rhok2, \rho_{3}(k,w) : rhok3
    ! \rho_{0}(k+q,w) : rhokq0, \rho_{1}(k+q,w) : rhokq1
    ! \rho_{2}(k+q,w) : rhokq2, \rho_{3}(k+q,w) : rhokq3
    !
    use para
    implicit none
    !use mpi
    include 'mpif.h'

    integer*4 num_cpu, cpuid, ierr
    complex*16, allocatable, dimension(:,:,:,:) :: G0k_matrix,G0kq_matrix
    complex*16, dimension(2,2) :: rhok0,rhok1,rhok2,rhok3,rhokq0,rhokq1,rhokq2,rhokq3
    real*8, dimension(0:3) :: dk_tmp,dk_tmp1
    integer*4 :: ikx,iky,iqx,iqy,ikpx,ikpy,ikp,i,j,iq
    real*8 :: kx,ky,kpx,kpy,qx_min_shape,qx_max_shape,qy_min_shape,qy_max_shape
    real*8, allocatable :: rho_jdos(:),rho_jdos_part(:),rho_jsdos(:),rho_jsdos_part(:)
    real*8 :: tmp0,tmp1,tmp2,tmp3
    integer*4, allocatable :: q12(:,:)

    allocate(G0k_matrix(2,2,1:n_kx,1:n_ky))
    allocate(G0kq_matrix(2,2,1:n_kx+n_qx,1:n_ky+n_qy))
    allocate(q12(2,1:Nq))
    allocate(rho_jdos(1:Nq),rho_jsdos(1:Nq),rho_jdos_part(1:Nq),rho_jsdos_part(1:Nq))

    if (cpuid == 0 ) then
        write(*,'(a)') '=> The integral of JDOS/JSDOS QPI start:'
        do ikx = 1,n_kx
            do iky = 1,n_ky
                kx = (ikx - 1.0/2)*delta_kx + kx_min
                ky = (iky - 1.0/2)*delta_ky + ky_min
                dk_tmp = dk(kx,ky)
                G0k_matrix(:,:,ikx,iky) = matinv(complex(omega,eta)*sigma_0 - &
                & (dk_tmp(0)*sigma_0 + dk_tmp(1)*sigma_x + dk_tmp(2)*sigma_y + dk_tmp(3)*sigma_z))
            end do
        end do
    end if
    call MPI_Bcast(G0k_matrix,size(G0k_matrix),MPI_double_complex,0,MPI_COMM_WORLD,ierr)

    iq = 0
    do i = 1,n_qx
        do j = 1,n_qy
            iq = iq + 1
            q12(1,iq) = i
            q12(2,iq) = j
        end do
    end do

    rho_jdos_part = 0d0
    rho_jdos = 0d0
    rho_jsdos_part = 0d0
    rho_jsdos = 0d0

    do iq = cpuid + 1, Nq, num_cpu
        iqx = q12(1,iq)
        iqy = q12(2,iq)
        do ikx = 1,n_kx
            do iky = 1,n_ky
                ikpx = ikx + iqx
                ikpy = iky + iqy
                kpx = (ikx - 1.0/2)*delta_kx + kx_min + (iqx - 1.0/2)*delta_qx + qx_min
                kpy = (iky - 1.0/2)*delta_ky + ky_min + (iqy - 1.0/2)*delta_qy + qy_min
                dk_tmp1 = dk(kpx,kpy)
                G0kq_matrix(:,:,ikpx,ikpy) = matinv(complex(omega,eta)*sigma_0-&
                & (dk_tmp1(0)*sigma_0+dk_tmp1(1)*sigma_x+dk_tmp1(2)*sigma_y+dk_tmp1(3)*sigma_z))

                rhok0 = matmul(sigma_0,G0k_matrix(:,:,ikx,iky))
                rhok1 = matmul(sigma_x,G0k_matrix(:,:,ikx,iky))
                rhok2 = matmul(sigma_y,G0k_matrix(:,:,ikx,iky))
                rhok3 = matmul(sigma_z,G0k_matrix(:,:,ikx,iky))

                rhokq0 = matmul(sigma_0,G0kq_matrix(:,:,ikpx,ikpy))
                rhokq1 = matmul(sigma_x,G0kq_matrix(:,:,ikpx,ikpy))
                rhokq2 = matmul(sigma_y,G0kq_matrix(:,:,ikpx,ikpy))
                rhokq3 = matmul(sigma_z,G0kq_matrix(:,:,ikpx,ikpy))

                tmp0 = imagpart(rhok0(1,1)+rhok0(2,2))*imagpart(rhokq0(1,1)+rhokq0(2,2))
                tmp1 = imagpart(rhok1(1,1)+rhok1(2,2))*imagpart(rhokq1(1,1)+rhokq1(2,2))
                tmp2 = imagpart(rhok2(1,1)+rhok2(2,2))*imagpart(rhokq2(1,1)+rhokq2(2,2))
                tmp3 = imagpart(rhok3(1,1)+rhok3(2,2))*imagpart(rhokq3(1,1)+rhokq3(2,2))

                rho_jdos_part(iq)= rho_jdos_part(iq) + tmp0
                rho_jsdos_part(iq)= rho_jsdos_part(iq) + tmp0 + tmp1 + tmp2 + tmp3
            end do
        end do
        rho_jdos_part(iq) = rho_jdos_part(iq)*delta_kx*delta_ky
        rho_jsdos_part(iq) = rho_jsdos_part(iq)*delta_kx*delta_ky
    end do

    call MPI_REDUCE(rho_jdos_part,rho_jdos,size(rho_jdos_part),MPI_real8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    call MPI_REDUCE(rho_jsdos_part,rho_jsdos,size(rho_jsdos_part),MPI_real8,MPI_SUM,0,MPI_COMM_WORLD,ierr)

    if (cpuid == 0) then
        write(*,'(a)') '   the integral was done.'
        write(*,'(a)') "=> Now it's ready to write qpi_jdos.dat and qpi_jdos.gnu>"
        write(*,'(a)') "   and qpi_jsdos.dat and qpi_jsdos.gnu>."
        open(11,file = 'qpi_jdos.dat')
        open(12,file = 'qpi_jsdos.dat')
        do iq = 1, Nq
            if (iq==1) then
                qx_min_shape = qx_min+(q12(1,iq)-1.0/2)*delta_qx
                qy_min_shape = qy_min+(q12(2,iq)-1.0/2)*delta_qy
            end if
            if (iq==Nq) then
                qx_max_shape = qx_min+(q12(1,iq)-1.0/2)*delta_qx
                qy_max_shape = qy_min+(q12(2,iq)-1.0/2)*delta_qy
            end if
            write(11,*) qx_min+(q12(1,iq)-1.0/2)*delta_qx, qy_min+(q12(2,iq)-1.0/2)*delta_qy, rho_jdos(iq)
            write(12,*) qx_min+(q12(1,iq)-1.0/2)*delta_qx, qy_min+(q12(2,iq)-1.0/2)*delta_qy, rho_jsdos(iq)
            if (mod(iq, n_qy)==0) write(11, *)' '
            if (mod(iq, n_qy)==0) write(12,*) ' '
        end do
        close(11)
        close(12)

        open(13,file = 'qpi_jdos.gnu')
        write(13,'(a)') 'set encoding iso_8859_1'
        write(13,'(a)') 'set terminal png truecolor enhanced font ",50" size 1920, 1680'
        write(13,'(a)') "set output 'qpi_jdos.png'"
        write(13,'(a)') 'unset ztics'
        write(13,'(a)') 'unset key'
        write(13,'(a)') 'set pm3d'
        write(13,'(a)') 'set border lw 6'
        write(13,'(a)') 'set size ratio -1'
        write(13,'(a)') 'set view map'
        write(13,'(a)') 'set xtics'
        write(13,'(a)') 'set ytics'
        write(13,'(a)') 'set xlabel "q_x"'
        write(13,'(a)') 'set ylabel "q_y"'
        write(13,'(a)') 'set ylabel offset 1,0'
        write(13,'(a)') 'set colorbox'
        write(13,'(a)') 'unset cbtics'
        write(13,'(a,f8.6,a,f8.6,a)') 'set xrange [', qx_min_shape, ':', qx_max_shape, ']'
        write(13,'(a,f8.6,a,f8.6,a)') 'set yrange [', qy_min_shape, ':', qy_max_shape, ']'
        write(13,'(a)') 'set pm3d interpolate 2,2'
        write(13,'(a)') "splot 'qpi_jdos.dat' u 1:2:3 w pm3d"
        close(13)

        open(14,file = 'qpi_jsdos.gnu')
        write(14,'(a)') 'set encoding iso_8859_1'
        write(14,'(a)') 'set terminal png truecolor enhanced font ",50" size 1920, 1680'
        write(14,'(a)') "set output 'qpi_jsdos.png'"
        write(14,'(a)') 'unset ztics'
        write(14,'(a)') 'unset key'
        write(14,'(a)') 'set pm3d'
        write(14,'(a)') 'set border lw 6'
        write(14,'(a)') 'set size ratio -1'
        write(14,'(a)') 'set view map'
        write(14,'(a)') 'set xtics'
        write(14,'(a)') 'set ytics'
        write(14,'(a)') 'set xlabel "q_x"'
        write(14,'(a)') 'set ylabel "q_y"'
        write(14,'(a)') 'set ylabel offset 1,0'
        write(14,'(a)') 'set colorbox'
        write(14,'(a)') 'unset cbtics'
        write(14,'(a,f8.6,a,f8.6,a)') 'set xrange [', qx_min_shape, ':', qx_max_shape, ']'
        write(14,'(a,f8.6,a,f8.6,a)') 'set yrange [', qy_min_shape, ':', qy_max_shape, ']'
        write(14,'(a)') 'set pm3d interpolate 2,2'
        write(14,'(a)') "splot 'qpi_jsdos.dat' u 1:2:3 w pm3d"
        close(14)
    end if

    !call mpi_barrier(MPI_COMM_WOLRD, ierr)

    deallocate(G0k_matrix)
    deallocate(G0kq_matrix)
    deallocate(q12)
    deallocate(rho_jdos,rho_jdos_part,rho_jsdos,rho_jsdos_part)

    return
end subroutine

subroutine Spin_Orbit_QPI(cpuid,num_cpu)
    !
    ! For spin-orbit scattering
    ! The scattering potential: V_{k,k'} = V_{0}[I + ic(k \times k')\dot \sigma]
    ! where 'c' is the effective spin-orbit coupling parameter and denotes.
    ! strength of spin-orbit sacattering.
    ! V_{k,k'} was divided into two part: V_{k,k'} = \sum_{j} u_{j}(k) v_{j}(k')
    ! In 2D, V_{k,k'} = V_{0}[ I + ic(k_{x}k'_{y} - k_{y}k'_{x})\sigma_z ].
    ! u(k) = (I k_{x}I k_{y}I), v(k) = U*u(k).
    ! U(1,1) = V_{0}I, U(2,3) = icV_{0}\sigma_z, U(3,2) = -icV_{0}\sigma_z.
    ! ==> More details see Physical Review B 95, 115307(2017) APPENDIX A
    !
    use para
    implicit none
    include 'mpif.h'
    
    integer*4 num_cpu, cpuid, ierr
    complex*16, allocatable, dimension(:,:,:,:) :: G0k_matrix,G0kq_matrix,kxkyG0k_matrix
    complex*16, allocatable, dimension(:,:,:,:) :: kxG0k_matrix,kyG0k_matrix
    complex*16, allocatable, dimension(:,:,:,:) :: kx2G0k_matrix,ky2G0k_matrix
    complex*16, dimension(2,2) :: int_G0k,T_matrix,G_matrix
    complex*16, dimension(2,2) :: int_kxG0k,int_kyG0k,int_kx2G0k,int_ky2G0k,int_kxkyG0k
    complex*16, dimension(6,6) :: U_matrix, Mw, Mw_tmp, Mw_new, Mw_old, I_matrix
    real*8, dimension(0:3) :: dk_tmp,dk_tmp1
    integer*4 :: ikx,iky,iqx,iqy,ikpx,ikpy,ikp,i,j,iq
    real*8 :: kx,ky,kpx,kpy,qx_min_shape,qx_max_shape,qy_min_shape,qy_max_shape
    real*8, dimension(2,6) :: uk
    real*8, dimension(6,2) :: vkq
    real*8, allocatable :: rho(:),rho_part(:)
    integer*4, allocatable :: q12(:,:)
    
    allocate(G0k_matrix(2,2,1:n_kx,1:n_ky),kxkyG0k_matrix(2,2,1:n_kx,1:n_ky))
    allocate(kxG0k_matrix(2,2,1:n_kx,1:n_ky),kyG0k_matrix(2,2,1:n_kx,1:n_ky))
    allocate(kx2G0k_matrix(2,2,1:n_kx,1:n_ky),ky2G0k_matrix(2,2,1:n_kx,1:n_ky))
    allocate(G0kq_matrix(2,2,1:n_kx+n_qx,1:n_ky+n_qy))
    allocate(q12(2,1:Nq))
    allocate(rho(1:Nq),rho_part(1:Nq))
    
    if (cpuid == 0 ) then
        write(*,'(a)') '=> The integral of Spin-Orbit QPI start:'
        do ikx = 1,n_kx
            do iky = 1,n_ky
                kx = (ikx - 1.0/2)*delta_kx + kx_min
                ky = (iky - 1.0/2)*delta_ky + ky_min
                dk_tmp = dk(kx,ky)
                G0k_matrix(:,:,ikx,iky) = matinv(complex(omega,eta)*sigma_0 - &
                & (dk_tmp(0)*sigma_0 + dk_tmp(1)*sigma_x + dk_tmp(2)*sigma_y + dk_tmp(3)*sigma_z))
                kxG0k_matrix(:,:,ikx,iky) = kx*G0k_matrix(:,:,ikx,iky)
                kyG0k_matrix(:,:,ikx,iky) = ky*G0k_matrix(:,:,ikx,iky)
                kx2G0k_matrix(:,:,ikx,iky) = kx*kx*G0k_matrix(:,:,ikx,iky)
                ky2G0k_matrix(:,:,ikx,iky) = ky*ky*G0k_matrix(:,:,ikx,iky)
                kxkyG0k_matrix(:,:,ikx,iky) = kx*ky*G0k_matrix(:,:,ikx,iky)
            end do
        end do
        
        ! \sum_{k} G_{0}(k,w): int_G0k
        int_G0k(1,1) = sum(G0k_matrix(1,1,:,:))
        int_G0k(1,2) = sum(G0k_matrix(1,2,:,:))
        int_G0k(2,1) = sum(G0k_matrix(2,1,:,:))
        int_G0k(2,2) = sum(G0k_matrix(2,2,:,:))
        int_G0k = int_G0k*delta_kx*delta_ky

        ! \sum_{k} k_{x}G_{0}(k,w): int_kxG0k
        int_kxG0k(1,1) = sum(kxG0k_matrix(1,1,:,:))
        int_kxG0k(1,2) = sum(kxG0k_matrix(1,2,:,:))
        int_kxG0k(2,1) = sum(kxG0k_matrix(2,1,:,:))
        int_kxG0k(2,2) = sum(kxG0k_matrix(2,2,:,:))
        int_kxG0k = int_kxG0k*delta_kx*delta_ky
        
        ! \sum_{k} k_{y}G_{0}(k,w): int_kyG0k
        int_kyG0k(1,1) = sum(kyG0k_matrix(1,1,:,:))
        int_kyG0k(1,2) = sum(kyG0k_matrix(1,2,:,:))
        int_kyG0k(2,1) = sum(kyG0k_matrix(2,1,:,:))
        int_kyG0k(2,2) = sum(kyG0k_matrix(2,2,:,:))
        int_kyG0k = int_kyG0k*delta_kx*delta_ky
        
        ! \sum_{k} k_{x}^{2}G_{0}(k,w): int_kx2G0k
        int_kx2G0k(1,1) = sum(kx2G0k_matrix(1,1,:,:))
        int_kx2G0k(1,2) = sum(kx2G0k_matrix(1,2,:,:))
        int_kx2G0k(2,1) = sum(kx2G0k_matrix(2,1,:,:))
        int_kx2G0k(2,2) = sum(kx2G0k_matrix(2,2,:,:))
        int_kx2G0k = int_kx2G0k*delta_kx*delta_ky
        
        ! \sum_{k} k_{y}^{2}G_{0}(k,w): int_ky2G0k
        int_ky2G0k(1,1) = sum(ky2G0k_matrix(1,1,:,:))
        int_ky2G0k(1,2) = sum(ky2G0k_matrix(1,2,:,:))
        int_ky2G0k(2,1) = sum(ky2G0k_matrix(2,1,:,:))
        int_ky2G0k(2,2) = sum(ky2G0k_matrix(2,2,:,:))
        int_ky2G0k = int_ky2G0k*delta_kx*delta_ky
        
        ! \sum_{k} k_{x}k_{y}G_{0}(k,w): int_kxkyG0k
        int_kxkyG0k(1,1) = sum(kxkyG0k_matrix(1,1,:,:))
        int_kxkyG0k(1,2) = sum(kxkyG0k_matrix(1,2,:,:))
        int_kxkyG0k(2,1) = sum(kxkyG0k_matrix(2,1,:,:))
        int_kxkyG0k(2,2) = sum(kxkyG0k_matrix(2,2,:,:))
        int_kxkyG0k = int_kxkyG0k*delta_kx*delta_ky
        
        ! \sum_{k} M(k,w)_matrix,  M(6,6)
        Mw_tmp = 0d0        
        Mw_tmp(1:2,1:2) = int_G0k
        Mw_tmp(1:2,3:4) = int_kxG0k
        Mw_tmp(1:2,5:6) = int_kyG0k
        Mw_tmp(3:4,1:2) = int_kxG0k
        Mw_tmp(3:4,3:4) = int_kx2G0k
        Mw_tmp(3:4,5:6) = int_kxkyG0k
        Mw_tmp(5:6,1:2) = int_kyG0k
        Mw_tmp(5:6,3:4) = int_kxkyG0k
        Mw_tmp(5:6,5:6) = int_ky2G0k
        
        ! U matrix
        U_matrix = 0d0
        U_matrix(1:2,1:2) = V0*sigma_0
        U_matrix(3:4,5:6) = complex(0,c)*V0*sigma_z
        U_matrix(5:6,3:4) = -complex(0,c)*V0*sigma_z
        
        ! M(w)
        Mw = 0d0
        Mw = matmul(U_matrix,Mw_tmp)
        
        Mw_new = 0d0
        Mw_old = 0d0
        I_matrix = 0d0

        DO i = 1,6
            Mw_old(i,i) = 1d0
            I_matrix(i,i) =1d0
        END DO
        
        DO i = 1,50
            Mw_new = matmul(Mw,Mw_old)
            Mw_old = Mw_new
        END DO
        Mw = Mw_new + I_matrix ! Unity Matrix
    end if
   
    call MPI_Bcast(Mw,size(Mw),MPI_double_complex,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(G0k_matrix,size(G0k_matrix),MPI_double_complex,0,MPI_COMM_WORLD,ierr)
    
    iq = 0
    do i = 1,n_qx
        do j = 1,n_qy
            iq = iq + 1
            q12(1,iq) = i
            q12(2,iq) = j
        end do
    end do
    
    rho_part = 0d0
    rho = 0d0
    
    do iq = cpuid + 1,Nq, num_cpu
        iqx = q12(1,iq)
        iqy = q12(2,iq)
        do ikx = 1,n_kx
            do iky = 1,n_ky
                ikpx = ikx + iqx
                ikpy = iky + iqy
                kx = (ikx - 1.0/2)*delta_kx + kx_min
                ky = (iky - 1.0/2)*delta_ky + ky_min
                kpx = (ikx - 1.0/2)*delta_kx + kx_min + (iqx - 1.0/2)*delta_qx + qx_min
                kpy = (iky - 1.0/2)*delta_ky + ky_min + (iqy - 1.0/2)*delta_qy + qy_min
                
                dk_tmp1 = dk(kpx,kpy)
                G0kq_matrix(:,:,ikpx,ikpy) = matinv(complex(omega,eta)*sigma_0-&
                &(dk_tmp1(0)*sigma_0+dk_tmp1(1)*sigma_x+dk_tmp1(2)*sigma_y+dk_tmp1(3)*sigma_z))
                
                uk(1:2,1:2) = sigma_0
                uk(1:2,3:4) = kx*sigma_0
                uk(1:2,5:6) = ky*sigma_0
                
                vkq(1:2,1:2) = sigma_0
                vkq(3:4,1:2) = kpx*sigma_0
                vkq(5:6,1:2) = kpy*sigma_0
                
                T_matrix = matmul(matmul(uk,Mw),vkq)
                
                G_matrix = matmul(matmul(G0k_matrix(:,:,ikx,iky),T_matrix),G0kq_matrix(:,:,ikpx,ikpy))
                rho_part(iq)= rho_part(iq) + imagpart(G_matrix(1,1)+G_matrix(2,2))
            end do
        end do
        rho_part(iq) = rho_part(iq)*delta_kx*delta_ky
    end do

    call MPI_REDUCE(rho_part,rho,size(rho_part),MPI_real8,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    
    if (cpuid == 0) then
        write(*,'(a)') '   the integral was done.'
        write(*,'(a)') "=> Now it's ready to write qpi_SO.dat and qpi_SO.gnu>."
        open(10,file = 'qpi_SO.dat')
        do iq = 1, Nq
            if (iq==1) then
                qx_min_shape = qx_min+(q12(1,iq)-1.0/2)*delta_qx
                qy_min_shape = qy_min+(q12(2,iq)-1.0/2)*delta_qy
            end if
            if (iq==Nq) then
                qx_max_shape = qx_min+(q12(1,iq)-1.0/2)*delta_qx
                qy_max_shape = qy_min+(q12(2,iq)-1.0/2)*delta_qy
            end if
            write(10,*) qx_min+(q12(1,iq)-1.0/2)*delta_qx, qy_min+(q12(2,iq)-1.0/2)*delta_qy, rho(iq)
            if (mod(iq, n_qy)==0) write(10, *)' '
        end do
        close(10)

        open(11,file = 'qpi_SO.gnu')
        write(11,'(a)') 'set encoding iso_8859_1'
        write(11,'(a)') 'set terminal png truecolor enhanced font ",50" size 1920, 1680'
        write(11,'(a)') "set output 'qpi_SO.png'"
        write(11,'(a)') 'unset ztics'
        write(11,'(a)') 'unset key'
        write(11,'(a)') 'set pm3d'
        write(11,'(a)') 'set border lw 6'
        write(11,'(a)') 'set size ratio -1'
        write(11,'(a)') 'set view map'
        write(11,'(a)') 'set xtics'
        write(11,'(a)') 'set ytics'
        write(11,'(a)') 'set xlabel "q_x"'
        write(11,'(a)') 'set ylabel "q_y"'
        write(11,'(a)') 'set ylabel offset 1,0'
        write(11,'(a)') 'set colorbox'
        write(11,'(a)') 'unset cbtics'
        write(11,'(a,f8.6,a,f8.6,a)') 'set xrange [', qx_min_shape, ':', qx_max_shape, ']'
        write(11,'(a,f8.6,a,f8.6,a)') 'set yrange [', qy_min_shape, ':', qy_max_shape, ']'
        write(11,'(a)') 'set pm3d interpolate 2,2'
        write(11,'(a)') "splot 'qpi_SO.dat' u 1:2:3 w pm3d"
        close(11)
    end if

    !call mpi_barrier(MPI_COMM_WOLRD, ierr)

    deallocate(G0k_matrix,kxG0k_matrix,kyG0k_matrix)
    deallocate(kx2G0k_matrix,ky2G0k_matrix,kxkyG0k_matrix)
    deallocate(G0kq_matrix)
    deallocate(q12)
    deallocate(rho,rho_part)

    return
end subroutine

subroutine timestamp()
    implicit none
    integer*4 y,m,d,h,n,s
    character(len=9),parameter, dimension(12) :: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
    integer*4 values(8)
    
    call date_and_time(values = values)
    y = values(1)
    m = values(2)
    d = values(3)
    h = values(5)
    n = values(6)
    s = values(7)
    write( *, '(a, i2.2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3)' )&
    &'=> The current time is: ',d, trim ( month(m) ), y, h, ':', n,':',s, '.'
    return
end subroutine

