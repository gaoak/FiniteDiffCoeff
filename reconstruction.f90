    subroutine recon_coef_u(r,k,coef)
    !reconstruction coeficient for uniform grids
    implicit none
    real*8 coef(0:k-1)
    integer r,k
    integer l,m,q,j
    real*8 dentemp,numtemp,temp
    do j=0,k-1
        coef(j)=0.d0
    enddo
    do j=0,k-1
        do m=j+1,k
            if(m==r+1) cycle
            dentemp = 1.d0
            numtemp = 1.d0
            do l=0,k
                if(l==m) cycle
                dentemp = dentemp*dble(m-l)
                if(l==r+1) cycle
                numtemp = numtemp*dble(r+1-l)
            enddo
          coef(j)=coef(j)+numtemp/dentemp
        enddo
        if(j<=r) then
            dentemp = 1.d0
            do l=0,k
                if(l==r+1) cycle
                dentemp = dentemp*dble(r+1-l)
            enddo
            numtemp = 0.d0
            do l=0,k
                if(l==r+1)cycle
                temp = 1.d0
                do q=0,k
                    if(q==r+1)cycle
                    if(q==l)cycle
                    temp = temp*dble(r+1-q)
                enddo
                numtemp = numtemp+temp
            enddo
            coef(j)=coef(j)+numtemp/dentemp
        endif
    enddo
    end subroutine recon_coef_u
    
    subroutine recon_coef_a(r,k,coef,stencil)
    !reconstruction coefficient for arbitrary grids
    implicit none
    integer r,k
    real*8 coef(0:k-1), stencil(-r:k-r) !stencil: -r-0.5, ..., s+1-0.5
    integer l,m,q,j
    real*8 dentemp,numtemp,temp
    do j=0,k-1
        coef(j)=0.d0
    enddo
    do j=0,k-1
        do m=j+1,k
            if(m==r+1) cycle
            dentemp = 1.d0
            numtemp = 1.d0
            do l=0,k
                if(l==m) cycle
                dentemp = dentemp*(stencil(m-r)-stencil(l-r))
                if(l==r+1) cycle
                numtemp = numtemp*(stencil(1)-stencil(l-r))
            enddo
          coef(j)=coef(j)+numtemp/dentemp
        enddo
        if(j<=r) then
            dentemp = 1.d0
            do l=0,k
                if(l==r+1) cycle
                dentemp = dentemp*(stencil(1)-stencil(l-r))
            enddo
            numtemp = 0.d0
            do l=0,k
                if(l==r+1)cycle
                temp = 1.d0
                do q=0,k
                    if(q==r+1)cycle
                    if(q==l)cycle
                    temp = temp*(stencil(1)-stencil(q-r))
                enddo
                numtemp = numtemp+temp
            enddo
            coef(j) = coef(j)+numtemp/dentemp
        endif
        coef(j) = coef(j)*(stencil(j-r+1)-stencil(j-r))
    enddo
    end subroutine recon_coef_a

!!!!call of these two subroutinew   
!    program reconstruction_coefficients
!    implicit none
!    integer k,r,j
!    real*8 cu(0:12),ca(0:12),stencil(12)
!    do j=1,12
!       stencil(j) = dble(j)*0.00000001d0
!    enddo
!    open(10,file='data.txt',status='replace')
!    write(10,"(' k,  r,')", Advance='NO')
!    do j=0,5
!        write(10,"('       j=',I1,',      ')", Advance='NO')j
!    enddo
!    write(10,*)
!    do k=1,6
!        do r=-1,k-1
!            call recon_coef_u(r,k,cu)
!            call recon_coef_a(r,k,ca,stencil)
!            write(10,"(' ',I1,', ',I2,',')", Advance='NO')k,r
!            do j=0,k-1
!                write(10,"(E16.7,',')", Advance='NO') (cu(j)-ca(j))/cu(j)
!            enddo
!            write(10,*)
!        enddo
!    enddo
!    close(10);
!    end program reconstruction_coefficients
!!!!!!!!!!!!!!!!!!!!!

    subroutine div_diff(k,u,stencil,divdiff)
    ! evaluate divid differences from 0 degree to k-1
    implicit none
    integer k
    real*8 u(1-k:k-1), stencil(1-k:k-1), divdiff(1-k:k-1,0:k-1)
    integer j,m
    do j = 1-k,k-1
        divdiff(j,0) = u(j)
    enddo
    do m = 1,k-1
        do j = 1-k,k-1-m
            divdiff(j,m) = (divdiff(j+1,m-1)-divdiff(j,m-1))/(stencil(j+m)-stencil(j))
        enddo
    enddo
    end subroutine div_diff
    subroutine div_diff_V(k,dV,l,stencil,s,divdiff)
    ! evaluate V's divid differences from l degree to k-1, kth order accuracy
    ! evaluate from points s
    implicit none
    integer k,l,s
    real*8 dV(s:k-1-l), stencil(s:k-1), divdiff(s:k-1-l,l:k-1)
    integer j,m
    do j = s,k-1-l
        divdiff(j,l) = dV(j)
    enddo
    do m = l+1,k-1
        do j = s,k-1-m
            divdiff(j,m) = (divdiff(j+1,m-1)-divdiff(j,m-1))/(stencil(j+m)-stencil(j))
        enddo
    enddo
    end subroutine div_diff_V
    subroutine div_diff_V_u(k,dV,l,s,divdiff)
    ! evaluate V's divid differences from l degree to k-1, kth order accuracy
    ! evaluate from points s
    implicit none
    integer k,l,s
    real*8 dV(s:k-1-l), divdiff(s:k-1-l,l:k-1)
    integer j,m
    do j = s,k-1-l
        divdiff(j,l) = dV(j)
    enddo
    do m = l+1,k-1
        do j = s,k-1-m
            divdiff(j,m) = (divdiff(j+1,m-1)-divdiff(j,m-1))/dble(m)
        enddo
    enddo
    end subroutine div_diff_V_u
!!!!!vilidation of this subroutine
!    program Newton_interpolation
!    implicit none
!    integer k,j
!    real*8 u(-20:20),au(-20:20),stencil(-20:20),detx,x,temp
!    real*8, external::func,inter_Newton
!    detx = 0.1d0
!    do j=-20,20
!        stencil(j) = dble(j)*detx
!        u(j) = func(stencil(j))
!    enddo
!    k = 4
!    do while(1)
!        write(*,"('x = ')",Advance = 'NO')
!        read(*,*) x
!        if(x>100.d0) exit
!        write(*,"('Error = ',F10.5,'(O^',I2,')')") (inter_Newton(k,u(1-k:k-1),stencil(1-k:k-1),x)-func(x))/(detx**k),k
!    enddo
!    end program Newton_interpolation
!    
!    real*8 function func(x)
!    real*8 x
!    func = 1.d0+2.d0*x+3.d0*x**2+4.d0*x**3
!    return
!    end function func
!    real*8 function inter_Newton(k,u,stencil,x)
!    ! k is the order of accuracy
!    ! only k points are used in interpolating u(0:k-1)
!    integer k 
!    real*8 u(1-k:k-1),stencil(1-k:k-1),x
!    integer j, m
!    real*8 divdiff(1-k:k-1,0:k-1), temp
!    call div_diff(k,u,stencil,divdiff)
!    inter_Newton = 0.d0
!    do j=0,k-1
!        temp = 1.d0
!        do m=0,j-1
!            temp = temp*(x-stencil(m))
!        enddo
!        inter_Newton = inter_Newton + temp*divdiff(0,j)
!    enddo
!    return
!    end function inter_Newton
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine ENO_recon_a(k,u_bar,stencil,N, u_recon_m,u_recon_p)! not finished yet
    ! node points stars from 0 to N
    ! cell number is N starting from 0 to N-1
    implicit none
    integer k,j, r, m, N
    real*8 u_bar(-k:k-1+N),stencil(-k:k+N)
    real*8 u_recon_p(0:N),u_recon_m(0:N)
    real*8 coef(0:k-1,-1:k-1)
    do j = -1,k-1
        call recon_coef_u(j,k,coef(0:k-1,j))
    enddo
    ! reconstruction
    r = k/2
    do j = 0,N-1
        call ENO_select(k,u_bar(j+1-k:k-1+j),stencil(j+1-k:k+j),r)
        u_recon_m(j+1) = 0.d0
        u_recon_p(j) = 0.d0
        do m = 0,k-1
            u_recon_m(j+1) = u_recon_m(j+1) + u_bar(j-r+m)*coef(m,r)
            u_recon_p(j) = u_recon_p(j) + u_bar(j-r+m)*coef(m,r-1)
        enddo
    enddo
    j = -1
111    call ENO_select(k,u_bar(j+1-k:k-1+j),stencil(j+1-k:k+j),r)
    if(j==-1) u_recon_m(j+1) = 0.d0
    if(j==N) u_recon_p(j) = 0.d0
    do m = 0,k-1
        if(j==-1) u_recon_m(j+1) = u_recon_m(j+1) + u_bar(j-r+m)*coef(m,r)
        if(j==N) u_recon_p(j) = u_recon_p(j) + u_bar(j-r+m)*coef(m,r-1)
    enddo
    if(j==-1) then
        j = N
        goto 111
    endif
    end subroutine ENO_recon_a
    
    subroutine ENO_recon_u(k,u_bar,N, u_recon_m,u_recon_p)
    ! node points stars from 0 to N
    ! cell number is N starting from 0 to N-1
    implicit none
    integer k,j, r, m, N
    real*8 u_bar(-k:k-1+N)
    real*8 u_recon_p(0:N),u_recon_m(0:N)
    real*8 coef(0:k-1,-1:k-1)
    do j = -1,k-1
        call recon_coef_u(j,k,coef(0:k-1,j))
    enddo
    ! reconstruction
    r = k/2
    do j = 0,N-1
        call ENO_select_u(k,u_bar(j+1-k:k-1+j),r)
        u_recon_m(j+1) = 0.d0
        u_recon_p(j) = 0.d0
        do m = 0,k-1
            u_recon_m(j+1) = u_recon_m(j+1) + u_bar(j-r+m)*coef(m,r)
            u_recon_p(j) = u_recon_p(j) + u_bar(j-r+m)*coef(m,r-1)
        enddo
    enddo
    j = -1
111    call ENO_select_u(k,u_bar(j+1-k:k-1+j),r)
    if(j==-1) u_recon_m(j+1) = 0.d0
    if(j==N) u_recon_p(j) = 0.d0
    do m = 0,k-1
        if(j==-1) u_recon_m(j+1) = u_recon_m(j+1) + u_bar(j-r+m)*coef(m,r)
        if(j==N) u_recon_p(j) = u_recon_p(j) + u_bar(j-r+m)*coef(m,r-1)
    enddo
    if(j==-1) then
        j = N
        goto 111
    endif
    end subroutine ENO_recon_u
    
    subroutine ENO_select(k,u_bar,stencil,r)
    ! select the stencil and do nothing
    ! k is the order of accuracy, using k cells or k+1 points
    implicit none
    integer k,r
    real*8 u_bar(1-k:k-1),stencil(1-k:k)
    integer j
    real*8 divdiff(1-k:k-1,1:k)
    call div_diff_V(k+1,u_bar,1,stencil,1-k,divdiff)
    r = 0
    do j = 2,k
        if(abs(divdiff(r-1,j))<abs(divdiff(r,j))) r = r-1
    enddo
    r = -r
    end subroutine ENO_select
    subroutine ENO_select_u(k,u_bar,r)
    ! select the stencil and do nothing
    ! k is the order of accuracy, using k cells or k+1 points
    implicit none
    integer k,r
    real*8 u_bar(1-k:k-1)
    integer j
    real*8 divdiff(1-k:k-1,1:k)
    call div_diff_V_u(k+1,u_bar,1,1-k,divdiff)
    r = 0
    do j = 2,k
        if(abs(divdiff(r-1,j))<abs(divdiff(r,j))) r = r-1
    enddo
    r = -r
    end subroutine ENO_select_u
!!!!validation of ENO
!    program ENO_reconstruction
!    implicit none
!    integer, parameter :: k = 2
!    integer j, r, m, N
!    real*8 u_bar(-512+1-k:k-1+512),stencil(-512+1-k:k+512),detx,len
!    real*8 u_star_p(-512+1-k:k+512),u_star_m(-512+1-k:k+512),u_star(-512+1-k:k+512)
!    real*8 coef(0:k-1,-1:k-1)
!    real*8 norm_L1p, norm_L1m, norm_Li
!    real*8 temp1p, temp1m, tempi
!    real*8, external :: func, funcf
!    len = 20.d0
!    N = 2
!1   detx = len/(2.d0*dble(N)+1.d0)
!    do j=-N-k,k+N+1
!        stencil(j) = (dble(j)-0.5d0)*detx
!    enddo ! mean values
!    do j=-N-k,k+N
!        u_bar(j) = (funcf(stencil(j+1)) - funcf(stencil(j)))/(stencil(j+1)-stencil(j))
!        u_star(j) = func(stencil(j))
!    enddo ! mean values
!    
!    call ENO_recon_a(k,u_bar(-N-k:k+N),stencil(-N-k:k+N+1),2*N+1, u_star_m(-N:N+1),u_star_p(-N:N+1))
!    
!    open(10,file='data.dat',status = 'replace')
!    write(10,*)'title = ENO reconstruction'
!    write(10,*)'variables = x, um, up, u'
!    temp1p = norm_L1p
!    temp1m = norm_L1m
!    tempi = norm_Li
!    norm_L1p = 0.d0
!    norm_L1m = 0.d0
!    norm_Li = 0.d0
!    do j = -N,N+1
!        write(10,"(4E18.6)")stencil(j),u_star_m(j),u_star_p(j),u_star(j)
!        if(dabs(u_star_m(j)-u_star(j))>norm_Li) norm_Li=dabs(u_star_m(j)-u_star(j))
!        if(dabs(u_star_p(j)-u_star(j))>norm_Li) norm_Li=dabs(u_star_p(j)-u_star(j))
!        norm_L1p = norm_L1p + dabs(u_star_p(j)-u_star(j))
!        norm_L1m = norm_L1m + dabs(u_star_m(j)-u_star(j))
!    enddo
!    norm_L1p = norm_L1p/dble(2*N)
!    norm_L1m = norm_L1m/dble(2*N)
!    close(10)
!    if(N==2) write(*,"(3A15)") 'detx,','L infinity norm,','order'
!    write(*,"(2(E14.6,','))",advance = 'no') detx,norm_Li
!    if(N==2) then
!        write(*,*)
!    else
!        write(*,"(F14.6,',')") (log(tempi)-log(norm_Li))/(log(dble(2*N+1))-log(dble(N+1)))
!    endif
!    N = N*2
!    if(N<64) goto 1
!    end program ENO_reconstruction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine WENO_recon_u(k,u_bar,N,u_recon_m,u_recon_p)
    ! k is number of cells in ENO, accuracy is 2*k-1
    ! node points stars from 0 to N
    ! cell number is N starting from 0 to N-1, adding k points at each boundarys
    implicit none
    integer k,j, r, m, N
    real*8 u_bar(-k:k-1+N)
    real*8 u_recon_p(0:N),u_recon_m(0:N), tempp, tempm
    real*8 coef(0:k-1,-1:k-1),omegam(0:k-1),omegap(0:k-1)
    do j = -1,k-1
        call recon_coef_u(j,k,coef(0:k-1,j))
    enddo
    ! reconstruction
    do j = 0,N-1
        call WENO_param_u(k,u_bar(j+1-k:k-1+j),omegam,omegap)
        u_recon_m(j+1) = 0.d0
        u_recon_p(j) = 0.d0
        do r = 0,k-1
            tempm = 0.d0
            tempp = 0.d0
            do m = 0,k-1
                tempm = tempm + u_bar(j-r+m)*coef(m,r)
                tempp = tempp + u_bar(j-r+m)*coef(m,r-1)
            enddo
            u_recon_m(j+1) = u_recon_m(j+1) + omegam(r)*tempm
            u_recon_p(j) = u_recon_p(j) + omegap(r)*tempp
        enddo
    enddo
    j = -1
111    call WENO_param_u(k,u_bar(j+1-k:k-1+j),omegam,omegap)
    if(j==-1) u_recon_m(j+1) = 0.d0
    if(j==N)  u_recon_p(j) = 0.d0
    do r = 0,k-1
        tempm = 0.d0
        tempp = 0.d0
        do m = 0,k-1
            if(j==-1) tempm = tempm + u_bar(j-r+m)*coef(m,r)
            if(j==N)  tempp = tempp + u_bar(j-r+m)*coef(m,r-1)
        enddo
        if(j==-1) u_recon_m(j+1) = u_recon_m(j+1) + omegam(r)*tempm
        if(j==N)  u_recon_p(j) = u_recon_p(j) + omegap(r)*tempp
    enddo
    if(j==-1) then
        j = N
        goto 111
    endif
    end subroutine WENO_recon_u
    
    subroutine WENO_param_u(k,u_bar,omegam,omegap)
    implicit none
    integer k,j
    real*8 u_bar(1-k:k-1),omegam(0:k-1),omegap(0:k-1)
    real*8 d(0:k-1),beta(0:k-1),epsilon,tempp,tempm
    epsilon = 1.d-1
    select case(k)
    case(2)
        d(0) = 2.d0/3.d0
        d(1) = 1.d0/3.d0
        beta(0) = (u_bar(1)-u_bar(0))**2
        beta(1) = (u_bar(0)-u_bar(-1))**2
    case(3)
        d(0) = 3.d-1
        d(1) = 6.d-1
        d(2) = 1.d-1
        beta(0) = 13.d0/12.d0*(u_bar(0)-2.d0*u_bar(1)+u_bar(2))**2+&
                   0.25d0*(3.d0*u_bar(0)-4.d0*u_bar(1)+u_bar(2))**2
        beta(1) = 13.d0/12.d0*(u_bar(-1)-2.d0*u_bar(0)+u_bar(1))**2+&
                   0.25d0*(u_bar(-1)-u_bar(1))**2
        beta(2) = 13.d0/12.d0*(u_bar(-2)-2.d0*u_bar(-1)+u_bar(0))**2+&
                   0.25d0*(3.d0*u_bar(0)-4.d0*u_bar(-1)+u_bar(-2))**2
    case(4)
        d(0) = 4.d0/35.d0
        d(1) = 18.d0/35.d0
        d(2) = 12.d0/35.d0
        d(3) = 1.d0/35.d0
        beta(0) = u_bar(0)*(2107.d0*u_bar(0)-9402.d0*u_bar(1)+7042.d0*u_bar(2)-1854.d0*u_bar(3))&
            +u_bar(1)*(11003.d0*u_bar(1)-17246.d0*u_bar(2)+4642.d0*u_bar(3))&
            +u_bar(2)*(7043.d0*u_bar(2)-3882.d0*u_bar(3))+547.d0*u_bar(3)*u_bar(3)
        beta(1) = u_bar(-1)*(547.d0*u_bar(-1)-2522.d0*u_bar(0)+1922.d0*u_bar(1)-494.d0*u_bar(2))&
            +u_bar(0)*(3443.d0*u_bar(0)-5966.d0*u_bar(1)+1602.d0*u_bar(2))&
            +u_bar(1)*(2843.d0*u_bar(1)-1642.d0*u_bar(2))+267.d0*u_bar(2)*u_bar(2)
        beta(2) = u_bar(-2)*(267.d0*u_bar(-2)-1642.d0*u_bar(-1)+1602.d0*u_bar(0)-494.d0*u_bar(1))&
            +u_bar(-1)*(2843.d0*u_bar(-1)-5966.d0*u_bar(0)+1922.d0*u_bar(1))&
            +u_bar(0)*(3443.d0*u_bar(0)-2522.d0*u_bar(1))+547.d0*u_bar(1)*u_bar(1)
        beta(3) = u_bar(-3)*(547.d0*u_bar(-3)-3882.d0*u_bar(-2)+4642.d0*u_bar(-1)-1854.d0*u_bar(0))&
            +u_bar(-2)*(7043.d0*u_bar(-2)-17246.d0*u_bar(-1)+7042.d0*u_bar(0))&
            +u_bar(-1)*(11003.d0*u_bar(-1)-9402.d0*u_bar(0))+2107.d0*u_bar(0)*u_bar(0)
    end select
    tempm = 0.d0
    tempp = 0.d0
    !write(*,*)beta
    do j = 0,k-1
        omegam(j) = d(j)/(epsilon+beta(j))**2
        omegap(j) = d(k-1-j)/(epsilon+beta(j))**2
        tempm = tempm+omegam(j)
        tempp = tempp+omegap(j)
    enddo
    do j = 0,k-1
        omegam(j) = omegam(j)/tempm
        omegap(j) = omegap(j)/tempp
    enddo
    end subroutine WENO_param_u
!!!!!validation of WENO reconstruction
!    program WENO_reconstruction
!    implicit none
!    integer, parameter :: k = 3
!    integer j, r, m, N
!    real*8 u_bar(-2048+1-k:k-1+2048),stencil(-2048+1-k:k+2048),detx,len
!    real*8 u_star_p(-2048+1-k:k+2048),u_star_m(-2048+1-k:k+2048),u_star(-2048+1-k:k+2048)
!    real*8 coef(0:k-1,-1:k-1)
!    real*8 norm_L1p, norm_L1m, norm_Li
!    real*8 temp1p, temp1m, tempi
!    real*8, external :: func, funcf
!    len = 20.d0
!    N = 2
!1   detx = len/(2.d0*dble(N)+1.d0)
!    do j=-N-k,k+N+1
!        stencil(j) = (dble(j)-0.5d0)*detx
!    enddo ! mean values
!    do j=-N-k,k+N
!        u_bar(j) = (funcf(stencil(j+1)) - funcf(stencil(j)))/(stencil(j+1)-stencil(j))
!        u_star(j) = func(stencil(j))
!    enddo ! mean values
!    
!    call WENO_recon_u(k,u_bar(-N-k:k+N),2*N+1, u_star_m(-N:N+1),u_star_p(-N:N+1))
!    
!    open(10,file='data.dat',status = 'replace')
!    write(10,*)'title = ENO reconstruction'
!    write(10,*)'variables = x, um, up, u'
!    temp1p = norm_L1p
!    temp1m = norm_L1m
!    tempi = norm_Li
!    norm_L1p = 0.d0
!    norm_L1m = 0.d0
!    norm_Li = 0.d0
!    do j = -N,N+1
!        write(10,"(4E18.6)")stencil(j),u_star_m(j),u_star_p(j),u_star(j)
!        if(dabs(u_star_m(j)-u_star(j))>norm_Li) norm_Li=dabs(u_star_m(j)-u_star(j))
!        if(dabs(u_star_p(j)-u_star(j))>norm_Li) norm_Li=dabs(u_star_p(j)-u_star(j))
!        norm_L1p = norm_L1p + dabs(u_star_p(j)-u_star(j))
!        norm_L1m = norm_L1m + dabs(u_star_m(j)-u_star(j))
!    enddo
!    norm_L1p = norm_L1p/dble(2*N)
!    norm_L1m = norm_L1m/dble(2*N)
!    close(10)
!    if(N==2) write(*,"(3A15)") 'detx,','L infinity norm,','order'
!    write(*,"(2(E14.6,','))",advance = 'no') detx,norm_Li
!    if(N==2) then
!        write(*,*)
!    else
!        write(*,"(F14.6,',')") (log(tempi)-log(norm_Li))/(log(dble(2*N+1))-log(dble(N+1)))
!    endif
!    N = N*2
!    if(N<16) goto 1
!    end program WENO_reconstruction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine lag_inter(k,cu)
    integer k,i,j
    real*8 cu(0:k-1),tempd,tempn
    do i=0,k-1
        tempd = 1.d0
        tempn = 1.d0
        do j=0,k-1
            if(j.eq.i)cycle
            tempd = tempd*dble(j-i)
            tempn = tempn*dble(j-k)
        enddo
        cu(i) = tempn/tempd
    enddo
    end

!    !!!!call of these two subroutinew   
!    subroutine reconstruction_coefficients()
!        implicit none
!        integer k,r,j
!        real*8 cu(0:12),ca(0:12),stencil(12)
!        do j=1,12
!           stencil(j) = dble(j)*0.00000001d0
!        enddo
!        open(10,file='data.txt',status='replace')
!        write(10,"(' k,  r,')", Advance='NO')
!        do j=0,5
!            write(10,"('       j=',I1,',      ')", Advance='NO')j
!        enddo
!        write(10,*)
!        do k=1,11
!            write(10,"('case(',I1,')')")k
!            do r=-1,k-1
!                call recon_coef_u(r,k,cu)
!                call recon_coef_a(r,k,ca,stencil)
!                write(10,"('coef(:,',I2,')=(/')", Advance='NO')r
!                do j=0,k-2
!                    write(10,"(D26.18,',')", Advance='NO') cu(j)
!                enddo
!                write(10,"(D26.18,A2)")cu(k-1),'/)'
!            enddo
!        enddo
!        close(10);
!        end