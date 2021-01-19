!!!!call of these two subroutinew   
    program laglange_interp
    implicit none
    integer k,j
    real cu(0:120)
    open(10,file='data.txt',status='replace')
    do k=1,15
        write(10,"('case(',I2,')')")k
            call lag_inter(k,cu)
            write(10,"('coef(:)=(/')", Advance='NO')
            do j=0,k-2
                write(10,"(D26.18,',')", Advance='NO') cu(j)
            enddo
            write(10,"(D26.18,A2)")cu(k-1),'/)'
        
    enddo
    close(10);
    end 
!!!!!!!!!!!!!!!!!!!!!!