subroutine subroutine_sum(fsize,fvec,sum)
 
    integer fsize,i
    real*8 fvec(fsize)
    real*8 sum
 
    sum=0.0
 
    do i=1,fsize
      sum=sum+fvec(i)
    end do
 
end