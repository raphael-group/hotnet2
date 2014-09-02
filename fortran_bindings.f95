subroutine compute_sim(M,infmat,h,indices,p,q)

    implicit none

    integer, intent(in) :: p, q
    integer, intent(in) :: indices(q)
    double precision, intent(in) :: infmat(p,p), h(q)
    double precision, intent(out) :: M(q,q)

    integer :: i, j

    do j=1,q
        do i=1,q
            M(i,j) = infmat(indices(i),indices(j))*h(j)
        end do
    end do

end subroutine compute_sim

subroutine compute_sim_classic(M,infmat,h,indices,p,q)

    implicit none

    integer, intent(in) :: p, q
    integer, intent(in) :: indices(q)
    double precision, intent(in) :: infmat(p,p), h(q)
    double precision, intent(out) :: M(q,q)

    integer :: i, j

    do j=1,q
        do i=1,q
            M(i,j) = min(infmat(indices(i),indices(j)), infmat(indices(j), indices(i)))*max(h(i), h(j))
        end do
    end do

end subroutine compute_sim_classic
