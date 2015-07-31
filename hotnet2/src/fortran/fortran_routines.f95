!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Instructions
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! To compile for use with Python/F2PY, run either the included setup.py script
! or the following command:
!
!    f2py -c fortran_routines.f95 -m fortran_routines
!
! Both approaches, which are generally equivalent, require a Fortran compiler
! and NumPy, which contains F2PY.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Routines
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
! compute_sim(M,infmat,h,indices,p,q)
!
!   condenses the similarity matrix given the influence matrix, heat scores,
!   and indices indicating for which genes the influence matrix and heat scores
!   overlap; for HotNet2
!

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

!
! compute_sim_classic(M,infmat,h,indices,p,q)
!
!   condenses the similarity matrix given the influence matrix, heat scores,
!   and indices indicating for which genes the influence matrix and heat scores
!   overlap; for HotNet
!

subroutine compute_sim_classic(M,infmat,h,indices,p,q)

    implicit none

    integer, intent(in) :: p, q
    integer, intent(in) :: indices(q)
    double precision, intent(in) :: infmat(p,p), h(q)
    double precision, intent(out) :: M(q,q)

    integer :: i, j

    do j=1,q
        do i=j,q
            M(i,j) = min(infmat(indices(i),indices(j)), infmat(indices(j),indices(i)))*max(h(i), h(j))
            M(j,i) = M(i,j)
        end do
    end do

end subroutine compute_sim_classic

!
! condense_graph(B,A,V,k,m,n)
!
!   condenses the graph given by the weighted adjacency matrix A into a graph
!   given by the weighted adjacency matrix B by contracting the vertices in V
!

subroutine condense_graph(B,A,V,k,m,n)

    implicit none

    integer, intent(in) :: m, n
    integer, intent(in) :: V(m), k(n+1)
    double precision, intent(in) :: A(m,m)

    double precision, intent(out) :: B(n,n)

    integer :: i, j

    ! m: number of vertices; dimension of A
    ! n: number of SCCs; dimension of B
    ! V: vertices of the components; component(k(i):k(i+1)-1) are the vertices
    !    in component i
    ! k: vertex indices of the components; k(i) is the index of the first
    !    vertex in component i
    ! A: weighted adjacency matrix
    ! B: condensed weighted adjacency matrix

    do i=1,n
        do j=1,n
            if (i/=j) then
                B(i,j) = minimum_slice(V(k(i):k(i+1)-1),V(k(j):k(j+1)-1),k(i+1)-k(i),k(j+1)-k(j))
            else
                B(i,j) = 0.d0
            end if
        end do
    end do

contains

    function minimum_slice(columns,rows,p,q)

        integer, intent(in) :: p, q
        integer, intent(in) :: columns(p), rows(q)

        double precision :: minimum_slice

        double precision :: weight
        integer :: ii, jj

        minimum_slice = huge(0.d0)

        do jj=1,q
            do ii=1,p
                weight = A(columns(ii),rows(jj))
                if (weight>0.d0) then
                    if (weight<minimum_slice) then
                        minimum_slice = weight
                    end if
                end if
            end do
        end do

        if (minimum_slice==huge(0.d0)) then
            minimum_slice = 0.d0
        end if

    end function minimum_slice

end subroutine condense_graph

!
! remove_edges(B,A,weight,m,n)
!
!   removes all edges from weighted adjacency matrix A with weights strictly
!   greater than the given weight
!

subroutine remove_edges(B,A,weight,m,n)

    implicit none

    integer, intent(in) :: m, n
    double precision, intent(in) :: A(m,n), weight

    double precision, intent(out) :: B(m,n)

    integer :: i, j

    do j=1,n
        do i=1,m
            if (A(i,j)<=weight) then
                B(i,j) = A(i,j)
            else
                B(i,j) = 0.d0
            end if
        end do
    end do

end subroutine remove_edges

!
! slice_array(B,A,columns,rows,m,n,p,q)
!
!   finds the slice of the adjacency matrix A given by the columns and rows
!

subroutine slice_array(B,A,columns,rows,m,n,p,q)

    implicit none

    integer, intent(in) :: m, n, p, q
    integer, intent(in) :: columns(p), rows(q)
    double precision, intent(in) :: A(m,n)

    double precision, intent(out) :: B(p,q)

    integer :: i, j

    do j=1,q
        do i=1,p
            B(i,j) = A(columns(i),rows(j))
        end do
    end do

end subroutine slice_array

!
! strongly_connected_components(scc_vertices,A,n)
!
!   finds the strongly connected components of the graph given by the weighted
!   adjacency matrix A; adapted from NetworkX's source code with many variable
!   names unchanged
!

subroutine strongly_connected_components(scc_vertices,A,n)

    implicit none

    integer, intent(in) :: n
    double precision, intent(in) :: A(n,n)

    integer, intent(out) :: scc_vertices(n)

    integer :: i, j, k, p, q, r, v, w, source

    integer :: preorder_collection(n), preorder(n), lowlink(n), scc_queue(n), queue(n)
    logical :: done, scc_found(n)

    i=0
    j=0
    p=0
    q=0
    r=0
    scc_found = .false.

    do source=1,n
        if (scc_found(source) .eqv. .false.) then
            queue(1) = source
            p = p+1
            do while (p>0)

                v = queue(p)
                if (included(v,preorder_collection(1:j),j) .eqv. .false.) then
                    j = j+1
                    preorder_collection(j) = v
                    i = i+1
                    preorder(v) = i
                end if
                done = .true.
                do w=1,n
                    if (A(v,w)/=0.d0) then
                        if (included(w,preorder_collection(1:j),j) .eqv. .false.) then
                            p = p+1
                            queue(p) = w
                            done = .false.
                            exit
                        end if
                    end if
                end do

                if (done .eqv. .true.) then
                    lowlink(v) = preorder(v)
                    do w=1,n
                        if (A(v,w)/=0.d0) then
                            if (scc_found(w) .eqv. .false.) then
                                if (preorder(w) > preorder(v)) then
                                    lowlink(v) = min(lowlink(v),lowlink(w))
                                else
                                    lowlink(v) = min(lowlink(v),preorder(w))
                                end if
                            end if
                        end if
                    end do
                    p = p-1
                    if (lowlink(v)==preorder(v)) then
                        scc_found(v) = .true.
                        r = r+1
                        scc_vertices(v) = r
                        do while (q>0 .and. preorder(scc_queue(q))>preorder(v))
                            k = scc_queue(q)
                            q = q-1
                            scc_found(k) = .true.
                            scc_vertices(k) = r
                        end do
                    else
                        q = q+1
                        scc_queue(q) = v
                    end if
                end if

            end do
        end if
    end do

    contains

        function included(s,t,m)

            integer, intent(in) :: s, m
            integer, intent(in) :: t(m)

            logical :: included
            integer :: ii

            included = .false.

            do ii=1,m
                if (s==t(ii)) then
                    included = .true.
                    exit
                end if
            end do

        end function included

end subroutine strongly_connected_components