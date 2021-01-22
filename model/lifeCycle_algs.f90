module lifecycle_algs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MODULE ALGS contains algorithms for solving the dynamic programming problem,
! as well as other algorithms coded by DB/TC to process the simulation.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    USE share
    USE nrtype
    USE nrutil, ONLY : assert_eq, imaxloc, iminloc, nrerror, swap
    USE OMP_LIB
    IMPLICIT NONE
    ! TODO: One day maybe switch from NR and custom functions to a better
    ! controlled library (like Math77). But this is far down the road.

    CONTAINS

    SUBROUTINE amoeba(stateholder, p, y, ftol, pol_opts, vfunc) ! %<
    ! Nelder-Mead NR algorithm (for choice of a', D')
        IMPLICIT NONE

        INTEGER(I4B) :: iter
        REAL(DP), INTENT(IN) :: ftol
        REAL(DP), DIMENSION(3), INTENT(INOUT) :: y
        REAL(DP), DIMENSION(3, 2), INTENT(INOUT) :: p
        LOGICAL, DIMENSION(2), INTENT(IN) :: pol_opts
        REAL(DP), DIMENSION(3) :: yholder
        REAL(DP), DIMENSION(3, 2) :: pholder
        REAL(DP), DIMENSION(8), INTENT(IN) :: stateholder
        INTEGER(I4B), PARAMETER :: ITMAX=5000000
        REAL(DP), PARAMETER :: TINY=1.0e-10
        INTEGER(I4B) :: ihi, ndim
        REAL(DP), DIMENSION(size(p, 2)) :: psum
        REAL(DP), EXTERNAL :: vfunc
        call amoeba_private
        CONTAINS
    !BL
        SUBROUTINE amoeba_private
        IMPLICIT NONE
        INTEGER(I4B) :: i, ilo, inhi
        REAL(DP) :: rtol, ysave, ytry, ytmp
        ndim=assert_eq(size(p, 2), size(p, 1)-1, size(y)-1,'amoeba')
        iter=0
        psum(:)=sum(p(:,:), dim=1)

        pholder=p
        yholder=y

        !write(*,*) "11", stateholder
        do
            ilo=iminloc(y(:))
            ihi=imaxloc(y(:))
            ytmp=y(ihi)
            y(ihi)=y(ilo)
            inhi=imaxloc(y(:))
            y(ihi)=ytmp

            !write(*,*) y(ihi), y(ilo)

            rtol=2.0_dp*abs(y(ihi)-y(ilo))/(abs(y(ihi))+abs(y(ilo))+TINY)
            if (rtol < ftol) then
                call swap(y(1), y(ilo))
                call swap(p(1,:), p(ilo,:))
                RETURN
            end if
            if (iter >= ITMAX) then
            write(*,*) "stateholder", stateholder
            write(*,*) "p", p
            write(*,*) "pholder", pholder
            write(*,*) "y", y
            write(*,*) "yholder", yholder
            call nrerror('ITMAX exceeded in amoeba')
            end if
            !write(*,*) "2", stateholder

            ytry=amotry(-1.0_dp, stateholder)
            iter=iter+1
            if (ytry <= y(ilo)) then

            !write(*,*) stateholder

                ytry=amotry(2.0_dp, stateholder)
                iter=iter+1
            else if (ytry >= y(inhi)) then
                ysave=y(ihi)

                !write(*,*) "3", stateholder

                ytry=amotry(0.5_dp, stateholder)
                iter=iter+1
                if (ytry >= ysave) then
                    p(:,:)=0.5_dp*(p(:,:)+spread(p(ilo,:), 1, size(p, 1)))
                    do i=1, ndim+1

                    !write(*,*) "4", stateholder

                        if (i /= ilo) y(i)=vfunc(p(i,:), stateholder, pol_opts)
                        !write(*,*) y(i)
                    end do
                    iter=iter+ndim
                    psum(:)=sum(p(:,:), dim=1)
                end if
            end if
        end do
        END SUBROUTINE amoeba_private
    !BL
        FUNCTION amotry(fac, stateholder)
        IMPLICIT NONE
        REAL(DP), INTENT(IN) :: fac
        REAL(DP) :: amotry
        REAL(DP) :: fac1, fac2, ytry
        REAL(DP), DIMENSION(8) :: stateholder
        REAL(DP), DIMENSION(size(p, 2)) :: ptry
        fac1=(1.0_dp-fac)/ndim
        fac2=fac1-fac
        ptry(:)=psum(:)*fac1-p(ihi,:)*fac2
        !write(*,*) "private", stateholder
        ytry=vfunc(ptry, stateholder, pol_opts)
        !write(*,*) "done?"
        if (ytry < y(ihi)) then
            y(ihi)=ytry
            psum(:)=psum(:)-p(ihi,:)+ptry(:)
            p(ihi,:)=ptry(:)
        end if
        amotry=ytry
        END FUNCTION amotry
    END SUBROUTINE amoeba ! %>

    FUNCTION brentnew(ax, bx, cx, nxt_pol, value, tol, stateholder, xmin) ! %<
    ! Brent minimization NR algorithm (for choice of a')
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: ax, bx, cx, tol
    LOGICAL, INTENT(IN) :: nxt_pol
    REAL(DP), DIMENSION(8), INTENT(IN) :: stateholder
    REAL(DP), INTENT(OUT) :: xmin
    REAL(DP) :: brentnew
    EXTERNAL :: value
    REAL(DP) :: value



    INTEGER(I4B), PARAMETER :: ITMAX=1000
    REAL(DP), PARAMETER :: CGOLD=0.3819660_dp, ZEPS=1.0e-3_dp*epsilon(ax)
    INTEGER(I4B) :: iter
    REAL(DP) :: a, b, d, e, etemp, fu, fv, fw, fx, p, q, r2, tol1, tol2, u, v, w, x, xm
    a=min(ax, cx)
    b=max(ax, cx)
    v=bx
    w=v
    x=v
    e=0.0
    fx=value(x, stateholder, nxt_pol)
    fv=fx
    fw=fx
    do iter=1, ITMAX
        xm=0.5_dp*(a+b)
        tol1=tol*abs(x)+ZEPS
        tol2=2.0_dp*tol1
        if (abs(x-xm) <= (tol2-0.5_dp*(b-a))) then
            xmin=x
            brentnew=fx
            RETURN
        end if
        if (abs(e) > tol1) then
            r2=(x-w)*(fx-fv)
            q=(x-v)*(fx-fw)
            p=(x-v)*q-(x-w)*r2
            q=2.0_dp*(q-r2)
            if (q > 0.0) p=-p
            q=abs(q)
            etemp=e
            e=d
            if (abs(p) >= abs(0.5_dp*q*etemp) .or. &
                p <= q*(a-x) .or. p >= q*(b-x)) then
                e=merge(a-x, b-x, x >= xm )
                d=CGOLD*e
            else
                d=p/q
                u=x+d
                if (u-a < tol2 .or. b-u < tol2) d=sign(tol1, xm-x)
            end if
        else
            e=merge(a-x, b-x, x >= xm )
            d=CGOLD*e
        end if
        u=merge(x+d, x+sign(tol1, d), abs(d) >= tol1 )
        fu=value(u, stateholder, nxt_pol)
        if (fu <= fx) then
            if (u >= x) then
                a=x
            else
                b=x
            end if
            call shft(v, w, x, u)
            call shft(fv, fw, fx, fu)
        else
            if (u < x) then
                a=u
            else
                b=u
            end if
            if (fu <= fw .or. w == x) then
                v=w
                fv=fw
                w=u
                fw=fu
            else if (fu <= fv .or. v == x .or. v == w) then
                v=u
                fv=fu
            end if
        end if
    end do
    call nrerror('brentnew: exceed maximum iterations')
    CONTAINS
    !BL
    SUBROUTINE shft(a, b, c, d)
    REAL(DP), INTENT(OUT) :: a
    REAL(DP), INTENT(INOUT) :: b, c
    REAL(DP), INTENT(IN) :: d
    a=b
    b=c
    c=d
    END SUBROUTINE shft
    END FUNCTION brentnew ! %>

    FUNCTION brentmindist(ax, bx, cx, value, numhh, tol, xmin) ! %<
    ! Brent minimization again (this is customized to solve the market clearing
    ! price vector)
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: ax, bx, cx, tol
    INTEGER, INTENT(IN) :: numhh
    REAL(DP), INTENT(OUT) :: xmin
    REAL(DP) :: brentmindist
    EXTERNAL :: value
    REAL(DP) :: value



    INTEGER(I4B), PARAMETER :: ITMAX=25
    REAL(DP), PARAMETER :: CGOLD=0.3819660_dp, ZEPS=1.0e-3_dp*epsilon(ax)
    INTEGER(I4B) :: iter
    REAL(DP) :: a, b, d, e, etemp, fu, fv, fw, fx, p, q, r2, tol1, tol2, u, v, w, x, xm
    a=min(ax, cx)
    b=max(ax, cx)
    v=bx
    w=v
    x=v
    e=0.0
    fx=value(x, numhh, 1)
    fv=fx
    fw=fx
    do iter=1, ITMAX
        xm=0.5_dp*(a+b)
        tol1=tol*abs(x)+ZEPS
        tol2=2.0_dp*tol1
        if (abs(x-xm) <= (tol2-0.5_dp*(b-a))) then
            xmin=x
            brentmindist=fx
            RETURN
        end if
        if (abs(e) > tol1) then
            r2=(x-w)*(fx-fv)
            q=(x-v)*(fx-fw)
            p=(x-v)*q-(x-w)*r2
            q=2.0_dp*(q-r2)
            if (q > 0.0) p=-p
            q=abs(q)
            etemp=e
            e=d
            if (abs(p) >= abs(0.5_dp*q*etemp) .or. &
                p <= q*(a-x) .or. p >= q*(b-x)) then
                e=merge(a-x, b-x, x >= xm )
                d=CGOLD*e
            else
                d=p/q
                u=x+d
                if (u-a < tol2 .or. b-u < tol2) d=sign(tol1, xm-x)
            end if
        else
            e=merge(a-x, b-x, x >= xm )
            d=CGOLD*e
        end if
        u=merge(x+d, x+sign(tol1, d), abs(d) >= tol1 )
        fu=value(u, numhh, iter+1)
        if (fu <= fx) then
            if (u >= x) then
                a=x
            else
                b=x
            end if
            call shft(v, w, x, u)
            call shft(fv, fw, fx, fu)
        else
            if (u < x) then
                a=u
            else
                b=u
            end if
            if (fu <= fw .or. w == x) then
                v=w
                fv=fw
                w=u
                fw=fu
            else if (fu <= fv .or. v == x .or. v == w) then
                v=u
                fv=fu
            end if
        end if
    end do
    call nrerror('brentnew: exceed maximum iterations')
    CONTAINS
    !BL
    SUBROUTINE shft(a, b, c, d)
    REAL(DP), INTENT(OUT) :: a
    REAL(DP), INTENT(INOUT) :: b, c
    REAL(DP), INTENT(IN) :: d
    a=b
    b=c
    c=d
    END SUBROUTINE shft
    END FUNCTION brentmindist ! %>

    FUNCTION ergodic(mat, states)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ! %<
        !
        !    Computes the stationary distribution of an ergodic matrix
        !    using state space reduction. We write the transition
        !    matrix as P = [(A C) / (B D)], A being n-1xn-1, then
        !    iterate downwards. Got this from the R. Feres notes
        !
        !    Modified: 06/22/2017
        !
        !    PARAMETERS
        !
        !    mat: The Markov transition matrix for which a stationary
        !    distribution is desired.
        !
        !    states: The number of discrete states in the trans. matrix
        !    (kinda superfluous)
        !    
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IMPLICIT NONE
        INTEGER, intent(in) :: states
        INTEGER :: iter, j, j2, dim_1, dim_2
        REAL(DP), DIMENSION(states, states), intent(in) :: mat
        REAL(DP), DIMENSION(states, states) :: P
        REAL(DP) :: inverse, total
        REAL(DP), DIMENSION(states) :: ergodic

        iter  = SIZE(mat, 1)
        P = mat

        do while (iter > 1)
            dim_1 = iter - 1
            dim_2 = iter - 1
            inverse = SUM(P(iter, 1:dim_1))
            ! Storing resulting vectors in P
            P(1:dim_1, iter) = P(1:dim_1, iter)/inverse
            ! Using the censored Markov identity P_n = A + B(1-D)C
            do while (dim_2 > 0)
                P(1:dim_1, dim_2) = P(1:dim_1, dim_2) + &
                P(iter, dim_2)*P(1:dim_1, iter)
                dim_2 = dim_2 - 1
            end do
            iter = iter - 1
        end do

        ergodic(1) = 1
        j=2
        !Back out the stationary probabilities
        do while (j <= size(P, 1))
            j2 = j - 1
            !write(*,*) DOT_PRODUCT(ergodic(1:j2), P(1:j2, j))
            ergodic(j) = DOT_PRODUCT(ergodic(1:j2), P(1:j2, j))
            j = j + 1
        end do
        total = sum(ergodic)
        ergodic = ergodic/total
        RETURN

    END FUNCTION ! %>

    SUBROUTINE weightPrimel(primeChoice, nodes, weightout, primeL, primeH)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ! %<
        !
        !    Function, given a value, find the interpolation rule to
        !    obtain that value on a state grid.
        !
        !    Modified: 05/21/2019
        !
        !    PARAMETERS
        !
        !    primeChoice: The chosen value for which an interpolation rule
        !    is needed.
        !
        !    nodes: The full grid over which interpolation is needed.
        !
        !    weightout: Weight assigned to the two closest grid points
        !    for linear interpolation purposes.
        !
        !    primeL/primeH: the lower and upper grid points used for
        !    linear interpolation purposes.
        !    
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        IMPLICIT NONE
        REAL(8), INTENT(IN) :: primeChoice
        REAL(8), DIMENSION(:), INTENT(IN) :: nodes
        REAL(8), INTENT(INOUT) :: weightout
        INTEGER, INTENT(INOUT) :: primeL, primeH
        INTEGER :: GridSize, nearestNode

        Gridsize = size(nodes,1)
        nearestNode = minloc(abs(primeChoice-nodes), 1)
        if (primeChoice==nodes(nearestNode)) then
            if (nearestNode<Gridsize) then ! Is it at amin?
                primeL=nearestNode
                primeH=nearestNode+1
                if (primeH>Gridsize) then
                write(*,*) "error1", primeH, primeL, primeChoice, &
                    nearestnode, nodes(nearestNode)
                end if
                weightout=1.0
            else ! Is it at amax?
                primeL=nearestNode-1
                primeH=nearestNode
                if (primeL<1) then
                write(*,*) "error2", primeH, primeL, primeChoice, &
                    nearestnode, nodes(nearestNode)
                end if
                weightout=0.0
            end if
        else
            if (primeChoice-nodes(nearestNode)>0) then
                primel=nearestNode
                primeh=nearestNode+1
                if (primeH>Gridsize) then
                write(*,*) "error3", primeH, primeL, primeChoice, &
                    nearestnode, nodes(nearestNode)
                end if
                weightout=1-(primeChoice-nodes(primel))/(nodes(primeh)-nodes(primel))
            else
                primel=nearestNode-1
                primeh=nearestNode
                if (primeL<1) then
                write(*,*) "error4", primeH, primeL, primeChoice, &
                    nearestnode, nodes(nearestNode)
                end if
                weightout=1-(primeChoice-nodes(primel))/(nodes(primeh)-nodes(primel))
            end if
        end if
        
        if (weightout < 0.0) then
            write(*,*) "undef", nodes(primel:primeh), primeChoice, (primeChoice-nodes(primel))/(nodes(primeh)-nodes(primel))
        end if
    END SUBROUTINE weightPrimeL ! %>

    SUBROUTINE weightPrimevec(primeChoice, nodes, weightout, primeL, primeH)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ! %<
        !
        !    Vector wrapper for weightPrimel.
        !
        !    Modified: 08/09/2019
        !
        !    PARAMETERS
        !
        !    primeChoice: The chosen value for which an interpolation rule
        !    is needed.
        !
        !    nodes: The full grid over which interpolation is needed.
        !
        !    weightout: Weight assigned to the two closest grid points
        !    for linear interpolation purposes.
        !
        !    primeL/primeH: the lower and upper grid points used for
        !    linear interpolation purposes.
        !    
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        IMPLICIT NONE
        REAL(8), DIMENSION(:), INTENT(IN) :: primeChoice
        REAL(8), DIMENSION(:), INTENT(IN) :: nodes
        REAL(8), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: weightout
        INTEGER, DIMENSION(:), ALLOCATABLE,  INTENT(INOUT) :: primeL, primeH
        INTEGER :: i, choiceSize

        choiceSize = size(primeChoice,1)
        allocate(weightout(choiceSize), primeL(choiceSize), primeH(choiceSize))

        do i=1, choiceSize
            call weightPrimeL(primeChoice(i), nodes, weightout(i), primeL(i), primeH(i))
        end do

    END SUBROUTINE weightPrimevec ! %>

    REAL(DP) FUNCTION cdfnormal (x) ! %<
    !
    !*******************************************************************************
    !
    !! cdfnormal evaluates the Normal 01 CDF.
    !
    !
    !  Reference:
    !
    !    A G Adams,
    !    Areas Under the Normal Curve,
    !    Algorithm 39,
    !    Computer j.,
    !    Volume 12, pages 197-198, 1969.
    !
    !  Modified:
    !
    !    10 February 1999
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, real X, the argument of the CDF.
    !
    !    Output, real CDF, the value of the CDF.
    !
      implicit none
    !
      real, parameter :: a1 = 0.398942280444E+00
      real, parameter :: a2 = 0.399903438504E+00
      real, parameter :: a3 = 5.75885480458E+00
      real, parameter :: a4 = 29.8213557808E+00
      real, parameter :: a5 = 2.6813121679E+00
      real, parameter :: a6 = 48.6959930692E+00
      real, parameter :: a7 = 5.92885724438E+00
      real, parameter :: b0 = 0.398942280385E+00
      real, parameter :: b1 = 3.8052E-08
      real, parameter :: b2 = 1.00000615302E+00
      real, parameter :: b3 = 3.98064794E-04
      real, parameter :: b4 = 1.98615381364E+00
      real, parameter :: b5 = 0.151679116635E+00
      real, parameter :: b6 = 5.29330324926E+00
      real, parameter :: b7 = 4.8385912808E+00
      real, parameter :: b8 = 15.1508972451E+00
      real, parameter :: b9 = 0.742380924027E+00
      real, parameter :: b10 = 30.789933034E+00
      real, parameter :: b11 = 3.99019417011E+00
      !real cdfnormal
      real q
      real(DP), intent(in) :: x
      real y
    !
    !  |X| <= 1.28.
    !

      if ( abs ( x ) <= 1.28 ) then

        y = 0.5E+00 * x**2

        q = 0.5E+00 - abs ( x ) * ( a1 - a2 * y / ( y + a3 - a4 / ( y + a5 &
          + a6 / ( y + a7 ) ) ) )

    !
    !  1.28 < |X| <= 12.7
    !
      else if ( abs ( x ) <= 12.7E+00 ) then

        y = 0.5E+00 * x**2

        q = exp ( - y ) * b0 / ( abs ( x ) - b1 &
          + b2 / ( abs ( x ) + b3 &
          + b4 / ( abs ( x ) - b5 &
          + b6 / ( abs ( x ) + b7 &
          - b8 / ( abs ( x ) + b9 &
          + b10 / ( abs ( x ) + b11 ) ) ) ) ) )
    !
    !  12.7 < |X|
    !
      else

        q = 0.0E+00

      end if
    !
    !  Take account of negative X.
    !
      if ( x < 0.0E+00 ) then
        cdfnormal = q
      else
        cdfnormal = 1.0E+00 - q
      end if

      return
    end ! %>

    subroutine pol_linworking(state, aArray, DArray,cArray,vArray,choiceArray,&
        achoicelin, Dchoicelin, cchoicelin, rentallin, choicelin, welfare)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ! %<
        !
        !    Given a current state vector, this subroutine interpolates
        !    (a', D') choice and constructs the policy functions corresponding
        !    to that state.
        !
        !    Modified: 05/21/2019
        !    
        !    PARAMETERS
        !
        !    state: Vector of the individual states as well as the period
        !    during the transition process.
        !
        !    *Array: Policy arrays and adjustment indicator considered with this
        !    function.
        !
        !    *choicelin: Outputted voluntary equity, durable and consumption
        !    choices given policy arrays.
        !
        !    rentallin: Coarsening of the three-valued choice indicator
        !    into two: whether one rents or not.
        !
        !    welfare: Linearly interpolated utility (currently antiquated).
        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        IMPLICIT NONE
        REAL(DP) :: weightal, weightDl
        REAL(DP) :: nearestanode, nearestDnode
        INTEGER :: al, Dl, ah, Dh
        REAL(DP), dimension(:), INTENT(IN) :: state
        REAL(DP), dimension(:,:,:,:,:,:), INTENT(IN) :: aArray,DArray,cArray,vArray,choiceArray
        REAL(DP) :: astate, Dstate, zstate, uestate, t, hpstate,&
            achoicelin, Dchoicelin, cchoicelin, rentallin, choicelin, welfare

        ! TODO: Fix the order here along with everything in simulation file.
        zstate=state(1)
        uestate=state(2)
        Dstate=state(3)
        astate=state(4)
        t=state(5)
        hpstate=state(6)

        if (astate>amax) then
        astate=amax-1e-10
        elseif (astate<amin) then
        astate=amin+1e-10
        end if

        if (dstate>dmax) then
        dstate=dmax-1e-10
        !elseif (dstate<0.0) then
        !dstate=1e-10
        end if

        call weightprimeL(astate, anodes, weightal, al, ah)
        call weightprimeL(Dstate, Dnodes, weightDl, Dl, Dh)

        ! Explicitly coding rental indicator as one or zero instead of interpolating,
        ! experimental.
        if ((weightal>0.5) .AND. (weightDl>0.5)) then
            choicelin=choiceArray(zstate, uestate, Dh, ah, t, hpstate)
        elseif ((weightal<=0.5) .AND. (weightDl>0.5)) then
            choicelin=choiceArray(zstate, uestate, Dh, al, t, hpstate)
        elseif ((weightal>0.5) .and. (weightDl<=0.5)) then
            choicelin=choiceArray(zstate, uestate, Dl, ah, t, hpstate)
        else
            choicelin=choiceArray(zstate, uestate, Dl, al, t, hpstate)
        end if
 

        ! Still interpolating over this
        achoicelin=weightal*weightDl*aArray(zstate, uestate, Dl, al, t, hpstate)+&
            (1-weightal)*weightDl*aArray(zstate, uestate, Dl, ah, t, hpstate)+&
            weightal*(1-weightDl)*aArray(zstate, uestate, Dh, al, t, hpstate)+&
            (1-weightal)*(1-weightDl)*aArray(zstate, uestate, Dh, ah, t, hpstate)
        Dchoicelin=weightal*weightDl*DArray(zstate, uestate, Dl, al, t, hpstate)+&
            (1-weightal)*weightDl*DArray(zstate, uestate, Dl, ah, t, hpstate)+&
            weightal*(1-weightDl)*DArray(zstate, uestate, Dh, al, t, hpstate)+&
            (1-weightal)*(1-weightDl)*DArray(zstate, uestate, Dh, ah, t, hpstate)
        cchoicelin=weightal*weightDl*cArray(zstate, uestate, Dl, al, t, hpstate)+&
            (1-weightal)*weightDl*cArray(zstate, uestate, Dl, ah, t, hpstate)+&
            weightal*(1-weightDl)*cArray(zstate, uestate, Dh, al, t, hpstate)+&
            (1-weightal)*(1-weightDl)*cArray(zstate, uestate, Dh, ah, t, hpstate)
        welfare=weightal*weightDl*vArray(zstate, uestate, Dl, al, t, hpstate)+&
           (1-weightal)*weightDl*vArray(zstate, uestate, Dl, ah, t, hpstate)+&
           weightal*(1-weightDl)*vArray(zstate, uestate, Dh, al, t, hpstate)+&
           (1-weightal)*(1-weightDl)*vArray(zstate, uestate, Dh, ah, t, hpstate)

        ! Interpolation check: if D', D clearly different, consider adjustment
        ! to have taken place.
        if ((ABS(Dchoicelin/Dnodes(Dh) - 1) > 0.2) .AND. (choicelin == 2)) choicelin=1
 
        if (choicelin < 3) then
            rentallin=0
            if (Dchoicelin <= Dmin .AND. choicelin == 1) Dchoicelin = Dmin
        else
            rentallin=1
        end if

       if (achoicelin < borrowconstraint-1e-3) then
           write(*,*) "error"
           write(*,*) state
           write(*,*) aArray(zstate, uestate, Dl, al, t, hpstate), aArray(zstate, uestate, Dh, al, t, hpstate), aArray(zstate, uestate, Dl, ah, t, hpstate), aArray(zstate, uestate, Dh, ah, t, hpstate)
           write(*,*) choiceArray(zstate, uestate, Dl, al, t, hpstate), choiceArray(zstate, uestate, Dh, al, t, hpstate), choiceArray(zstate, uestate, Dl, ah, t, hpstate), choiceArray(zstate, uestate, Dh, ah, t, hpstate)
           write(*,*) dArray(zstate, uestate, Dl, al, t, hpstate), dArray(zstate, uestate, Dh, al, t, hpstate), dArray(zstate, uestate, Dl, ah, t, hpstate), dArray(zstate, uestate, Dh, ah, t, hpstate)
           write(*,*) weightal, weightdl, achoicelin, dchoicelin, choicelin
           call EXIT(1)
       else if (achoicelin < 0) then
           achoicelin = 0.0
       end if
    end subroutine pol_linworking ! %>


end module lifecycle_algs
