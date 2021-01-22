module lifecycle_algs
    USE share
    USE nrtype
    USE nrutil, ONLY : assert_eq, imaxloc, iminloc, nrerror, swap
    USE lifeCycle_vfuncs
    USE OMP_LIB
    IMPLICIT NONE

    CONTAINS

    SUBROUTINE amoeba_price(p, y, ftol, numHH, simSteady)
        IMPLICIT NONE

        INTEGER(I4B) :: iter
        REAL(DP), INTENT(IN) :: ftol
        REAL(DP), DIMENSION(2), INTENT(INOUT) :: y
        REAL(DP), DIMENSION(2,1), INTENT(INOUT) :: p
        REAL(DP), DIMENSION(2) :: yholder
        REAL(DP), DIMENSION(2,1) :: pholder
        INTEGER(I4B), INTENT(IN) :: numHH
        INTEGER(I4B), PARAMETER :: ITMAX=5000000
        REAL(DP), PARAMETER :: TINY=1.0e-10
        INTEGER(I4B) :: ihi, ndim
        REAL(DP), DIMENSION(size(p,2)) :: psum
        REAL(DP), EXTERNAL :: simSteady
        call amoeba_private
        CONTAINS
    !BL
        SUBROUTINE amoeba_private
        IMPLICIT NONE
        INTEGER(I4B) :: i, ilo, inhi
        REAL(DP) :: rtol, ysave, ytry, ytmp
        ndim=assert_eq(size(p, 2), size(p, 1)-1, size(y)-1,'amoeba')
        iter=0
        psum(:)=sum(p(:,1))

        pholder=p
        yholder=y

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
            write(*,*) "p", p
            write(*,*) "pholder", pholder
            write(*,*) "y", y
            write(*,*) "yholder", yholder
            call nrerror('ITMAX exceeded in amoeba')
            end if

            ytry=amotry(-1.0_dp)
            iter=iter+1
            if (ytry <= y(ilo)) then

                ytry=amotry(2.0_dp)
                iter=iter+1
            else if (ytry >= y(inhi)) then
                ysave=y(ihi)


                ytry=amotry(0.5_dp)
                iter=iter+1
                if (ytry >= ysave) then
                    p(:,:)=0.5_dp*(p(:,:)+spread(p(ilo,:), 1, size(p, 1)))
                    do i=1, ndim+1


                        if (i /= ilo) y(i)=simSteady(p(i,:),numHH)
                        !write(*,*) y(i)
                    end do
                    iter=iter+ndim
                    psum(:)=sum(p(:,1))
                end if
            end if
        end do
        END SUBROUTINE amoeba_private
    !BL
        FUNCTION amotry(fac)
        IMPLICIT NONE
        REAL(DP), INTENT(IN) :: fac
        REAL(DP) :: amotry
        REAL(DP) :: fac1, fac2, ytry
        REAL(DP), DIMENSION(size(p,2)) :: ptry
        fac1=(1.0_dp-fac)/ndim
        fac2=fac1-fac
        ptry(:)=psum(:)*fac1-p(ihi,1)*fac2
        ytry=simSteady(ptry,numHH)
        !write(*,*) "done?"
        if (ytry < y(ihi)) then
            y(ihi)=ytry
            psum(:)=psum(:)-p(ihi,1)+ptry(:)
            p(ihi,1)=ptry(1)
        end if
        amotry=ytry
        END FUNCTION amotry
    END SUBROUTINE amoeba_price




    SUBROUTINE amoeba(stateholder, p, y, ftol, valfuncadjust)
        IMPLICIT NONE

        INTEGER(I4B) :: iter
        REAL(DP), INTENT(IN) :: ftol
        REAL(DP), DIMENSION(3), INTENT(INOUT) :: y
        REAL(DP), DIMENSION(3, 2), INTENT(INOUT) :: p
        REAL(DP), DIMENSION(3) :: yholder
        REAL(DP), DIMENSION(3, 2) :: pholder
        REAL(DP), DIMENSION(5), INTENT(IN) :: stateholder
        INTEGER(I4B), PARAMETER :: ITMAX=5000000
        REAL(DP), PARAMETER :: TINY=1.0e-10
        INTEGER(I4B) :: ihi, ndim
        REAL(DP), DIMENSION(size(p, 2)) :: psum
        REAL(DP), EXTERNAL :: valfuncadjust
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

                        if (i /= ilo) y(i)=valfuncadjust(p(i,:), stateholder)
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
        REAL(DP), DIMENSION(5) :: stateholder
        REAL(DP), DIMENSION(size(p, 2)) :: ptry
        fac1=(1.0_dp-fac)/ndim
        fac2=fac1-fac
        ptry(:)=psum(:)*fac1-p(ihi,:)*fac2
        !write(*,*) "private", stateholder
        ytry=valfuncadjust(ptry, stateholder)
        !write(*,*) "done?"
        if (ytry < y(ihi)) then
            y(ihi)=ytry
            psum(:)=psum(:)-p(ihi,:)+ptry(:)
            p(ihi,:)=ptry(:)
        end if
        amotry=ytry
        END FUNCTION amotry
    END SUBROUTINE amoeba





    FUNCTION brentnew(ax, bx, cx, value, tol, aholder, Dholder, zholder, hpholder, tholder, xmin)
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: ax, bx, cx, tol, aholder, Dholder, zholder, hpholder, tholder
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
    fx=value(x, aholder, Dholder, zholder, hpholder, tholder)
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
        fu=value(u, aholder, Dholder, zholder, hpholder, tholder)
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
    END FUNCTION brentnew




    FUNCTION brentmindist(ax, bx, cx, value, numhh, tol, xmin)
    IMPLICIT NONE
    REAL(DP), INTENT(IN) :: ax, bx, cx, tol
    INTEGER :: numhh
    REAL(DP), INTENT(OUT) :: xmin
    REAL(DP) :: brentmindist
    EXTERNAL :: value
    REAL(DP) :: value



    INTEGER(I4B), PARAMETER :: ITMAX=200
    REAL(DP), PARAMETER :: CGOLD=0.3819660_dp, ZEPS=1.0e-3_dp*epsilon(ax)
    INTEGER(I4B) :: iter
    REAL(DP) :: a, b, d, e, etemp, fu, fv, fw, fx, p, q, r2, tol1, tol2, u, v, w, x, xm
    a=min(ax, cx)
    b=max(ax, cx)
    v=bx
    w=v
    x=v
    e=0.0
    fx=value(x, numhh)
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
        fu=value(u, numhh)
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
    END FUNCTION brentmindist

    ! Computes the stationary distribution of an ergodic matrix
    ! using state space reduction. We write the transition
    ! matrix as P = [(A C) / (B D)], A being n-1xn-1, then
    ! iterate downwards. Got this from the R. Feres notes
    FUNCTION ergodic(mat, states)
        IMPLICIT NONE
        INTEGER, intent(in) :: states
        INTEGER :: iter, j, j2, dim1, dim2
        REAL(DP), DIMENSION(states, states), intent(in) :: mat
        REAL(DP), DIMENSION(states, states) :: P
        REAL(DP) :: inverse, total
        REAL(DP), DIMENSION(states) :: ergodic

        iter  = SIZE(mat, 1)
        P = mat

        do while (iter > 1)
            dim1 = iter - 1
            dim2 = iter - 1
            inverse = SUM(P(iter, 1:dim1))
            ! Storing resulting vectors in P
            P(1:dim1, iter) = P(1:dim1, iter)/inverse
            ! Using the censored Markov identity P_n = A + B(1-D)C
            do while (dim2 > 0)
                P(1:dim1, dim2) = P(1:dim1, dim2) + &
                P(iter, dim2)*P(1:dim1, iter)
                dim2 = dim2 - 1
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

    END FUNCTION



       REAL(DP) FUNCTION cdfnormal (x)
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
    end

    subroutine pol_linworking(astate, Dstate, zstate, hpstate, t, achoicelin, Dchoicelin, cchoicelin, rentallin, choicelin, welfare)
        IMPLICIT NONE
        REAL(DP) :: weightal, weightDl
        REAL(DP) :: nearestanode, nearestDnode
        REAL(DP) :: al, ah, Dl, Dh
        REAL(DP) :: astate, Dstate, zstate, hpstate, t, achoicelin, Dchoicelin, cchoicelin, rentallin, choicelin, welfare

        if (astate>amax) then
        astate=amax-.00000000001
        elseif (astate<amin) then
        astate=amin+.000000000001
        end if

        nearestanode=minloc(abs(astate-anodes), 1)
        if (astate==anodes(nearestanode)) then
            if (nearestanode<agridsize) then
                al=nearestanode
                ah=nearestanode+1
                weightal=1.0
            else
                al=nearestanode-1
                ah=nearestanode
                weightal=0.0
            end if
        else
            if (astate-anodes(nearestanode)>0) then
                al=nearestanode
                ah=nearestanode+1
                if (ah>agridsize) then
                write(*,*) "error", ah, al, astate, nearestanode, anodes(nearestanode)
                end if
                weightal=1-(astate-anodes(al))/(anodes(ah)-anodes(al))
            else
                al=nearestanode-1
                if (al<1) then
                write(*,*) "error", ah, al, astate, nearestanode, anodes(nearestanode)
                end if
                ah=nearestanode
                weightal=1-(astate-anodes(al))/(anodes(ah)-anodes(al))
            end if
        end if


        nearestDnode=minloc(abs(Dstate-Dnodes), 1)
        if (Dstate==Dnodes(nearestDnode)) then
            if (nearestDnode<Dgridsize) then
                Dl=nearestDnode
                Dh=nearestDnode+1
                weightDl=1.0
            else
                Dl=nearestDnode-1
                Dh=nearestDnode
                weightDl=0.0
            end if
        else
            if (Dstate-Dnodes(nearestDnode)>0) then
                Dl=nearestDnode
                Dh=nearestDnode+1
                weightDl=1-(Dstate-Dnodes(Dl))/(Dnodes(Dh)-Dnodes(Dl))
            else
                Dl=nearestDnode-1
                Dh=nearestDnode
                weightDl=1-(Dstate-Dnodes(Dl))/(Dnodes(Dh)-Dnodes(Dl))
            end if
        end if

        !write(*,*) weightal, weightDl, al, Dl, estate, Dchoice(al, Dl, estate)

        !write(*,*) "r  ", "r  ", al, ah, Dl, Dh, zstate, epsstate, t

       ! Explicitly coding rental indicator as one or zero instead of interpolating,
       ! experimental.
        if ((weightal>0.5) .AND. (weightDl>0.5)) then
            choicelin=choiceindicator(ah, Dh, zstate, hpstate, t)

            if (choiceindicator(ah, Dh, zstate, hpstate, t)==1) then
                !achoicelin=weightal*weightDl*achoiceadjust(al, Dl, zstate, hpstate, t)+(1-weightal)*weightDl*achoiceadjust(ah, Dl, zstate, hpstate, t)+weightal*(1-weightDl)*achoiceadjust(al, Dh, zstate, hpstate, t)+(1-weightal)*(1-weightDl)*achoiceadjust(ah, Dh, zstate, hpstate, t)
                !Dchoicelin=weightal*weightDl*Dchoiceadjust(al, Dl, zstate, hpstate, t)+(1-weightal)*weightDl*Dchoiceadjust(ah, Dl, zstate, hpstate, t)+weightal*(1-weightDl)*Dchoiceadjust(al, Dh, zstate, hpstate, t)+(1-weightal)*(1-weightDl)*Dchoiceadjust(ah, Dh, zstate, hpstate, t)
                !cchoicelin=weightal*weightDl*cchoiceadjust(al, Dl, zstate, hpstate, t)+(1-weightal)*weightDl*cchoiceadjust(ah, Dl, zstate, hpstate, t)+weightal*(1-weightDl)*cchoiceadjust(al, Dh, zstate, hpstate, t)+(1-weightal)*(1-weightDl)*cchoiceadjust(ah, Dh, zstate, hpstate, t)
                rentallin=0

            elseif  (choiceindicator(ah, Dh, zstate, hpstate, t)==2) then
                !achoicelin=weightal*weightDl*achoicenoadjust(al, Dl, zstate, hpstate, t)+(1-weightal)*weightDl*achoicenoadjust(ah, Dl, zstate, hpstate, t)+weightal*(1-weightDl)*achoicenoadjust(al, Dh, zstate, hpstate, t)+(1-weightal)*(1-weightDl)*achoicenoadjust(ah, Dh, zstate, hpstate, t)
                !Dchoicelin=weightal*weightDl*Dchoicenoadjust(al, Dl, zstate, hpstate, t)+(1-weightal)*weightDl*Dchoicenoadjust(ah, Dl, zstate, hpstate, t)+weightal*(1-weightDl)*Dchoicenoadjust(al, Dh, zstate, hpstate, t)+(1-weightal)*(1-weightDl)*Dchoicenoadjust(ah, Dh, zstate, hpstate, t)
                !cchoicelin=weightal*weightDl*cchoicenoadjust(al, Dl, zstate, hpstate, t)+(1-weightal)*weightDl*cchoicenoadjust(ah, Dl, zstate, hpstate, t)+weightal*(1-weightDl)*cchoicenoadjust(al, Dh, zstate, hpstate, t)+(1-weightal)*(1-weightDl)*cchoicenoadjust(ah, Dh, zstate, hpstate, t)
                rentallin=0

            else
                !achoicelin=weightal*weightDl*achoicerent(al, Dl, zstate, hpstate, t)+(1-weightal)*weightDl*achoicerent(ah, Dl, zstate, hpstate, t)+weightal*(1-weightDl)*achoicerent(al, Dh, zstate, hpstate, t)+(1-weightal)*(1-weightDl)*achoicerent(ah, Dh, zstate, hpstate, t)
                !Dchoicelin=weightal*weightDl*Dchoicerent(al, Dl, zstate, hpstate, t)+(1-weightal)*weightDl*Dchoicerent(ah, Dl, zstate, hpstate, t)+weightal*(1-weightDl)*Dchoicerent(al, Dh, zstate, hpstate, t)+(1-weightal)*(1-weightDl)*Dchoicerent(ah, Dh, zstate, hpstate, t)
                !cchoicelin=weightal*weightDl*cchoicerent(al, Dl, zstate, hpstate, t)+(1-weightal)*weightDl*cchoicerent(ah, Dl, zstate, hpstate, t)+weightal*(1-weightDl)*cchoicerent(al, Dh, zstate, hpstate, t)+(1-weightal)*(1-weightDl)*cchoicerent(ah, Dh, zstate, hpstate, t)
                rentallin=1
            endif

        elseif ((weightal<=0.5) .AND. (weightDl>0.5)) then
            choicelin=choiceindicator(al, Dh, zstate, hpstate, t)

            if (choiceindicator(al, Dh, zstate, hpstate, t)==1) then
                !achoicelin=weightal*weightDl*achoiceadjust(al, Dl, zstate, hpstate, t)+(1-weightal)*weightDl*achoiceadjust(ah, Dl, zstate, hpstate, t)+weightal*(1-weightDl)*achoiceadjust(al, Dh, zstate, hpstate, t)+(1-weightal)*(1-weightDl)*achoiceadjust(ah, Dh, zstate, hpstate, t)
                !Dchoicelin=weightal*weightDl*Dchoiceadjust(al, Dl, zstate, hpstate, t)+(1-weightal)*weightDl*Dchoiceadjust(ah, Dl, zstate, hpstate, t)+weightal*(1-weightDl)*Dchoiceadjust(al, Dh, zstate, hpstate, t)+(1-weightal)*(1-weightDl)*Dchoiceadjust(ah, Dh, zstate, hpstate, t)
                !cchoicelin=weightal*weightDl*cchoiceadjust(al, Dl, zstate, hpstate, t)+(1-weightal)*weightDl*cchoiceadjust(ah, Dl, zstate, hpstate, t)+weightal*(1-weightDl)*cchoiceadjust(al, Dh, zstate, hpstate, t)+(1-weightal)*(1-weightDl)*cchoiceadjust(ah, Dh, zstate, hpstate, t)
                rentallin=0

            elseif  (choiceindicator(al, Dh, zstate, hpstate, t)==2) then
                !achoicelin=weightal*weightDl*achoicenoadjust(al, Dl, zstate, hpstate, t)+(1-weightal)*weightDl*achoicenoadjust(ah, Dl, zstate, hpstate, t)+weightal*(1-weightDl)*achoicenoadjust(al, Dh, zstate, hpstate, t)+(1-weightal)*(1-weightDl)*achoicenoadjust(ah, Dh, zstate, hpstate, t)
                !Dchoicelin=weightal*weightDl*Dchoicenoadjust(al, Dl, zstate, hpstate, t)+(1-weightal)*weightDl*Dchoicenoadjust(ah, Dl, zstate, hpstate, t)+weightal*(1-weightDl)*Dchoicenoadjust(al, Dh, zstate, hpstate, t)+(1-weightal)*(1-weightDl)*Dchoicenoadjust(ah, Dh, zstate, hpstate, t)
                !cchoicelin=weightal*weightDl*cchoicenoadjust(al, Dl, zstate, hpstate, t)+(1-weightal)*weightDl*cchoicenoadjust(ah, Dl, zstate, hpstate, t)+weightal*(1-weightDl)*cchoicenoadjust(al, Dh, zstate, hpstate, t)+(1-weightal)*(1-weightDl)*cchoicenoadjust(ah, Dh, zstate, hpstate, t)
                rentallin=0

            else
                !achoicelin=weightal*weightDl*achoicerent(al, Dl, zstate, hpstate, t)+(1-weightal)*weightDl*achoicerent(ah, Dl, zstate, hpstate, t)+weightal*(1-weightDl)*achoicerent(al, Dh, zstate, hpstate, t)+(1-weightal)*(1-weightDl)*achoicerent(ah, Dh, zstate, hpstate, t)
                !Dchoicelin=weightal*weightDl*Dchoicerent(al, Dl, zstate, hpstate, t)+(1-weightal)*weightDl*Dchoicerent(ah, Dl, zstate, hpstate, t)+weightal*(1-weightDl)*Dchoicerent(al, Dh, zstate, hpstate, t)+(1-weightal)*(1-weightDl)*Dchoicerent(ah, Dh, zstate, hpstate, t)
                !cchoicelin=weightal*weightDl*cchoicerent(al, Dl, zstate, hpstate, t)+(1-weightal)*weightDl*cchoicerent(ah, Dl, zstate, hpstate, t)+weightal*(1-weightDl)*cchoicerent(al, Dh, zstate, hpstate, t)+(1-weightal)*(1-weightDl)*cchoicerent(ah, Dh, zstate, hpstate, t)
                rentallin=1
            endif

        elseif ((weightal>0.5) .and. (weightDl<=0.5)) then
            choicelin=choiceindicator(ah, Dl, zstate, hpstate, t)

            if (choiceindicator(ah, Dl, zstate, hpstate, t)==1) then
                !achoicelin=weightal*weightDl*achoiceadjust(al, Dl, zstate, hpstate, t)+(1-weightal)*weightDl*achoiceadjust(ah, Dl, zstate, hpstate, t)+weightal*(1-weightDl)*achoiceadjust(al, Dh, zstate, hpstate, t)+(1-weightal)*(1-weightDl)*achoiceadjust(ah, Dh, zstate, hpstate, t)
                !Dchoicelin=weightal*weightDl*Dchoiceadjust(al, Dl, zstate, hpstate, t)+(1-weightal)*weightDl*Dchoiceadjust(ah, Dl, zstate, hpstate, t)+weightal*(1-weightDl)*Dchoiceadjust(al, Dh, zstate, hpstate, t)+(1-weightal)*(1-weightDl)*Dchoiceadjust(ah, Dh, zstate, hpstate, t)
                !cchoicelin=weightal*weightDl*cchoiceadjust(al, Dl, zstate, hpstate, t)+(1-weightal)*weightDl*cchoiceadjust(ah, Dl, zstate, hpstate, t)+weightal*(1-weightDl)*cchoiceadjust(al, Dh, zstate, hpstate, t)+(1-weightal)*(1-weightDl)*cchoiceadjust(ah, Dh, zstate, hpstate, t)
                rentallin=0

            elseif  (choiceindicator(ah, Dl, zstate, hpstate, t)==2) then
                !achoicelin=weightal*weightDl*achoicenoadjust(al, Dl, zstate, hpstate, t)+(1-weightal)*weightDl*achoicenoadjust(ah, Dl, zstate, hpstate, t)+weightal*(1-weightDl)*achoicenoadjust(al, Dh, zstate, hpstate, t)+(1-weightal)*(1-weightDl)*achoicenoadjust(ah, Dh, zstate, hpstate, t)
                !Dchoicelin=weightal*weightDl*Dchoicenoadjust(al, Dl, zstate, hpstate, t)+(1-weightal)*weightDl*Dchoicenoadjust(ah, Dl, zstate, hpstate, t)+weightal*(1-weightDl)*Dchoicenoadjust(al, Dh, zstate, hpstate, t)+(1-weightal)*(1-weightDl)*Dchoicenoadjust(ah, Dh, zstate, hpstate, t)
                !cchoicelin=weightal*weightDl*cchoicenoadjust(al, Dl, zstate, hpstate, t)+(1-weightal)*weightDl*cchoicenoadjust(ah, Dl, zstate, hpstate, t)+weightal*(1-weightDl)*cchoicenoadjust(al, Dh, zstate, hpstate, t)+(1-weightal)*(1-weightDl)*cchoicenoadjust(ah, Dh, zstate, hpstate, t)
                rentallin=0

            else
                !achoicelin=weightal*weightDl*achoicerent(al, Dl, zstate, hpstate, t)+(1-weightal)*weightDl*achoicerent(ah, Dl, zstate, hpstate, t)+weightal*(1-weightDl)*achoicerent(al, Dh, zstate, hpstate, t)+(1-weightal)*(1-weightDl)*achoicerent(ah, Dh, zstate, hpstate, t)
                !Dchoicelin=weightal*weightDl*Dchoicerent(al, Dl, zstate, hpstate, t)+(1-weightal)*weightDl*Dchoicerent(ah, Dl, zstate, hpstate, t)+weightal*(1-weightDl)*Dchoicerent(al, Dh, zstate, hpstate, t)+(1-weightal)*(1-weightDl)*Dchoicerent(ah, Dh, zstate, hpstate, t)
                !cchoicelin=weightal*weightDl*cchoicerent(al, Dl, zstate, hpstate, t)+(1-weightal)*weightDl*cchoicerent(ah, Dl, zstate, hpstate, t)+weightal*(1-weightDl)*cchoicerent(al, Dh, zstate, hpstate, t)+(1-weightal)*(1-weightDl)*cchoicerent(ah, Dh, zstate, hpstate, t)
                rentallin=1
            endif

        else
            choicelin=choiceindicator(al, Dl, zstate, hpstate, t)

            if (choiceindicator(al, Dl, zstate, hpstate, t)==1) then
                !achoicelin=weightal*weightDl*achoiceadjust(al, Dl, zstate, hpstate, t)+(1-weightal)*weightDl*achoiceadjust(ah, Dl, zstate, hpstate, t)+weightal*(1-weightDl)*achoiceadjust(al, Dh, zstate, hpstate, t)+(1-weightal)*(1-weightDl)*achoiceadjust(ah, Dh, zstate, hpstate, t)
                !Dchoicelin=weightal*weightDl*Dchoiceadjust(al, Dl, zstate, hpstate, t)+(1-weightal)*weightDl*Dchoiceadjust(ah, Dl, zstate, hpstate, t)+weightal*(1-weightDl)*Dchoiceadjust(al, Dh, zstate, hpstate, t)+(1-weightal)*(1-weightDl)*Dchoiceadjust(ah, Dh, zstate, hpstate, t)
                !cchoicelin=weightal*weightDl*cchoiceadjust(al, Dl, zstate, hpstate, t)+(1-weightal)*weightDl*cchoiceadjust(ah, Dl, zstate, hpstate, t)+weightal*(1-weightDl)*cchoiceadjust(al, Dh, zstate, hpstate, t)+(1-weightal)*(1-weightDl)*cchoiceadjust(ah, Dh, zstate, hpstate, t)
                rentallin=0

            elseif  (choiceindicator(al, Dl, zstate, hpstate, t)==2) then
                !achoicelin=weightal*weightDl*achoicenoadjust(al, Dl, zstate, hpstate, t)+(1-weightal)*weightDl*achoicenoadjust(ah, Dl, zstate, hpstate, t)+weightal*(1-weightDl)*achoicenoadjust(al, Dh, zstate, hpstate, t)+(1-weightal)*(1-weightDl)*achoicenoadjust(ah, Dh, zstate, hpstate, t)
                !Dchoicelin=weightal*weightDl*Dchoicenoadjust(al, Dl, zstate, hpstate, t)+(1-weightal)*weightDl*Dchoicenoadjust(ah, Dl, zstate, hpstate, t)+weightal*(1-weightDl)*Dchoicenoadjust(al, Dh, zstate, hpstate, t)+(1-weightal)*(1-weightDl)*Dchoicenoadjust(ah, Dh, zstate, hpstate, t)
                !cchoicelin=weightal*weightDl*cchoicenoadjust(al, Dl, zstate, hpstate, t)+(1-weightal)*weightDl*cchoicenoadjust(ah, Dl, zstate, hpstate, t)+weightal*(1-weightDl)*cchoicenoadjust(al, Dh, zstate, hpstate, t)+(1-weightal)*(1-weightDl)*cchoicenoadjust(ah, Dh, zstate, hpstate, t)
                rentallin=0

            else
                !achoicelin=weightal*weightDl*achoicerent(al, Dl, zstate, hpstate, t)+(1-weightal)*weightDl*achoicerent(ah, Dl, zstate, hpstate, t)+weightal*(1-weightDl)*achoicerent(al, Dh, zstate, hpstate, t)+(1-weightal)*(1-weightDl)*achoicerent(ah, Dh, zstate, hpstate, t)
                !Dchoicelin=weightal*weightDl*Dchoicerent(al, Dl, zstate, hpstate, t)+(1-weightal)*weightDl*Dchoicerent(ah, Dl, zstate, hpstate, t)+weightal*(1-weightDl)*Dchoicerent(al, Dh, zstate, hpstate, t)+(1-weightal)*(1-weightDl)*Dchoicerent(ah, Dh, zstate, hpstate, t)
                !cchoicelin=weightal*weightDl*cchoicerent(al, Dl, zstate, hpstate, t)+(1-weightal)*weightDl*cchoicerent(ah, Dl, zstate, hpstate, t)+weightal*(1-weightDl)*cchoicerent(al, Dh, zstate, hpstate, t)+(1-weightal)*(1-weightDl)*cchoicerent(ah, Dh, zstate, hpstate, t)
                rentallin=1
            endif

        endif


        achoicelin=weightal*weightDl*achoice(al, Dl, zstate, hpstate, t)+(1-weightal)*weightDl*achoice(ah, Dl, zstate, hpstate, t)+weightal*(1-weightDl)*achoice(al, Dh, zstate, hpstate, t)+(1-weightal)*(1-weightDl)*achoice(ah, Dh, zstate, hpstate, t)
        Dchoicelin=weightal*weightDl*Dchoice(al, Dl, zstate, hpstate, t)+(1-weightal)*weightDl*Dchoice(ah, Dl, zstate, hpstate, t)+weightal*(1-weightDl)*Dchoice(al, Dh, zstate, hpstate, t)+(1-weightal)*(1-weightDl)*Dchoice(ah, Dh, zstate, hpstate, t)
        cchoicelin=weightal*weightDl*cchoice(al, Dl, zstate, hpstate, t)+(1-weightal)*weightDl*cchoice(ah, Dl, zstate, hpstate, t)+weightal*(1-weightDl)*cchoice(al, Dh, zstate, hpstate, t)+(1-weightal)*(1-weightDl)*cchoice(ah, Dh, zstate, hpstate, t)
       ! rentallin=weightal*weightDl*rentalindicator(al, Dl, zstate, hpstate, t)+(1-weightal)*weightDl*rentalindicator(ah, Dl, zstate, hpstate, t)+weightal*(1-weightDl)*rentalindicator(al, Dh, zstate, hpstate, t)+(1-weightal)*(1-weightDl)*rentalindicator(ah, Dh, zstate, hpstate, t)
        welfare=weightal*weightDl*EV(al, Dl, zstate, hpstate, t)+(1-weightal)*weightDl*EV(ah, Dl, zstate, hpstate, t)+weightal*(1-weightDl)*EV(al, Dh, zstate, hpstate, t)+(1-weightal)*(1-weightDl)*EV(ah, Dh, zstate, hpstate, t)

    end subroutine pol_linworking


end module lifecycle_algs
