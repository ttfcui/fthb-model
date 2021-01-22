module lifecycle_calibrate
    USE share
    USE lifecycle_algs
    IMPLICIT NONE

    CONTAINS

    subroutine income_array
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ! %<
        !
        !    Master subroutine for discretizing the AR(1), generating
        !    social security payments and recording the set of possible
        !    income realizations.
        !
        !    Modified: 10/23/2019
        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        IMPLICIT NONE
        REAL(8), DIMENSION(zgridsize) :: predictedlifetimeearnings, retirementincome
        INTEGER :: i, j, k, t, bet_stat
        REAL(8) :: cdfout, cdfinv
        REAL(8) :: averagelifetimeearnings, alpha_init, beta_init, beta, bound_init
        REAL(8), DIMENSION(3) :: parammodel, paramcheck
        CHARACTER(4) :: rhostr, sigmastr
        CHARACTER(255) :: matlabcall
   

      ! call tauchen
        call rouwenhorst

        do j=1, zgridsize
            Probz(j,:)=Probz(j,:)/sum(Probz(j,:))
        end do
        do j=1, zgridsize
            do k=1, zgridsize
                Probzcum(j, k)=sum(Probz(j, 1:k))
            end do
        end do

    ! Following line pops out stationary distribution.
        if (poismean <= 0) then
            Probinit = ergodic(Probz(1:zgridsize,1:zgridsize), zgridsize)
            do j=1, zgridsize-1
                Probinitcum(j) = SUM(Probinit(1:j))
            end do 
            Probinitcum(zgridsize) = 1.0
            write(0,*) NEW_LINE('A') // "// Stationary distribution grid:"
        else
        ! This is *not* the stat. dist. but instead a parametric distribution
        ! estimated over the young adult SCF income dist (ages 22 to 25)
        ! TODO: This - 2 stuff is to reduce the initial variation in the
        ! income process. Smarter way to do this?

            ! alpha, beta are corresponding shape parameters for the beta distribution
            OPEN(30, FILE='input_data/init_incomes.txt', STATUS='OLD')
            READ(30,*) alpha_init, beta_init
            write(0,*) NEW_LINE('A')
            write(0,*) alpha_init, beta_init
            Probinitcum(zgridsize-4:zgridsize) = 1.0
            do j=zgridsize-5,1,-1
                call cdfbet(1, cdfout, cdfinv, exp(znodes(j))/exp(znodes(zgridsize-5)), &
                            1.0-exp(znodes(j))/exp(znodes(zgridsize-5)), &
                            alpha_init, beta_init, bet_stat, bound_init)
                !call gamma_inc(j+1.0, poismean, cdfinv, cdfout, 1)
                write(0,*) cdfout
                Probinitcum(j) = cdfout
                if (j < zgridsize) Probinit(j+1) = Probinitcum(j+1) - Probinit(j)
            end do
        end if
        write(0,'(8F16.10)') Probinitcum

        ! Chi/permanent portion of earnings proccess
        ! This is a cubic of log income over age from the
        ! Blundell, Pistaferri, Preston data, normalized;
        ! See BPP_calibration.do for more details
        if (no_lifeinc) then
            ageearnings=0
        else
            OPEN(31, FILE='input_data/ageearnings.txt', STATUS='OLD')
            do i=1,TRetire
                READ(31, '(F16.10)') ageearnings(i)
            end do
        end if

        ! death probabilities (starting at age 60 or 65)
        if (no_death) then
            deathprob=0
        else
            OPEN(32, FILE='input_data/deathprob.txt', STATUS='OLD')
            do i=1,TDie-Tretire
                READ(32, '(F5.3)') deathprob(i)
            end do
        end if

        if (no_utilscale) then
            utilscale=1
        else
            OPEN(33, FILE='input_data/utilscale.txt', STATUS='OLD')
            do i=1,TDie
                READ(33, '(F16.10)') utilscale(i)
            end do
        end if

        ! Unemployment probabilities
        ALLOCATE(empprob(Tretire+1))
        if (no_unemp) then
            empprob = 1.0
        else
            OPEN(36, FILE='input_data/empprobs.txt', STATUS='OLD')
            do i=1,CEILING(Tretire/5.0)
                READ(36, '(F16.10)') empprob(5*(i-1)+1)
                empprob(5*(i-1)+2:MIN(5*i, Tretire)) = empprob(5*(i-1)+1)
            end do
        end if
        empprob(Tretire+1) = 1.0
        write(0,*) NEW_LINE('A') // "// Employment probabilities:"
        write(0,*) empprob

        ! Stochastic trans. cost shocks (uniform version)
        if (transgridsize == 1) then
            transcst(1) = F
            transcstProbs = 1.0
        else
        transcstProbs = 1.0/transgridsize
        do j=1, transgridsize
            ! The "0.5" means we're taking midpoints of bins over support
            transcst(j) = F*(j-0.5)/transgridsize
        end do
        end if

        !retirement regression needs to first run matlab program simulateearningsprocess
        write(rhostr, '(F4.2)') rho_z
        write(sigmastr, '(F4.2)') sigma_z
        matlabcall = 'matlab -nodisplay -desktop -r "cd ../matlab; sigma_z=' // &
            sigmastr // ';rho_z=' // rhostr // ';simulateearningsprocess;&
            exit;exit;"'
        write(*,*) matlabcall
        CALL SYSTEM(trim(matlabcall))
        averagelifetimeearnings=0.0
        OPEN(34, FILE='input_data/lifeearn.txt', STATUS='OLD')
        READ(34,*) paramCheck(:), beta
        parammodel = (/ rho_z, sigma_z, REAL(Tretire) /)
        if (maxval(abs(parammodel - paramcheck)) >= 1e-7) THEN
            write(*,*) parammodel
            write(*,*) paramCheck
            write(*,*) 'ERROR: input parameters of simulateearningsprocess.m&
                       & does not match those of current model.'
            call EXIT(3)
        end if
        predictedlifetimeearnings=exp(beta*znodes)/exp(averagelifetimeearnings)
          
        ! Other economy-wide parameters during policy period
        if (agg_policies) then
            OPEN(37, FILE='input_data/polparams.txt', STATUS='OLD')
            do i=1,size(polParam, 1)
                READ(37, '(I6.2, F16.10)') polParam(i), polRead(i)
            end do
            agg_policies = .FALSE.
        end if

        write(0,*) NEW_LINE('A') // "// Incomes post retirement (as a f'n of Z):"
        do j=1, zgridsize
        ! This function is from Guvenen and Smith 2014
            if (predictedlifetimeearnings(j)<=0.3) then
                retirementincome(j)=0.9*predictedlifetimeearnings(j)*exp(averagelifetimeearnings)
            elseif (predictedlifetimeearnings(j)>0.3 .and. predictedlifetimeearnings(j)<=2) then
                retirementincome(j)=(0.27+0.32*(predictedlifetimeearnings(j)-0.3))*exp(averagelifetimeearnings)
            elseif (predictedlifetimeearnings(j)>2 .and. predictedlifetimeearnings(j)<=4.1) then
                retirementincome(j)=(0.81+0.15*(predictedlifetimeearnings(j)-2))*exp(averagelifetimeearnings)
            else
                retirementincome(j)=1.13*exp(averagelifetimeearnings)
            end if

                write(0,*) j, retirementincome(j)
            end do

        ! Populate income process array
        do j=1, zgridsize
            do t=1, Tretire
                income(j, t)=exp(znodes(j)+ageearnings(t))
            end do
            income(j, Tretire+1:Tdie)=retirementincome(j)
        end do
        ! Lump sum retirement account transfer upon retirement
        income(:, Tretire+1) = income(:,Tretire+1) + &
            ret_wealth*income(:,Tretire)

        ! Progressive income tax (not double charged on retirees)
        ! TODO: in the future a mortgage interest deduction can also be applied
        ! here (so also varies over a)
        do j=1, zgridsize
            do t=1, Tretire
                inctax(j, t) = tax1*(income(j, t) - 0.0)**(1-tax2)
            end do
            inctax(j, Tretire+1:Tdie) = 0.0
        end do

        do t=1, Tdie
            do j=1, zgridsize
                if (income(j, t) - inctax(j, t) <= 0) write(*,*) "errrr"
            end do
        end do

        CLOSE(31)
        CLOSE(32)
        CLOSE(37)

    end subroutine ! %>

    subroutine tauchen
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! !%<
        !
        !    Tauchen algorithm for discretizing an AR(1),
        !    originally found in BGLV code.
        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        IMPLICIT NONE
        INTEGER :: i, j, k
        REAL(8), parameter :: sigma_eta_init=(sigma_z**2.0/(1-rho_z**2.0))**.5
        REAL(8), parameter :: mu = 0.0
        REAL(8), parameter :: zmin=-2.5*sigma_eta_init
        REAL(8), parameter :: zmax=2.5*sigma_eta_init
        REAL(8) :: w
   
        ! Shocks grid - equal intervals a la Tauchen discretization
        write(0,*) NEW_LINE('A') // "// Income shocks grid:"
        do i=1, zgridsize
            znodes(i)= mu+zmin+((zmax-zmin)/(zgridsize-1))*(1.0*i-1.0)
            write(0,'(F16.10)') znodes(i)
        end do
   
        ! create transition matrix for log idiosyncratic labor shock using Tauchen 86
        w=znodes(2)-znodes(1)
        do j=1, zgridsize
            Probz(j, 1)=cdfnormal((znodes(1)-rho_z*znodes(j)-(1.0-rho_z)*mu+w/2)/(sigma_z))
            Probz(j, zgridsize)=1-cdfnormal((znodes(zgridsize)-rho_z*znodes(j)-(1.0-rho_z)*mu-w/2)/(sigma_z))
            do k=2, zgridsize-1
                Probz(j, k)=cdfnormal((znodes(k)-rho_z*znodes(j)-(1.0-rho_z)*mu+w/2)/(sigma_z))&
                           &-cdfnormal((znodes(k)-rho_z*znodes(j)-(1.0-rho_z)*mu-w/2)/(sigma_z))
            end do
        end do
        call prob_trunc

    end subroutine tauchen ! %>

    subroutine rouwenhorst
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ! %<
        !
        !    Rouwenhorst (1996) algorithm for discretizing an AR(1),
        !    as coded by Kopecky and Suen and adapted to Fortran by TC.
        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        IMPLICIT NONE 
        INTEGER :: i
        REAL(8), PARAMETER :: q = (rho_z+1.0)/2.0
        REAL(8), PARAMETER :: mu = 0.0 !-0.06
        REAL(8), PARAMETER :: nu = sigma_z*(((zgridsize-1)/(1-rho_z**2.0))**0.5)
        REAL(8), DIMENSION(:,:), ALLOCATABLE :: Ptemp
        REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: P_iter

        write(0,*) NEW_LINE('A') // "// Income shocks grid:"
        do i=1,zgridsize
            znodes(i) = mu - nu + (2.0*nu/(zgridsize-1))*(1.0*i-1.0)
            write(0,'(F16.10)') znodes(i)
        end do

        ALLOCATE(Ptemp(2,2))
        Ptemp = reshape((/ q, 1.0-q, 1.0-q, q /), (/2,2/))
        do i=3,zgridsize
            ALLOCATE(P_iter(i,i,4))
            P_iter=0
            P_iter(1:i-1,1:i-1,1) = Ptemp
            P_iter(1:i-1,2:i,2) = Ptemp
            P_iter(2:i,1:i-1,3) = Ptemp
            P_iter(2:i,2:i,4) = Ptemp
            DEALLOCATE(Ptemp)
            ALLOCATE(Ptemp(i,i))
            Ptemp = q*P_iter(:,:,1) + (1.0-q)*P_iter(:,:,2) + &
                    (1.0-q)*P_iter(:,:,3) + q*P_iter(:,:,4)
            DEALLOCATE(P_iter)
        end do
        Probz = Ptemp
        call prob_trunc

    end subroutine rouwenhorst ! %>

    subroutine prob_trunc
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ! %<
        !
        !    BGLV code for truncating the transition matrix
        !    to ignore AR(1) state transitions with low probability
        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        IMPLICIT NONE
        REAL(8) :: foundmin, foundmax
        INTEGER :: j, k
    !Truncating the transistion probabilities of transitory shock. Want to
    !minimize number of gridpoints have to compute conditional expectation
        minexpectationz=1
        maxexpectationz=zgridsize
        do j=1, zgridsize
            foundmin=0.0
            foundmax=0.0
            do k=1, zgridsize
                if (Probz(j, k)>.01) then
                    foundmin=1.0
                elseif (foundmin==0.0) then
                    minexpectationz(j)=minexpectationz(j)+1
                end if
            end do
            do k=0, zgridsize-1
                if (Probz(j, zgridsize-k)>.01) then
                    foundmax=1.0
                elseif (foundmax==0.0) then
                    maxexpectationz(j)=maxexpectationz(j)-1
                end if
            end do
        end do

        do j=1, zgridsize
            do k=1, zgridsize
                if (k<minexpectationz(j)) then
                    !write(*,*) Probz(j, k)
                    Probz(j, k)=0.0
                end if
                if (k>maxexpectationz(j)) then
                    !write(*,*) Probz(j, k)
                    Probz(j, k)=0.0
                end if
            end do
        end do

    end subroutine ! %>

    subroutine bequests(price, init_period)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ! %<
        !
        !    Calculate bequest utilities for all state spaces
        !    in an agent's last year of life.
        !
        !    Modified: 08/09/2019
        !
        !    PARAMETERS
        !
        !    price: Real: records the price per unit of housing services.
        !
        !    init_period: index on the overall policy array where the
        !    calculation is to be stored. E.g. to solve the model in GE,
        !    the initial period needs to be "2" (1 is for stationary eq.
        !    values).
        !
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        IMPLICIT NONE
        REAL(8), DIMENSION(:), INTENT(IN) :: price
        INTEGER, INTENT(IN) :: init_period
        INTEGER :: i, j, k
        REAL(8) :: finwealth, wholder, p_const
        REAL(8), DIMENSION(zgridsize) :: realwealth
        INTEGER :: pricelen

        pricelen=MAX(SIZE(price,1)-1,1)
        ! Consumption share constant for the price index, to the -1 power
        p_const =(elasticity2**elasticity2)*((1.0-elasticity2)**(1.0-elasticity2))

    ! Put in bequest motive in last period

        write(0,*) NEW_LINE('A') // "// Sample bequest utility in final period:"
        do i=1, agridsize
            do j=1, Dgridsize
                do k=1, pricelen
                    ! This is the a term
                    finwealth=anodes(i)-(1-delta+maint*delta-dtau)*(1-theta)*Dnodes(j)*exp(price(k))
                    ! This is the a term + d term
                    if (finwealth<0) then
                        wholder=(1+rborrow)*finwealth+(1-delta+maint*delta-dtau)*Dnodes(j)*exp(price(k))
                    !    EV(i, j,:, k, Tdie+1)=1/((1-elasticity))*((finwealth*(1+rborrow)+(1-delta-dtau)*Dnodes(j)*exp(price(k))*(1-F)))**(1-elasticity)
                    else
                        wholder=(1+r)*finwealth+(1-delta+maint*delta-dtau)*Dnodes(j)*exp(price(k))
                    !    EV(i, j,:, k, Tdie+1)=1/((1-elasticity))*((finwealth*(1+r)+(1-delta-dtau)*Dnodes(j)*exp(price(k))*(1-F)))**(1-elasticity)
                    end if
                    ! Bequest utility is in real terms
                    realwealth(:)=(income(:, Tdie)+balancer(k)+wholder+beq_base)*p_const/&
                        ((rent*exp(price(k)))**(1.0-elasticity2))
                    if (elasticity .ne. 1.0) then
                    EV(:,2,j,i,Tdie+1,init_period+k-1)=1.0/(1.0-elasticity)*psi*(realwealth(:)&
                        )**(1.0-elasticity)
                    else
                    EV(:,2,j,i,Tdie+1,init_period+k-1)=psi*log(realwealth(:))
                    end if
                    if ((j == 5) .and. (k == 1) .and. ((i==20) .or. (i==agridsize))) then
                        write(0,*) i, j, k, EV(:,2,j,i,Tdie+1,init_period+k-1)
                    end if
                end do
            end do
        end do
        EV(:,1,:,:,:,:)=0.0  ! Can't be unemployed in retirement

        write(0,*) NEW_LINE('A') // "// EV sum over", pricelen, "periods:"
        write(0,*) sum(EV(:,:,1:Dgridsize-10,10:agridsize,:,init_period:init_period+pricelen-1))
        !EV=max(-10000000000.0, EV)
        !EV=0


    end subroutine ! %>

    ! Generate initial income shock
    subroutine gen_init_shock(numHH, out_array)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ! %<
        !
        !    Call the Fortran compiler's internal random number
        !    generator to create the sequence of shocks for simulated
        !    individuals in period 1 of the simulation.
        !
        !    Modified: 06/22/2017
        !
        !    PARAMETERS
        !
        !    numHH: Marks the number of agents simulated.
        !
        !    out_array: Output matrix storing the random numbers.
        !    
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        IMPLICIT NONE
        INTEGER, INTENT(IN) :: numHH
        REAL(8), dimension(:), intent(INOUT) :: out_array
        INTEGER :: i

        do i=1, numHH
            out_array(i) = minloc(abs(Probinitcum - (real(i)/numHH)), 1)
            if (real(i)/numHH > Probinitcum(out_array(i))) then
                out_array(i) = out_array(i) + 1
            end if
        end do
        !do i=1, numHH
        !    if (real(i)/numHH<Probinitcum(1)) then
        !        init(1, i)=1
        !    else
        !        do j=1, zgridsize-zcut-1
        !            if (real(i)/numHH > Probinitcum(j) .AND. real(i)/numHH<=Probinitcum(j+1)) then
        !                init(1, i)=j+1
        !                exit
        !            end if
        !        end do
        !    end if
        !end do

    end subroutine ! %>

    subroutine gen_life_shocks(numHH, shocks, out_array)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ! %<
        !
        !    Call the Fortran compiler's internal random number
        !    generator to create the sequence of shocks for simulated
        !    individuals over the lifecycle.
        !
        !    Modified: 10/23/2019
        !
        !    PARAMETERS
        !
        !    numHH: Marks the number of agents simulated.
        !
        !    shocks: Period 1 shocks generated in gen_init_shock.
        !
        !    out_array: Output matrix storing the random numbers.
        !    
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



        IMPLICIT NONE
        INTEGER, INTENT(IN) :: numHH
        REAL(8), dimension(:,:), intent(IN) :: shocks
        REAL(8), dimension(:,:), intent(INOUT) :: out_array
        INTEGER :: shock_time, t, i

        call gen_init_shock(numHH, out_array(:, 1))
        shock_time = SIZE(out_array, 2)

        do i=1, numHH
            do t=2,shock_time
                out_array(i,t) = minloc(abs(Probzcum(out_array(i,t-1), :) -&
                                              shocks(i,t-1)), 1)
                if (shocks(i,t-1) > Probzcum(out_array(i,t-1),out_array(i,t))) then
                    out_array(i,t) = out_array(i,t) + 1
                end if
            end do
        end do

    end subroutine ! %>

    ! These initial values are from the create_initial_conditions Matlab files
    ! Each if block is a quantile of the income distribution, so quartiles for now.
    ! The shock variable determines if agent is house owner. Else agent rents.
    subroutine cohort_calibrate(numHH, price, init, fthb, rentflag, c_shock) !%<
        IMPLICIT NONE
        REAL(8), intent(IN) :: price
        REAL(8), dimension(:), intent(INOUT) :: fthb, rentflag
        REAL(8), dimension(:,:), intent(INOUT) :: init, c_shock
        REAL(8), dimension(:,:), ALLOCATABLE :: distMat
        REAL(8), dimension(:), ALLOCATABLE :: q_ecdf, h_ecdf, a_coefs, houseprob
        REAL(8), DIMENSION(2) :: init_unemp
        ! Starting values, statuses to be imputed
        INTEGER, INTENT(IN) :: numHH
        REAL(8) :: medianInc, bound_init
        INTEGER :: i, j

        ! Generate initial distribution of income states (unless overriden)
        call gen_init_shock(numHH, init(1, :))
        if (custom_incs) then
            init(1, :) = 8
        end if

        ! Steady-state unemployment rates (TBD)
        if (no_unemp) then
            init_unemp = 1.0
        else
            init_unemp = ergodic(RESHAPE((/ empprob(1), 1.0-unempprob,&
                1.0-empprob(1), unempprob /), (/ 2, 2 /)), 2)
        end if

        ! Calibration_mat is expressed in units of 19-23 median income, but
        ! everything else in model is expressed in lifetime income. This
        ! is the conversion.
        ! See SCF_calibration.do for generation of these quantiles
        medianInc=income((zgridsize)/2, 1)
        ALLOCATE(distMat(6,6))
        ALLOCATE(h_ecdf(5), q_ecdf(79), a_coefs(3),&
                 houseprob(2))

        distMat = 0.0
        OPEN(34, FILE='input_data/initial_grid_calibration.txt', STATUS='OLD')
        do i=1, 6
            READ(34,*) distMat(:,i)
        end do

        ! An elaborate procedure: we have 6 different bins (below median/below 3rd q/above
        ! income dist X homeowner/renter). For each of these a gamma dist of assets
        ! and a uniform dist of housing fit in SCF. The homeownership probs
        ! are taken over just income bins. These are made for ages 22-25.
        do i=1, numHH
            ! %<
            if (c_shock(4, i) < init_unemp(1)) then
                init(2, i) = 2.0
            else
                init(2, i) = 1.0
            end if

            if ((income(init(1, i), 1) <= distMat(1,2))) then
                if (c_shock(1, i) > distMat(2,2)) then
                    rentflag(i) = 1
                    init(3, i) = 0.0
                    call gamma_inc_inv(distMat(3, 1), init(4,i), 0.0, &
                                       c_shock(2,i), 1.0-c_shock(2,i), bound_init)
                    init(4, i) = init(4, i)*distMat(4, 1)
                    !init(4, i) = 0
                else
                    fthb(i) = 0
                    rentflag(i) = 0
                    ! Draw from the uniform distribution
                    init(3, i) = distMat(5, 2) + &
                        c_shock(3, i)*(distMat(6, 2) - distMat(5, 2))
                    init(3, i) = max(init(3, i), Dmin)
                    ! Get the quantile from the gamma CDF (the lower incomplete gamma ratio f'n)
                    ! for a probability. Then scale by scale param.
                    call gamma_inc_inv(distMat(3, 2), init(4,i), 0.0, &
                                       c_shock(2,i), 1.0-c_shock(2,i), bound_init)
                    init(4, i) = init(4, i)*distMat(4, 2)
                    !init(4, i) = 0
                end if
            else if (income(init(1, i), 1) <= distMat(1,4)) then
                if (c_shock(1, i) > distMat(2,4)) then
                rentflag(i) = 1
                init(3, i) = 0.0
                    call gamma_inc_inv(distMat(3, 3), init(4,i), 0.0, &
                                       c_shock(2,i), 1.0-c_shock(2,i), bound_init)
                    init(4, i) = init(4, i)*distMat(4, 3)
                    !init(4,i) = 0
                else
                    fthb(i) = 0
                    rentflag(i) = 0
                    init(3, i) = distMat(5, 4) + &
                        c_shock(3, i)*(distMat(6, 4) - distMat(5, 4))
                    init(3, i) = max(init(3, i), Dmin)
                    call gamma_inc_inv(distMat(3, 4), init(4,i), 0.0, &
                                       c_shock(2,i), 1.0-c_shock(2,i), bound_init)
                    init(4, i) = init(4, i)*distMat(4, 4)
                    !init(4,i) = 0
                end if
            else
                if (c_shock(1, i) > distMat(2,6)) then
                    rentflag(i) = 1
                    init(3, i) = 0.0
                    call gamma_inc_inv(distMat(3, 5), init(4,i), 0.0, &
                                       c_shock(2,i), 1.0-c_shock(2,i), bound_init)
                    init(4, i) = init(4, i)*distMat(4, 5)
                    !init(4,i) = 0
                else
                    fthb(i) = 0
                    rentflag(i) = 0
                    init(3, i) = distMat(5, 6) + &
                       c_shock(3,i)*(distMat(6, 6) - distMat(5, 6))
                    init(3, i) = max(init(3, i), Dmin)
                    call gamma_inc_inv(distMat(3, 6), init(4,i), 0.0, &
                                       c_shock(2,i), 1.0-c_shock(2,i), bound_init)
                    init(4, i) = init(4, i)*distMat(4, 6)
                    !init(4,i) = 0
                end if
            end if
            if (init(4, i) <= 0) write(*,*) 'error_a'
            ! %>
        end do

        ! init(3,:)=init(3,:) / (sum(Probinit*exp(znodes))*exp(ageearnings(1)))
        ! init(4,:)=init(4,:) / (sum(Probinit*exp(znodes))*exp(ageearnings(1)))
        if (custom_start) then
            init(3,:)=0.0
            init(4,:)=0.0!0.25*(rand()*0.5)
            where (init(3,:) > 0.0) fthb = 0
            where (init(3,:) == 0.0) rentflag = 1
        end if


        CLOSE(34)
        DEALLOCATE(q_ecdf, h_ecdf, a_coefs, distMat)

    end subroutine !%>

end module lifecycle_calibrate

