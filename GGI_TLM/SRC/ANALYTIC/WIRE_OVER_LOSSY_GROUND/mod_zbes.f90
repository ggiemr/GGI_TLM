! File downloaded from: dl.acm.org/doi/10.1145/1916461.1916471

! File mod_zbes.f90

! Author: Masao Kodama
! Address: 21-20, Gakuen 1 chome, Mastue-shi, Shimane-ken, 690-0825 Japan
! Email: mkodama@mable.ne.jp

! MODULE mod_zbes in which the present algorithm is described includes module
! subprograms computing the cylindrical functions:
! * the Bessel function of the first kind Jnu(z),
! * the Bessel function of the second kind Ynu(z), that is, the Neumann
!   function Nnu(z),
! * the Hankel function of the first kind Hnu(1)(z),
! * the Hankel function of the second kind Hnu(2)(z),
! of complex order nu and complex argument z. -pi<arg(z)<=pi.
! MODULE mod_zbes is written in Fortran 90.

! The directions for the algorithm

! In all the files in archive XXX.zip, the following notations are applied:
! zbessel1(znu,zz)=the Bessel function of the first kind Jnu(z).
! zbessel2(znu,zz)=the Bessel function of the second kind Ynu(z).
! zhankel1(znu,zz)=the Hankel function of the first kind Hnu(1)(z).
! zhankel2(znu,zz)=the Hankel function of the second kind Hnu(2)(z).
! znu= complex order nu.
! zz= complex argument z.
! zans=the complex value of the invoked cylindrical function, that is, one of
!      zbessel1, zbessel2, zhankel1 and zhankel2.
! kp=the kind type parameter for the real numbers and the complex numbers used
!    in MODULE mod_zbes.
! nregion=the integer which is equal to the suffix n of Region Rn defined in
!         Section 2 of Ref. (1).
! zunit=the imaginary unit=(0,1)
! pi=the circle ratio=3.1415926...
! epsilon0=EPSILON(1._kp),   epsilon1=MAX(1E-20_kp,epsilon0)
! ai_arg_m=1E3
! re_znu=REAL(znu),   ai_znu=AIMAG(znu)
! re_zz=REAL(zz),     ai_zz=AIMAG(zz)


! The module subprograms belonging to MODULE mod_zbes are classified into the
! following two groups.

! (A) The subprograms of the public attribute

! SUBROUTINE bessel1(znu,zz,zans,info)
! This SUBROUTINE calculates zans. zans=zbessel1(znu,zz).
! This invokes subprograms bes1_series, bes2_series, cylin_inte, num_region and
! abs1.

! SUBROUTINE bessel2(znu,zz,zans,info)
! This SUBROUTINE calculates zans. zans=zbessel2(znu,zz).
! This invokes subprograms bes1_series, bes2_series, cylin_inte, num_region and
! abs1.

! SUBROUTINE hankel1(znu,zz,zans,info)
! This SUBROUTINE calculates zans. zans=zhankel1(znu,zz).
! This invokes subprograms bes1_series, bes2_series, cylin_inte, num_region and
! abs1.

! SUBROUTINE hankel2(znu,zz,zans,info)
! This SUBROUTINE calculates zans. zans=zhankel2(znu,zz).
! This invokes subprograms bes1_series, bes2_series, cylin_inte, num_region and
! abs1.

! znu:  input,  complex, the order of the invoked function.
! zz:   input,  complex, the argument of the invoked function.
! zans: output, complex, the value of the invoked cylindrical function.
! info: output, integer, the output condition
!   info=0:  normal output; relative error of zans is less than 
!            MAX(epsilon1,1E3*epsilon0)
!   info=10: zans is not reliable because of one of the following reasons.
!            (1) There is a possibility of an overflow.
!            (2) There is a possibility of an underflow.
!            (3) The precision of zans is not sufficient.
!            (4) The output zans is indefinite theoretically.
!                Function zbessel1(zunit,0) is indefinite for example.
!   info=20: out of range. The cases of ABS(re_zz)>ai_arg_m for example.

! (B) The subprograms of the private attribute

! SUBROUTINE bes1_series(znu,zz,zsum,zlogbes,info)
! This SUBROUTINE calculates zsum and zlogbes by the use of the series expansion
! [Ref. (9), 9.1.10] when nregion is 1 or 2. zbessel1(znu,zz)=zsum*EXP(zlogbes).
! This invokes abs1. This is invoked by bessel1, bessel2, hankel1 and hankel2.
! The internal subprograms of bes1_series are zgamma and subgam, which calculate
! the gamma function of a complex argument.

! SUBROUTINE bes2_series(znu,zz,zbes2,info)
! This SUBROUTINE calculates zbes2 by the use of the series expansion when
! nregion=2 and re_znu>=0 [Ref. (3)]. zbes2=zbessel2(znu,zz)
! This is invoked by bessel1, bessel2, hankel1 and hankel2.
! The internal subprograms of bes2_series are bes2_srs_init, def_bessel1, multi,
! multi_z, divis and divis_z.

! SUBROUTINE cylin_inte(znu,zz,zbes1_t,zhan1_t,zhan2_t,nfv,info)
! This SUBROUTINE calculates zbes1_t, zhan1_t zhan2_t and nfv by the use of the
! method of numerical integration when nregion=3 [Ref. (1)].
! zbessel1(znu,zz)=zbes1_t*zex**nfv(1)
! zhankel1(znu,zz)=zhan1_t*zex**nfv(2)
! zhankel2(znu,zz)=zhan2_t*zex**nfv(3)
! zex=EXP(pi*znu*zunit)
! nfv is an integer array.
! This is invoked by bessel1, bessel2, hankel1 and hankel2. This invokes abs1.
! The internal subprograms of cylin_inte are integration, mposition, abs2,
! root1, zfun_stand, dfzero2, re_zfun2, newton, cal_phi, zdgaus8d and zfun.
! Subroutines root1 and dfzero2 are made by modifying SUBROUTINE DFZERO, which
! is downloaded from the site http://gams.nist.gov/. SUBROUTINE DFZERO was made
! by Shampine, L. F. (SNLA) and Watts, H. A. (SNLA). SUBROUTINE DFZERO uses a 
! combination of bisection and the secant rule.
! Subroutine zdgaus8d is written by modifying subroutine DGAUS8, which was 
! produced by Jones, R. E. (SNLA). SUBROUTINE DGAUS8 can be downloaded from site
! http://gams.nist.gov/. SUBROUTINE DGAUS8 uses an adaptive 8-point Legendre
! -Gauss algorithm.

! INTEGER FUNCTION num_region(znu,zz)
! This SUBROUTINE computes nregion. nregion=num_region(znu,zz)
! This is invoked by bessel1, bessel2, hankel1 and hankel2.

! REAL FUNCTION abs1(za)
! FUNCTION abs1 calculates a rough absolute value of the complex argument za.
! This is invoked by bessel1, bessel2, hankel1, hankel2, num_region,
! bes1_series, bes2_series, bes2_srs_init, def_bessel1, cylin_inte, cal_phi and
! zdgaus8d.


! Comments

! The subroutines that the user's program can invoke directly are only the four
! SUBROUTINEs bessel1, bessel2, hankel1 and hankel2 of the public attribute.
! If one of SUBROUTINEs bessel1, bessel2, hankel1 and hankel2 is invoked by the
! user's program, this subroutine invokes FUNCTION num_region first and
! computes nregion, and invokes some of the SUBROUTINEs bes1_series,
! bes2_series and cylin_inte of the private attribute according to nregion, and
! calculates the invoked cylindrical function.

! In this algorithm, an arbitrary kp is available if the used processor accepts
! the kp. The user can select kp in the beginning of MODULE mod_zbes of this 
! file. The default value of kp is KIND(1D0). The data type of the actual 
! arguments znu, zz and zans must accord with the data type determined by kp.

! The parameter kp has the public attribute. All the variables and the 
! parameters except kp have the private attribute. 

! The variable epsilon1 is used as a standard of the relative errors of the
! cylindrical functions computed by the present algorithm, so that the relative
! relative errors of zans may be about 1E-20 even if the rounding errors in
! mathematical operations are negligible. 

! If an overflow or an underflow occurs when kp=KIND(1D0), it is very effective
! in avoiding the overflow or the underflow to use extended precision, in which
! kp=SELECTED_REAL_KIND(20,550).

! A complex variable za itself has a relative error of about epsilon0.
! Hence, the relative error of EXP(za) is about ABS(za)*epsilon0.
! See [Ref. (7), Section 3.1]. The error of SIN(za) or COS(za) is also about
! ABS(za)*epsilon0 except when SIN(za) or COS(za) is nearly equal to 0.
! Hence, in order to protect unlimited increase in errors of the intrinsic
! functions EXP, SIN and COS, we must put an upper limit to ABS(za),
! because the errors of the intrinsic functions have a bad influence on the
! accuracy of zans. In the present algorithm, the upper limit to ABS(za) is 
! ai_arg_m, so that the error of a cylindrical function calculated by the
! present algorithm is less than about 1E3*epsilon0 except the neighborhoods of
! the zeros of the function. When ABS(za) is larger than ai_arg_m, info is set
! at 20, and the calculation stops.

! When an overflow or an underflow of a variable occurs on the middle of
! numerical calculation of a cylindrical function, info is set at 10 and the
! calculation stops. Hence, the case of info=10 does not always designate the
! overflow or the underflow of the value of the invoked cylindrical function
! itself.

! Underflows may have bad influence on the precision of zans [Ref. (1), Section
! 4.2]. We must pay attention to the precision of zans even if ABS(zans) is
! much larger than TINY(1._kp). However if a Fortran system uses the 
! denormalized numbers [Ref. (5), 8.7.1], the influence of underflows greatly 
! decreases. 

! The algorithm is designed on the assumption that RADIX(1._kp)=2 [Ref. (1),
! Section 4.2]. Hence, it is desirable that RADIX(1._kp)=2.


! References

! (1) Masao Kodama, "Algorithm XXX: A module for calculating cylindrical
! functions of complex order and complex argument," ACM Transactions on
! Mathematical Software.

! (2) Hirondo Kuki, "Algorithm 421   Complex gamma function with error control
! [S14]," Communications of the ACM, vol. 15, no. 4, pp. 271-272, April 1972.

! (3) M. Kodama, M. Yamasato and S. Yamashiro, "Numerical calculation of the
! Neumann function Nnu(x) of complex order nu," IEICE Trans. Fundamentals, vol.
! E78-A, no. 6, pp. 727-736, June 1995.

! (4) Masao Kodama and Kengo Taira, "Polynomials approximating complex
! functions," IEICE Trans. Fundamentals, vol. E80-A, no. 4, pp. 778-781, April
! 1997.

! (5) Michael Metcalf and John Reid, "Fortran 90/95 Explained (Second Edition),"
! Oxford University Press, 1999.

! (6) D. E. Amos, "Algorithm 644   A portable package for Bessel functions of a
! complex argument and nonnegative order," ACM Transactions on Mathematical
! Software, vol. 12, no. 3, pp. 265-273, Sept. 1986.

! (7) Masao Kodama, "Algorithm 877: A subroutine package for cylindrical
! functions of complex order and nonnegative argument", ACM Transactions on
! Mathematical Software, vol. 34, no. 4, pp. 22:1-22:21, July 2008.

! (8) Amparo Gil, Javier Segura and Nico M. Temme, "Algorithm 819   AIZ, BIZ:
! two Fortran 77 routines for the computation of complex Airy functions,"
! ACM Transactions on Mathematical Software, vol. 28, no. 3, pp. 325-336,
! September 2002.

! (9) Milton Abramowitz and Irene A. Stegun, "Handbook of mathematical
! functions with formulas, graphs, and mathematical tables," Dover
! Publications, Inc., New York, 1970.

! (10) Thompson, I. J. and Barnett, A. R., "COULCC: A continued-fraction
! algorithm for Coulomb functions of complex order with complex arguments,"
! Computer Physics Communications, vol. 36, Issue 4, pp. 363-372.

    MODULE mod_zbes
! MODULE mod_besz declares the following named constants: kp, huge1, tiny1,
! ihuge, epsilon0, pi, alog2.
! MODULE mod_besz declares the following variables: epsilon1, alog_huge,
! alog_tiny, ai_znu_m, ai_arg_m, alog_eps1, znua, zza, znupi, znupi_i, ai_znu
! re_znu, re_zz, ai_zz, ai_znupi.
      IMPLICIT NONE
      PRIVATE
! The value of kp can be selected from the following three lines.
!      INTEGER, PARAMETER :: kp = kind(1.) ! For single precision
      INTEGER, PARAMETER :: kp = kind(1D0) ! For double precision and the default
! INTEGER,PARAMETER:: kp=SELECTED_REAL_KIND(20,550)   ! For extended precision
      REAL (kp), PARAMETER :: huge1 = huge(1._kp)/3, tiny1 = tiny(1._kp)
      INTEGER, PARAMETER :: ihuge = huge(1) - 2 ! =2.1E9
      REAL (kp), PARAMETER :: epsilon0 = epsilon(1._kp) ! =2.2E-16
      REAL (kp), PARAMETER :: ai_arg_m = 1E3 ! epsilon0*ai_arg_m=max relative error
      REAL (kp), PARAMETER :: pi = 3.14159265358979323846264_kp ! the circle ratio
! The values of variables declared in the following type declaration
! statements are defined in the function num_region.
      REAL (kp), SAVE :: epsilon1, alog_huge = -1, alog_tiny, ai_znu_m, &
        alog_eps1
      REAL, SAVE :: alog2 ! =LOG(2.)
      COMPLEX (kp) :: znua, zza, znupi, znupi_i
      REAL (kp) :: ai_znu, re_znu, re_zz, ai_zz, ai_znupi
      PUBLIC :: kp, bessel1, bessel2, hankel1, hankel2

    CONTAINS

      SUBROUTINE bessel1(znu,zz,zans,info)
        COMPLEX (kp), INTENT (IN) :: znu, zz
        COMPLEX (kp), INTENT (OUT) :: zans
        INTEGER, INTENT (OUT) :: info
! This SUBROUTINE computes zans. zans=zbessel1(znu,zz).
! The virtual arguments znu, zz, zans, info, and the outline of calculation
! in SUBROUTINE bessel1 are explained in the beginning of this file.
! The named constants kp, huge1, ai_arg_m, alog2 and variables alog_huge,
! ai_znu_m, zex, znupi, znupi_i, re_znu, ai_znu, re_zz, ai_zz, ai_znupi are
! declared in MODULE mod_zbes.
! This is invoked only by the user's program.
! This subroutine invokes subprograms bes1_series, bes2_series, cylin_inte,
! num_region and abs1.
        COMPLEX (kp) :: zlogbes, zsum, zbes1_t, zhan1_t, zhan2_t, zbes1a, &
          zbes2a, zsum1, zlogbes1, z1, zex
        REAL :: a1
        INTEGER :: nregion, nfv(3), info1

        zans = 0
! Determination of nregion
        nregion = num_region(znu,zz)
        IF (nregion==0) THEN
          info = 20
          RETURN
        END IF
        SELECT CASE (nregion)
        CASE (1) ! By the series expansion. ABS(zz)<1, ABS(zeps)>0.3
          CALL bes1_series(znu,zz,zsum,zlogbes,info)
          IF (info>1) RETURN
          IF (real(zlogbes)>alog_huge) THEN
            zans = huge1
            info = 10
            RETURN
          END IF
          zans = zsum*exp(zlogbes)
        CASE (2) ! By the series expansion. ABS(zz)<1, ABS(zeps)<=0.3
          IF (re_znu>=0) THEN !  ABS(zz)<1, ABS(zeps)<=0.3, re_znu>=0.
            CALL bes1_series(znu,zz,zsum,zlogbes,info)
            IF (info>1) RETURN
            IF (real(zlogbes)>alog_huge) THEN
              zans = huge1
              info = 10
              RETURN
            END IF
            zans = zsum*exp(zlogbes)
            RETURN
          END IF
! ABS(zz)<1, ABS(zeps)<=0.3, re_znu<0. See Eq. (8) of Ref. (1).
          CALL bes1_series(znua,zz,zsum1,zlogbes1,info)
          IF (info>1) RETURN
          IF (real(zlogbes1)>alog_huge) THEN
            zans = huge1
            info = 10
            RETURN
          END IF
          zbes1a = zsum1*exp(zlogbes1)
          CALL bes2_series(znua,zz,zbes2a,info1)
          IF (info1>11) RETURN
          IF (info1==10 .AND. abs1(znu-nint(re_znu))<epsilon1*10) THEN
            zans = zbes1a*cos(znupi)
            info = 0
            RETURN
          END IF
          IF (info1>1) THEN
            info = info1
            RETURN
          END IF
          zans = zbes1a*cos(znupi) + zbes2a*sin(znupi)
        CASE (3) ! By numerical integration. ABS(zz)>=1
          IF (ai_znu<-ai_znu_m) THEN
            zans = huge1
            info = 10
            RETURN
          END IF
          zex = exp(znupi_i)
          IF (re_znu>=0) THEN
            IF (re_zz>=0) THEN ! re_znu>=0, re_zz>=0.
              CALL cylin_inte(znu,zz,zbes1_t,zhan1_t,zhan2_t,nfv,info)
              IF (info>1) RETURN
              zans = zbes1_t*zex**nfv(1)
            ELSE ! re_znu>=0, re_zz<0
              CALL cylin_inte(znu,zza,zbes1_t,zhan1_t,zhan2_t,nfv,info)
              IF (info>1) RETURN
              IF (ai_zz>=0) THEN
! re_znu>=0, re_zz<0, ai_zz>=0. See Eq. (38a) of Ref. (1).
                IF (exponent(abs1(zbes1_t))*alog2-ai_znupi*(nfv( &
                    1)+1)>alog_huge) THEN
                  zans = huge1
                  info = 10
                  RETURN
                END IF
                zans = zbes1_t*zex**(nfv(1)+1)
              ELSE ! re_znu>=0, re_zz<0, ai_zz<0. See Eq. (38b) of Ref. (1).
                IF (exponent(abs1(zbes1_t))*alog2-ai_znupi*(nfv( &
                    1)-1)>alog_huge) THEN
                  zans = huge1
                  info = 10
                  RETURN
                END IF
                IF (ai_znu>ai_znu_m) THEN
                  zans = huge1
                  info = 10
                  RETURN
                END IF
                zans = zbes1_t*zex**(nfv(1)-1)
              END IF
            END IF
          ELSE ! re_znu<0
            IF (re_zz>=0) THEN ! re_znu<0, re_zz>=0
              CALL cylin_inte(znua,zz,zbes1_t,zhan1_t,zhan2_t,nfv,info)
              IF (info>1) RETURN
              nfv = -nfv
              IF (ai_znu>=0) THEN
! re_znu<0, ai_znu>=0, re_zz>=0. See Eq. (38c) of Ref. (1).
                IF (exponent(abs1(zhan1_t)/2)*alog2+abs(ai_znupi)-ai_znupi*nfv &
                    (2)>alog_huge) THEN
                  zans = huge1
                  info = 10
                  RETURN
                END IF
                zans = zbes1_t*zex**(1+nfv(1)) - zex**nfv(2)*(zex-1/zex)* &
                  zhan1_t/2
              ELSE ! re_znu<0, ai_znu<0, re_zz>=0. See Eq. (38d) of Ref. (1).
                IF (exponent(abs1(zhan2_t)/2)*alog2+abs(ai_znupi)-ai_znupi*nfv &
                    (3)>alog_huge) THEN
                  zans = huge1
                  info = 10
                  RETURN
                END IF
                zans = zbes1_t*zex**(nfv(1)-1) + zex**nfv(3)*(zex-1/zex)* &
                  zhan2_t/2
              END IF
            ELSE ! re_znu<0, re_zz<0
              CALL cylin_inte(znua,zza,zbes1_t,zhan1_t,zhan2_t,nfv,info)
              IF (info>1) RETURN
              nfv = -nfv
              IF (ai_zz>=0) THEN
! re_znu<0, re_zz<0, ai_zz>=0. See Eq. (38e) of Ref. (1).
                a1 = -ai_znupi*(nfv(3)+2)
                IF (exponent(abs1(zhan2_t)/2)*alog2+a1>alog_huge) THEN
                  zans = huge1
                  info = 10
                  RETURN
                END IF
                IF (a1>alog_huge) THEN
                  zans = huge1
                  info = 10
                  RETURN
                END IF
                zans = (zhan1_t*zex**nfv(2)+zhan2_t*zex**(nfv(3)+2))/2
              ELSE ! re_znu<0, re_zz<0, ai_zz<0. See Eq. (38f) of Ref. (1).
                a1 = ai_znupi*(2-nfv(2))
                IF (exponent(abs1(zhan1_t))*alog2+a1>alog_huge) THEN
                  zans = huge1
                  info = 10
                  RETURN
                END IF
                IF (a1>alog_huge) THEN
                  zans = huge1
                  info = 10
                  RETURN
                END IF
                IF (nfv(2)-2<0) THEN
                  z1 = (1/zex)**(2-nfv(2))
                ELSE
                  z1 = zex**(nfv(2)-2)
                END IF
                zans = (zhan1_t*z1+zhan2_t*zex**nfv(3))/2
              END IF
            END IF
          END IF
        END SELECT
      END SUBROUTINE bessel1

      SUBROUTINE bessel2(znu,zz,zans,info)
        COMPLEX (kp), INTENT (IN) :: znu, zz
        COMPLEX (kp), INTENT (OUT) :: zans
        INTEGER, INTENT (OUT) :: info
! This SUBROUTINE calculates zans. zans=zbessel2(znu,zz).
! The virtual arguments znu, zz, zans, info, and the outline of calculation
! in SUBROUTINE bessel2 are explained in the beginning of this file.
! The named constants kp, huge1, ai_arg_m, alog2 and variables alog_huge,
! ai_znu_m, zex, znupi, znupi_i, re_znu, ai_znu, re_zz, ai_zz, ai_znupi are
! declared in MODULE mod_zbes.
! This is invoked only by the user's program.
! This subroutine invokes subprograms bes1_series, bes2_series, cylin_inte,
! num_region and abs1.
        COMPLEX (kp) :: zarg1, zarg2, zbes1a, zbes2a, zlogbes1, zlogbes2, &
          zpart1, zpart2, zsum1, zsum2, zbes1_t, zhan1_t, zhan2_t, za, zb, zc, &
          zex
        REAL :: a1
        INTEGER :: nregion, nfv(3)

        zans = 0
! Determination of nregion
        nregion = num_region(znu,zz)
        IF (nregion==0) THEN
          info = 20
          RETURN
        END IF
        SELECT CASE (nregion)
        CASE (1) ! By the series expansion.  ABS(zz)<1, ABS(zeps)>0.3.
          CALL bes1_series(znu,zz,zsum1,zlogbes1,info)
          IF (info>1) RETURN
          IF (real(zlogbes1)>alog_huge) THEN
            zans = huge1
            info = 10
            RETURN
          END IF
          CALL bes1_series(znua,zz,zsum2,zlogbes2,info)
          IF (info>1) RETURN
          IF (ai_znu>=0) THEN ! ai_znu>=0. See Eq. (5a) of Ref. (1).
            zarg1 = zlogbes1 + 2*znupi_i
            zarg2 = zlogbes2 + znupi_i
            IF (real(zarg1)>alog_huge .OR. real(zarg2)>alog_huge) THEN
              zans = huge1
              info = 10
              RETURN
            END IF
            IF ((abs(aimag(zarg1))>ai_arg_m) .OR. (abs( &
                aimag(zlogbes1))>ai_arg_m)) THEN
              info = 20
              RETURN
            END IF
            zpart1 = zsum1*(exp(zarg1)+exp(zlogbes1))
            zpart2 = -2*zsum2*exp(zarg2)
            zex = exp(znupi_i)
            za = (zpart1+zpart2)/(zex**2-1)
            zans = cmplx(-aimag(za),real(za),kp)
          ELSE ! ai_znu<0. See Eq. (5b) of Ref. (1).
            zarg1 = zlogbes1 - 2*znupi_i
            zarg2 = zlogbes2 - znupi_i
            IF (real(zarg1)>alog_huge .OR. real(zarg2)>alog_huge) THEN
              zans = huge1
              info = 10
              RETURN
            END IF
            zpart1 = zsum1*(exp(zlogbes1)+exp(zarg1))
            zpart2 = -2*zsum2*exp(zarg2)
            za = exp(-znupi_i)
            za = (zpart1+zpart2)/(1-za**2)
            zans = cmplx(-aimag(za),real(za),kp)
          END IF
        CASE (2) ! By the series expansion.  ABS(zz)<1, ABS(zeps)<=0.3
          IF (re_znu>=0) THEN
            CALL bes2_series(znu,zz,zans,info)
            RETURN
          END IF ! re_znu<0. See Eq. (9) of Ref. (1).
          CALL bes1_series(znua,zz,zsum1,zlogbes1,info)
          IF (info>1) RETURN
          IF (real(zlogbes1)>alog_huge) THEN
            zans = huge1
            info = 10
            RETURN
          END IF
          zbes1a = zsum1*exp(zlogbes1)
          CALL bes2_series(znua,zz,zbes2a,info)
          IF (info>1) RETURN
          zans = -zbes1a*sin(znupi) + zbes2a*cos(znupi)
        CASE (3) ! By numerical integration.  ABS(zz)>=1
          IF (ai_znu<-ai_znu_m) THEN
            zans = huge1
            info = 10
            RETURN
          END IF
          zex = exp(znupi_i)
          IF (re_znu>=0) THEN
            IF (re_zz>=0) THEN ! re_znu>=0, re_zz>=0. See Eq. (39a) of Ref. (1).
              CALL cylin_inte(znu,zz,zbes1_t,zhan1_t,zhan2_t,nfv,info)
              IF (info>1) RETURN
              za = (zhan1_t*zex**nfv(2)-zhan2_t*zex**nfv(3))/2
              zans = cmplx(aimag(za),-real(za),kp)
            ELSE ! re_znu>=0, re_zz<0
              CALL cylin_inte(znu,zza,zbes1_t,zhan1_t,zhan2_t,nfv,info)
              IF (info>1) RETURN
              IF (ai_zz>=0) THEN
! re_znu>=0, re_zz<0, ai_zz>=0. See Eq. (39b) of Ref. (1).
                IF (exponent(abs1(zbes1_t))*alog2-ai_znupi*(nfv( &
                    1)+1)>alog_huge) THEN
                  zans = huge1
                  info = 10
                  RETURN
                END IF
                IF (exponent(abs1(zhan2_t))*alog2-ai_znupi*(nfv( &
                    3)-1)>alog_huge) THEN
                  zans = huge1
                  info = 10
                  RETURN
                END IF
                IF (nfv(3)-1>=0) THEN
                  zb = zex**(nfv(3)-1)
                ELSE
                  IF (ai_znu>ai_znu_m) THEN
                    info = 10
                    RETURN
                  END IF
                  zb = exp(-znupi_i)**(1-nfv(3))
                END IF
                za = zbes1_t*zex**(nfv(1)+1) + zhan2_t*zb
                zans = cmplx(-aimag(za),real(za),kp)
              ELSE ! re_znu>=0, re_zz<0, ai_zz<0. See Eq. (39c) of Ref. (1).
                IF (exponent(abs1(zbes1_t))*alog2-ai_znupi*(nfv( &
                    1)-1)>alog_huge) THEN
                  zans = huge1
                  info = 10
                  RETURN
                END IF
                IF (exponent(abs1(zhan1_t))*alog2-ai_znupi*(nfv( &
                    2)+1)>alog_huge) THEN
                  zans = huge1
                  info = 10
                  RETURN
                END IF
                IF (nfv(2)+1>=0) THEN
                  zb = zex**(nfv(2)+1)
                ELSE
                  IF (ai_znu>ai_znu_m) THEN
                    info = 10
                    RETURN
                  END IF
                  zb = exp(-znupi_i)**(-nfv(2)-1)
                END IF
                IF (nfv(1)-1>=0) THEN
                  zc = zex**(nfv(1)-1)
                ELSE
                  IF (ai_znu>ai_znu_m) THEN
                    info = 10
                    RETURN
                  END IF
                  zc = exp(-znupi_i)**(1-nfv(1))
                END IF
                za = zbes1_t*zc + zhan1_t*zb
                zans = cmplx(aimag(za),-real(za),kp)
              END IF
            END IF
          ELSE ! re_znu<0
            IF (re_zz>=0) THEN ! re_znu<0, re_zz>=0. See Eq. (39d) of Ref. (1).
              CALL cylin_inte(znua,zz,zbes1_t,zhan1_t,zhan2_t,nfv,info)
              IF (info>1) RETURN
              nfv = -nfv
              IF (exponent(abs1(zhan1_t))*alog2-ai_znupi*(nfv(2)-1)>alog_huge) &
                  THEN
                zans = huge1
                info = 10
                RETURN
              END IF
              IF (exponent(abs1(zhan2_t))*alog2-ai_znupi*(nfv(3)+1)>alog_huge) &
                  THEN
                zans = huge1
                info = 10
                RETURN
              END IF
              za = zhan1_t*zex**(nfv(2)-1) - zhan2_t*zex**(nfv(3)+1)
              zans = cmplx(aimag(za),-real(za),kp)/2
            ELSE ! re_znu<0, re_zz<0
              CALL cylin_inte(znua,zza,zbes1_t,zhan1_t,zhan2_t,nfv,info)
              IF (info>1) RETURN
              nfv = -nfv
              IF (ai_zz>=0) THEN
! re_znu<0, re_zz<0, ai_zz>=0. See Eq. (39e) of Ref. (1).
                a1 = -ai_znupi*(2+nfv(3))
                IF (exponent(abs1(zhan2_t)*2)*alog2+a1>alog_huge) THEN
                  zans = huge1
                  info = 10
                  RETURN
                END IF
                IF (a1>alog_huge) THEN
                  zans = huge1
                  info = 10
                  RETURN
                END IF
                za = zhan1_t*zex**nfv(2) + (zex**2+2)*zhan2_t*zex**nfv(3)
                zans = cmplx(-aimag(za),real(za),kp)/2
              ELSE ! re_znu<0, re_zz<0, ai_zz<0. See Eq. (39f) of Ref. (1).
                a1 = ai_znupi*(2-nfv(2))
                IF (exponent(abs1(zhan1_t))*alog2+a1>alog_huge) THEN
                  zans = huge1
                  info = 10
                  RETURN
                END IF
                IF (a1>alog_huge) THEN
                  zans = huge1
                  info = 10
                  RETURN
                END IF
                za = (2+(1/zex)**2)*zhan1_t*zex**nfv(2) + zhan2_t*zex**nfv(3)
                zans = cmplx(aimag(za),-real(za),kp)/2
              END IF
            END IF
          END IF
        END SELECT
      END SUBROUTINE bessel2

      SUBROUTINE hankel1(znu,zz,zans,info)
        COMPLEX (kp), INTENT (IN) :: znu, zz
        COMPLEX (kp), INTENT (OUT) :: zans
        INTEGER, INTENT (OUT) :: info
! This SUBROUTINE calculates zans. zans=zhankel1(znu,zz).
! The virtual arguments znu, zz, zans, info, and the outline of calculation
! in SUBROUTINE hankel1 are explained in the beginning of this file.
! The named constants kp, huge1, tiny1, alog2 and variables ai_znu_m, alog_huge,
! zex, znupi_i, re_znu, ai_znu, re_zz, ai_zz, ai_znupi are declared in MODULE
! mod_zbes.
! This is invoked only by the user's program.
! This subroutine invokes subprograms bes1_series, bes2_series, cylin_inte,
! num_region, abs1.
        COMPLEX (kp) :: zarg1, zarg2, zlogbes1, zlogbes2, zbes1, zbes2, &
          zbes1_t, zhan1_t, zhan2_t, zpart1, zpart2, zsum1, zsum2, za, z1, zex
        REAL :: a1
        INTEGER :: nregion, nfv(3)

        zans = 0
! Determination of nregion
        nregion = num_region(znu,zz)
        IF (nregion==0) THEN
          info = 20
          RETURN
        END IF
        SELECT CASE (nregion)
        CASE (1) ! By the series expansion. ABS(zz)<1, ABS(zeps)>0.3.
          CALL bes1_series(znu,zz,zsum1,zlogbes1,info)
          IF (info>1) RETURN
          CALL bes1_series(znua,zz,zsum2,zlogbes2,info)
          IF (info>1) RETURN
          IF (ai_znu>=0) THEN ! ai_znu>=0. See Eq. (6a) of Ref. (1).
            zarg1 = zlogbes1
            zarg2 = zlogbes2 + znupi_i
            IF (real(zarg1)>alog_huge .OR. real(zarg2)>alog_huge) THEN
              zans = huge1
              info = 10
              RETURN
            END IF
            zpart1 = -zsum1*exp(zarg1)
            zpart2 = zsum2*exp(zarg2)
            zex = exp(znupi_i)
            zans = 2*(zpart1+zpart2)/(zex**2-1)
          ELSE !  ai_znu<0. See Eq. (6b) of Ref. (1).
            zarg1 = zlogbes1 - 2*znupi_i
            zarg2 = zlogbes2 - znupi_i
            IF (real(zarg1)>alog_huge .OR. real(zarg2)>alog_huge) THEN
              zans = huge1
              info = 10
              RETURN
            END IF
            zpart1 = -zsum1*exp(zarg1)
            zpart2 = zsum2*exp(zarg2)
            za = exp(-znupi_i)
            zans = 2*(zpart1+zpart2)/(1-za**2)
            IF (abs1(zans)<tiny1*1E4) THEN
              info = 10
              RETURN
            END IF
          END IF
        CASE (2) ! By the series expansion. ABS(zz)<1, ABS(zeps)<=0.3
          IF (re_znu>=0) THEN ! re_znu>=0. See Eq. (10a) of Ref. (1).
            CALL bes1_series(znu,zz,zsum1,zlogbes1,info)
            IF (info>1) RETURN
            IF (real(zlogbes1)>alog_huge) THEN
              zans = huge1
              info = 10
              RETURN
            END IF
            zbes1 = zsum1*exp(zlogbes1)
            CALL bes2_series(znu,zz,zbes2,info)
            IF (info>1) RETURN
            zans = zbes1 + cmplx(-aimag(zbes2),real(zbes2),kp)
            RETURN
          END IF ! re_znu<0. See Eq. (10b) of Ref. (1).
          CALL bes1_series(znua,zz,zsum1,zlogbes1,info)
          IF (info>1) RETURN
          zarg1 = zlogbes1 - znupi_i
          IF (real(zarg1)>alog_huge) THEN
            zans = huge1
            info = 10
            RETURN
          END IF
          CALL bes2_series(znua,zz,zbes2,info)
          IF (info>1) RETURN
          zex = exp(znupi_i)
          za = zbes2/zex
          zans = exp(zarg1)*zsum1 + cmplx(-aimag(za),real(za),kp)
        CASE (3) ! By numerical integration. ABS(zz)>=1
          IF (ai_znu<-ai_znu_m) THEN
            zans = huge1
            info = 10
            RETURN
          END IF
          zex = exp(znupi_i)
          IF (re_znu>=0) THEN
            IF (re_zz>=0) THEN ! re_znu>=0, re_zz>=0
              CALL cylin_inte(znu,zz,zbes1_t,zhan1_t,zhan2_t,nfv,info)
              IF (info>1) RETURN
              zans = zhan1_t*zex**nfv(2)
              RETURN
            END IF
! re_znu>=0, re_zz<0
            IF (ai_znu>ai_znu_m) THEN
              zans = huge1
              info = 10
              RETURN
            END IF
            CALL cylin_inte(znu,zza,zbes1_t,zhan1_t,zhan2_t,nfv,info)
            IF (info>1) RETURN
            IF (ai_zz>=0) THEN
! re_znu>=0, re_zz<0, ai_zz>=0. See Eq. (40a) of Ref. (1).
              a1 = exponent(abs1(zhan2_t))*alog2 - ai_znupi*(nfv(3)-1)
              IF (a1>alog_huge) THEN
                zans = huge1
                info = 10
                RETURN
              END IF
              zans = -zhan2_t*zex**(nfv(3)-1)
            ELSE !  re_znu>=0, re_zz<0, ai_zz<0. See Eq. (40b) of Ref. (1).
              IF (exponent(abs1(zbes1_t)*2)*alog2-ai_znupi*(nfv( &
                  1)-1)>alog_huge) THEN
                zans = huge1
                info = 10
                RETURN
              END IF
              IF (exponent(abs1(zhan1_t))*alog2-ai_znupi*(nfv(2)+1)>alog_huge) &
                  THEN
                zans = huge1
                info = 10
                RETURN
              END IF
              zans = 2*zbes1_t*zex**(nfv(1)-1) + zhan1_t*zex**(nfv(2)+1)
            END IF
          ELSE ! re_znu<0
            IF (re_zz>=0) THEN ! re_znu<0, re_zz>=0. See Eq. (40c) of Ref. (1)
              CALL cylin_inte(znua,zz,zbes1_t,zhan1_t,zhan2_t,nfv,info)
              IF (info>1) RETURN
              nfv = -nfv
              a1 = exponent(abs1(zhan1_t))*alog2 - ai_znupi*(nfv(2)-1)
              IF (a1>alog_huge) THEN
                zans = huge1
                info = 10
                RETURN
              END IF
              zans = zhan1_t*zex**(nfv(2)-1)
            ELSE ! re_znu<0, re_zz<0.
              CALL cylin_inte(znua,zza,zbes1_t,zhan1_t,zhan2_t,nfv,info)
              IF (info>1) RETURN
              nfv = -nfv
              IF (ai_zz>=0) THEN
! re_znu<0, re_zz<0, ai_zz>=0. See Eq. (40d) of Ref. (1).
                zans = -zhan2_t*zex**nfv(3)
              ELSE ! re_znu<0, re_zz<0, ai_zz<0. See Eq. (40e) of Ref. (1).
                a1 = ai_znupi*(2-nfv(2))
                IF (exponent(abs1(zhan1_t))*alog2+a1>alog_huge) THEN
                  zans = huge1
                  info = 10
                  RETURN
                END IF
                IF (a1>alog_huge) THEN
                  zans = huge1
                  info = 10
                  RETURN
                END IF
                IF (nfv(2)-2<0) THEN
                  z1 = (1/zex)**(2-nfv(2))
                ELSE
                  z1 = zex**(nfv(2)-2)
                END IF
                zans = 2*zbes1_t*zex**nfv(1) + zhan1_t*z1
              END IF
            END IF
          END IF
        END SELECT
      END SUBROUTINE hankel1

      SUBROUTINE hankel2(znu,zz,zans,info)
        COMPLEX (kp), INTENT (IN) :: znu, zz
        COMPLEX (kp), INTENT (OUT) :: zans
        INTEGER, INTENT (OUT) :: info
! This SUBROUTINE calculates zans. zans=zhankel2(znu,zz).
! The virtual arguments znu, zz, zans, info, and the outline of calculation
! in SUBROUTINE hankel2 are explained in the beginning of this file.
! The named constants kp, huge1, tiny1, alog2 and variables ai_znu_m, alog_huge,
! zex, znupi_i, re_znu, ai_znu, re_zz, ai_zz, ai_znupi are declared in MODULE
! mod_zbes.
! This is invoked only by the user's program.
! This subroutine invokes subprograms bes1_series, bes2_series, cylin_inte,
! num_region, abs1.
        COMPLEX (kp) :: zarg1, zarg2, zlogbes1, zlogbes2, zbes1, zbes2, &
          zbes1_t, zhan1_t, zhan2_t, zpart1, zpart2, zsum1, zsum2, za, zb, zex
        REAL (kp) :: a1
        INTEGER :: nregion, nfv(3)

        zans = 0
! Determination of nregion
        nregion = num_region(znu,zz)
        IF (nregion==0) THEN
          info = 20
          RETURN
        END IF
        SELECT CASE (nregion)
        CASE (1) ! With the series expansion. ABS(zz)<1, ABS(zeps)>0.3.
          CALL bes1_series(znu,zz,zsum1,zlogbes1,info)
          IF (info>1) RETURN
          CALL bes1_series(-znu,zz,zsum2,zlogbes2,info)
          IF (info>1) RETURN
          IF (ai_znu>=0) THEN ! ai_znu>=0. See Eq. (7a) of Ref. (1).
            zarg1 = zlogbes1 + 2*znupi_i
            zarg2 = zlogbes2 + znupi_i
            IF (real(zarg1)>alog_huge .OR. real(zarg2)>alog_huge) THEN
              zans = huge1
              info = 10
              RETURN
            END IF
            zpart1 = zsum1*exp(zarg1)
            zpart2 = -zsum2*exp(zarg2)
            zex = exp(znupi_i)
            zans = 2*(zpart1+zpart2)/(zex**2-1)
            IF (abs1(zans)<tiny1*1E4) THEN
              info = 10
              RETURN
            END IF
          ELSE ! ai_znu<0. See Eq. (7b) of Ref. (1).
            zarg1 = zlogbes1
            zarg2 = zlogbes2 - znupi_i
            IF (real(zarg1)>alog_huge .OR. real(zarg2)>alog_huge) THEN
              zans = huge1
              info = 10
              RETURN
            END IF
            zpart1 = zsum1*exp(zarg1)
            zpart2 = -zsum2*exp(zarg2)
            za = exp(-znupi_i)
            zans = 2*(zpart1+zpart2)/(1-za**2)
          END IF
        CASE (2) ! With the series expansion. ABS(zz)<1, ABS(zeps)>=0.3.
          IF (re_znu>=0) THEN ! re_znu>=0. See Eq. (11a) of Ref. (1).
            CALL bes1_series(znu,zz,zsum1,zlogbes1,info)
            IF (info>1) RETURN
            IF (real(zlogbes1)>alog_huge) THEN
              zans = huge1
              info = 10
              RETURN
            END IF
            zbes1 = zsum1*exp(zlogbes1)
            CALL bes2_series(znu,zz,zbes2,info)
            IF (info>1) RETURN
            zans = zbes1 + cmplx(aimag(zbes2),-real(zbes2),kp)
            RETURN
          END IF ! re_znu<0. See Eq. (11b) of Ref. (1).
          CALL bes1_series(znua,zz,zsum1,zlogbes1,info)
          IF (info>1) RETURN
          zarg1 = zlogbes1 + znupi_i
          IF (real(zarg1)>alog_huge) THEN
            zans = huge1
            info = 10
            RETURN
          END IF
          CALL bes2_series(znua,zz,zbes2,info)
          IF (info>1) RETURN
          zex = exp(znupi_i)
          za = zex*zbes2
          zans = exp(zarg1)*zsum1 + cmplx(aimag(za),-real(za),kp)
        CASE (3) ! By numerical integration. ABS(zz)>=1.
          IF (ai_znu<-ai_znu_m) THEN
            zans = huge1
            info = 10
            RETURN
          END IF
          zex = exp(znupi_i)
          IF (re_znu>=0) THEN
            IF (re_zz>=0) THEN ! re_znu>=0, re_zz>=0.
              CALL cylin_inte(znu,zz,zbes1_t,zhan1_t,zhan2_t,nfv,info)
              IF (info>1) RETURN
              zans = zhan2_t*zex**nfv(3)
            ELSE ! re_znu>=0, re_zz<0.
              CALL cylin_inte(znu,zza,zbes1_t,zhan1_t,zhan2_t,nfv,info)
              IF (info>1) RETURN
              IF (ai_zz>=0) THEN
! re_znu>=0, re_zz<0, ai_zz>=0. See Eq. (41a) of Ref. (1).
                IF (exponent(abs1(2*zbes1_t))*alog2-ai_znupi*(nfv(1)+1)> &
                    alog_huge) THEN
                  zans = huge1
                  info = 10
                  RETURN
                END IF
                IF (exponent(abs1(zhan2_t))*alog2-ai_znupi*(nfv( &
                    3)-1)>alog_huge) THEN
                  zans = huge1
                  info = 10
                  RETURN
                END IF
                IF (nfv(3)-1>=0) THEN
                  zb = zex**(nfv(3)-1)
                ELSE
                  IF (ai_znu>ai_znu_m) THEN
                    info = 10
                    RETURN
                  END IF
                  zb = exp(-znupi_i)**(1-nfv(3))
                END IF
                zans = 2*zbes1_t*zex**(nfv(1)+1) + zhan2_t*zb
              ELSE ! re_znu>=0, re_zz<0, ai_zz<0. See Eq. (41b) of Ref. (1).
                a1 = exponent(abs1(zhan1_t))*alog2 - ai_znupi
                IF (a1>alog_huge) THEN
                  zans = huge1
                  info = 10
                  RETURN
                END IF
                zans = -zhan1_t*zex**(nfv(2)+1)
                IF (abs1(zans)<5*tiny1) info = 10
              END IF
            END IF
          ELSE ! re_znu<0
            IF (re_zz>=0) THEN ! re_znu<0, re_zz>=0. See Eq. (41c) of Ref. (1).
              CALL cylin_inte(znua,zz,zbes1_t,zhan1_t,zhan2_t,nfv,info)
              IF (info>1) RETURN
              nfv = -nfv
              a1 = exponent(abs1(zhan2_t))*alog2 - ai_znupi*(nfv(3)+1)
              IF (a1>alog_huge) THEN
                zans = huge1
                info = 10
                RETURN
              END IF
              zans = zhan2_t*zex**(nfv(3)+1)
              IF (abs1(zans)<10*tiny1) info = 10
            ELSE ! re_znu<0, re_zz<0
              CALL cylin_inte(znua,zza,zbes1_t,zhan1_t,zhan2_t,nfv,info)
              IF (info>1) RETURN
              nfv = -nfv
              IF (ai_zz>=0) THEN
! re_znu<0, re_zz<0, ai_zz>=0. See Eq. (41d) of Ref. (1).
                IF (exponent(abs1(zhan2_t))*alog2-ai_znupi*(nfv( &
                    3)+2)>alog_huge) THEN
                  zans = huge1
                  info = 10
                  RETURN
                END IF
                IF (-ai_znupi*(nfv(3)+2)>alog_huge) THEN
                  zans = huge1
                  info = 10
                  RETURN
                END IF
                zans = 2*zbes1_t*zex**nfv(1) + zhan2_t*zex**(nfv(3)+2)
              ELSE ! re_znu<0, re_zz<0, ai_zz<0. See Eq. (41e) of Ref. (1).
                zans = -zhan1_t*zex**nfv(2)
              END IF
            END IF
          END IF
        END SELECT
      END SUBROUTINE hankel2

      SUBROUTINE bes1_series(znu,zz,zsum,zlogbes,info)
        COMPLEX (kp), INTENT (IN) :: znu, zz
        COMPLEX (kp), INTENT (OUT) :: zsum, zlogbes
        INTEGER, INTENT (OUT) :: info
! This SUBROUTINE calculates zsum and zlogbes by means of the series expansion
! described in 9.1.10 of Ref. (9) when nregion is 1 or 2.
! zbessel1(znu,zz)=zsum*EXP(zlogbes).
! This calculates zbessel1(znu,zz) when nregion=1 and when nregion=2 and
! re_znu>=0.
! The argument info informs us of the output condition.
! The named constants kp, pi, tiny1, epsilon0 and a variable epsilon1 are 
! declared in MODULE mod_zbes.
! This is invoked by bessel1, bessel2, hankel1 and hankel2.
! This subroutine invokes internal subroutines zgamma, subgam.
        COMPLEX (kp) :: za, zloggam, zfct, zxx2
        INTEGER :: i

        zlogbes = 0
        info = 0
        zxx2 = (zz/2)**2
        za = 1
        zsum = za
        DO i = 1, 500
          za = -za*zxx2/(i*(znu+i))
          zsum = za + zsum
          IF (abs1(za)<epsilon1) EXIT
        END DO
        CALL zgamma(znu+1,zfct,zloggam,info)
        zsum = zsum/zfct
        IF (info>5) RETURN
        IF (abs1(zz)<tiny1) THEN
          IF (real(znu)>0) THEN
            zsum = 0
            zlogbes = 0
            RETURN
          END IF
          IF (abs1(znu)<tiny1) THEN
            zsum = 1
            zlogbes = 0
            RETURN
          END IF
          zsum = 1
          zlogbes = 20
          info = 10
          RETURN ! re_znu<=0,  ai_znu/=0
        END IF
        zlogbes = znu*log(zz/2) - zloggam
        IF (abs1(zlogbes)>ai_arg_m) info = 10
      CONTAINS

        SUBROUTINE zgamma(zx,zfct,zloggam,info)
! This subroutine and subroutine subgam compute the gamma function gamma(zx).
! These subroutines are written on the basis of Ref. (2).
! gamma(zx)=zfct*EXP(zloggam).
          COMPLEX (kp), INTENT (IN) :: zx
          COMPLEX (kp), INTENT (OUT) :: zfct, zloggam
          INTEGER, INTENT (OUT) :: info
          COMPLEX (kp) :: zx1, zx2, zlog_zsn, z1, zfct1
          REAL (kp), PARAMETER :: pi2 = pi*2
          INTEGER :: i2

          zfct = 1
          zloggam = 0
          info = 0
          i2 = 1
          IF (real(zx)<0.5) THEN
! gamma(zx)=pi/(SIN(pi*zx)*gamma(1-zx))
!          =2*pi*zunit/(EXP(zunit*pi*zx)*(1-EXP(-2*zunit*zx))*gamma(1-zx))
            CALL subgam(1-zx,zfct1,zloggam,info)
            IF (info/=0) RETURN
            zx1 = zx
            IF (abs1(zx1-nint(real(zx1)))<epsilon1) THEN
              info = 10
              zloggam = 0
              zfct = huge1
              RETURN
            END IF
            IF (aimag(zx)>0) THEN
              zx1 = conjg(zx)
              i2 = -1
            END IF
            zx2 = cmplx(-aimag(zx1),real(zx1),kp)*pi
            zlog_zsn = zx2
            z1 = cmplx(0,pi2,kp)/(1-exp(-2*zx2))
            IF (i2<0) THEN
              zlog_zsn = conjg(zlog_zsn)
              z1 = conjg(z1)
            END IF
            zfct = z1/zfct1
            zloggam = -zloggam - zlog_zsn ! SIN(pi*zx)/pi=EXP(zlog_zsn)/z1
            RETURN
          END IF
          CALL subgam(zx,zfct,zloggam,info)
        END SUBROUTINE zgamma

        SUBROUTINE subgam(zx,zfct,zgam,info)
! This subroutine is available only when REAL(zx)>=0.5.
! gamma(zx)=zfct*EXP(zgam).
          COMPLEX (kp), INTENT (IN) :: zx
          COMPLEX (kp), INTENT (OUT) :: zfct, zgam
          INTEGER, INTENT (OUT) :: info
          COMPLEX (kp) :: za, za1, zb, zf, zx1, zx2
          REAL (kp) :: coeff(7)
          REAL (kp), PARAMETER :: sq_pi2 = 2.5066282746310005024158_kp ! =SQRT(2*pi)
          REAL, PARAMETER :: dis = 20, dis2 = dis*dis
          INTEGER :: i
          DATA coeff/8.33333333333333333333E-2_kp, - &
            2.77777777777777777778E-3_kp, 7.93650793650793650794E-4_kp, &
            -5.95238095238095238095E-4_kp, 8.41750841750841750842E-4_kp, &
            -1.91752691752691752692E-3_kp, 6.41025641025641025641E-3_kp/

          info = 0
          zgam = 0
          zfct = 0
          zf = 1
          zx1 = zx
          DO
            IF (real(zx1)**2+aimag(zx1)**2>dis2) EXIT
            zf = zf*zx1
            zx1 = zx1 + 1
          END DO
          za = zx1
          zx2 = zx1**2
          zb = (zx1-.5_kp)*log(zx1) - zx1
          DO i = 1, 7
            za1 = coeff(i)/za
            zb = zb + za1
            IF (abs1(za1/zb)<epsilon1) EXIT
            IF (i>=7) THEN
              info = 20
              RETURN
            END IF
            za = za*zx2
          END DO
          zfct = sq_pi2/zf
          zgam = zb
        END SUBROUTINE subgam

      END SUBROUTINE bes1_series

      SUBROUTINE bes2_series(znu,zz,zbes2,info)
        COMPLEX (kp), INTENT (IN) :: znu, zz
        COMPLEX (kp), INTENT (OUT) :: zbes2
        INTEGER, INTENT (OUT) :: info
! This SUBROUTINE calculates zbes2 by the use of the series expansion when
! nregion=2. zbes2=zbessel2(znu,zz).
! This SUBROUTINE is written on the basis of Ref. (3).
! This is available only for re_znud>=0.  re_znud=REAL(znu).
! The internal subprograms of bes2_series are bes2_srs_init, def_bessel1, multi,
! multi_z, divis, divis_z.
! The named constants kp, huge1, ihuge, tiny1 are declared in MODULE mod_zbes.
! This is invoked by SUROUTINEs bessel1, bessel2, hankel1 and hankel2.
! This subroutine invokes subroutines bes2_srs_init and abs1.
        COMPLEX (kp) :: zeps, z1, z2, z3
        REAL (kp) :: a1, re_znud
        INTEGER :: i, nn

        zbes2 = 0
        info = 0
        IF (abs1(zz)<tiny1*2) THEN
          zbes2 = 1E10
          info = 10
          RETURN
        END IF
        re_znud = znu
        IF (abs(re_znud)>ihuge) THEN
          info = 20
          RETURN
        END IF
        nn = nint(re_znud)
        zeps = znu - nn
        IF (nn<=1) THEN !  n=0,1
          CALL bes2_srs_init(nn,zz,zbes2,info)
          IF (info>1) RETURN
        ELSE !  nn >= 2
! When nn >= 2, zbessel2(zeps,zz) and zbessel2(1+zeps,zz) are calculated first,
! and zbessel2(nn+zeps,zz) (n=2,3,...) are determined by the forward recurrence
! method. nn+zeps=znu. zeps is epsilon in Eq. (8) of Ref. (3).
          CALL bes2_srs_init(0,zz,z1,info)
          IF (info>1) RETURN
          CALL bes2_srs_init(1,zz,z2,info)
          IF (info>1) RETURN
          IF (nn>huge1*abs1(zz)/10) THEN
            info = 20
            RETURN
          END IF
          a1 = abs1(2*((nn-1)+zeps)/zz)*5
          z2 = z2/a1
          z1 = z1/a1
          DO i = 2, nn
            z3 = (2*((i-1)+zeps)/zz)*z2 - z1
            IF (abs1(z3)>huge1/a1) THEN
              zbes2 = huge1
              info = 10
              RETURN
            END IF
            z1 = z2
            z2 = z3
          END DO
          zbes2 = z3*a1
        END IF
      CONTAINS

        SUBROUTINE bes2_srs_init(nn,zz,zbes2,info)
          INTEGER, INTENT (IN) :: nn
          COMPLEX (kp), INTENT (IN) :: zz
          COMPLEX (kp), INTENT (OUT) :: zbes2
          INTEGER, INTENT (OUT) :: info
! This is invoked by bes2_series.
! This SUBROUTINE calculates zbes2 for nn=0,1 by the use of the series
! expansion.
! znu=nn+zeps. nn and zeps are defined by Eq. (8) of Ref. (3).
! zbes2=zbessel2(znu,zz).
! The named constants kp, huge1, pi, epsilon1 are declared in MODULE mod_zbes.
! The variable used in common with the host program: zeps.
! This is invoked by SUROUTINE bes2_series.
! This subroutine invokes subprograms def_bessel1 and abs1.
          COMPLEX (kp) :: zbes0, zbesp, zbesm, zcsbes, zepspi, zepsm, zeps2, &
            zdcos, zsin1, z0, z1, z2
          INTEGER :: i, nnm

          zbes2 = 0
          info = 0
          zepsm = -zeps
          nnm = -nn
          CALL def_bessel1(nn,zeps,zz,zbes0,zbesp,info)
          IF (info>1) RETURN
          CALL def_bessel1(nnm,zepsm,zz,z2,zbesm,info)
          IF (info>1) RETURN
          IF (nn>=1) zbesm = -zbesm
          zbes0 = zbes0 + zeps*zbesp ! zbes0=zbessel1(nn+zeps,zz)
! Computation of zdcos.
! zdcos is Df{cos(zeps*pi)} defined by Eq. (10) of Ref. (3).
          zepspi = pi*zeps
          zeps2 = zepspi**2
          z0 = -.5_kp
          z1 = z0
          DO i = 3, 55, 2
            z0 = -z0*zeps2/(i*(i+1))
            z1 = z1 + z0
            IF (abs1(z0)<epsilon1) EXIT
            IF (i>50) THEN
              info = 35
              RETURN
            END IF
          END DO
          zdcos = pi*zepspi*z1

          zcsbes = zdcos*zbes0
! zcsbes=zbessel1(nn+zeps,zz)*Df{cos(zeps*pi)}

! Computation of zsin1.
! zsin1=SIN(pi*zeps)/zeps. zsin1 is defined by Eq. (12) of Ref. (3).
          IF (abs1(zepspi)<1E-5) THEN
            zsin1 = pi*(1-zepspi**2/6)
          ELSE
            zsin1 = sin(zepspi)/zeps
          END IF
! The following statement is based on Eq. (9) of Ref. (3).
          zbes2 = (zcsbes+zbesp+zbesm)/zsin1
        END SUBROUTINE bes2_srs_init

        SUBROUTINE def_bessel1(nn,zeps,zz,zbes10,zbes11,info)
          INTEGER, INTENT (IN) :: nn
          COMPLEX (kp), INTENT (IN) :: zeps, zz
          COMPLEX (kp), INTENT (OUT) :: zbes10, zbes11
          INTEGER, INTENT (OUT) :: info
! znu=nn+zeps. nn and zeps are defined by Eq. (8) of Ref. (3).
! This SUBROUTINE calculates zbes10 and zbes11.
! zbessel1(nn+zeps,zz)=zbes10+zeps*zbes11. nn=-1,0,1
! zbessel1(nn     ,zz)=zbes10
! Def{zbessel1(nn+zeps,zz)}=zbes11
! Def{zbessel1(nn+zeps,zz)} is shown in Eq. (9) of Ref. (3).
! info: output condition.
!   info=0: normal output.
!   info=20: out of range.
!   info=35: failure in calculation.
! The named constants kp, tiny1, pi and variables epsilon1, alog_huge are
! declared in MODULE mod_zbes.
! This is invoked by SUROUTINE bes2_srs_init.
! This subroutine invokes subroutines abs1, multi, multi_z, divis, divis_z.
          COMPLEX (kp) :: za, zalogzz, za0, za1, zexp0, zexp1, zgamma1, zsig0, &
            zsig1, zd1, ze1, zb1, zs, zss, zeps3, zz2, z1, z2, z3, zd0, ze0, &
            zexp
          REAL (kp) :: amax, a0, a1, a2, gamma0, bmax, b2
          INTEGER :: i, i5, n1
          INTEGER, PARAMETER :: m3 = 26
          REAL (kp) :: da(m3)
          DATA da/0.5772156649015328606064_kp, -0.6558780715202538810765_kp, &
            -0.0420026350340952355289_kp, 0.1665386113822914895012_kp, &
            -0.0421977345555443367480_kp, -0.0096219715278769735622_kp, &
            0.0072189432466630995427_kp, -0.0011651675918590651118_kp, &
            -0.0002152416741149509743_kp, 0.0001280502823881161874_kp, &
            -0.0000201348547807882377_kp, -0.0000012504934821426730_kp, &
            0.0000011330272319816972_kp, -0.0000002056338416977615_kp, &
            0.0000000061160951044836_kp, 0.0000000050020076444708_kp, &
            -0.0000000011812745704992_kp, 0.0000000001043426711827_kp, &
            0.0000000000077822634465_kp, -0.0000000000036968056444_kp, &
            0.0000000000005100370461_kp, -0.0000000000000205832667_kp, &
            -0.0000000000000053481037_kp, 0.0000000000000012267886_kp, &
            -0.0000000000000001182577_kp, 0.0000000000000000013656_kp/

          info = 0
          zbes10 = 0
          zbes11 = 0
! Calculation of 1/Gamma(1+nn1+zeps). nn1=MAX(nn,0)=0,1. nn=-1,0,1
! gamma0+zeps*zgamma1=1/Gamma(1+nn1+zeps)
! Refer to Eqs. (22)-(25) of Ref. (3) and Ref. (4).
          gamma0 = 1
          zgamma1 = 0
          zeps3 = 1
          amax = 0
          DO i = 1, m3
            zs = da(i)*zeps3
            zgamma1 = zgamma1 + zs
            a1 = abs1(zs)
            amax = max(amax,a1)
            IF (a1/amax<epsilon1) EXIT
            IF (i==m3) THEN
              info = 35
              RETURN
            END IF
            zeps3 = zeps3*zeps
          END DO
          IF (nn>=1) THEN !  nn=nn1=1
            a0 = gamma0
            z1 = zgamma1
            CALL divis(zeps,a0,z1,1._kp,(1._kp,0._kp),gamma0,zgamma1)
          END IF

! Calculation of zsigma
! zsigma equals to sigma defined by Eqs. (17), (18a) of Ref. (3).
          zz2 = (zz/2)**2
          IF (nn>=0) THEN !  nn=0,1
            zsig0 = 1
            zsig1 = 0
            ze0 = 1
            n1 = 0
          ELSE !  nn=-1
            zsig0 = -zz2
            zsig1 = 1
            ze0 = -zz2
            n1 = 1
          END IF
          ze1 = 0
          amax = max(abs1(zsig0),tiny1)
          bmax = max(abs1(zsig1),tiny1)
          DO i = 1, 200
            n1 = n1 + 1
            a0 = nn + n1
            zb1 = 1
            zd0 = -ze0*zz2/n1
            zd1 = -ze1*zz2/n1
            CALL divis_z(zeps,zd0,zd1,a0,zb1,ze0,ze1)
            zsig0 = zsig0 + ze0
            zsig1 = zsig1 + ze1
            a2 = abs1(ze0)
            amax = max(a2,amax)
            b2 = abs1(ze1)
            bmax = max(b2,bmax)
            IF (a2/amax<epsilon1 .AND. b2/bmax<epsilon1) EXIT
            IF (i>=12) THEN
              info = 35
              RETURN
            END IF
          END DO
! zsigma=zsig0+zeps*zsig1
          CALL multi(zeps,cmplx(gamma0,0,kp),zgamma1,zsig0,zsig1,za0,za1)
! za0+zeps*za1=(zsig0+zeps*zsig1)/Gamma(1+nn1+zeps). nn1=MAX(nn,0)

! Calculation of zexp0, zexp1
! zexp0+zeps*zexp1=(zz/2)**(nn+zeps)
! Refer to Eqs. (20) and (21) of Ref. (3).
          z2 = (zz/2)**zeps
          z3 = z2 - 1
          zalogzz = log(zz/2)
          IF (abs1(z3)<0.5) THEN
            za = 1
            zss = zeps*zalogzz
            i5 = nint(aimag(zss)/(2*pi))
            zss = cmplx(real(zss),aimag(zss)-i5*2*pi,kp)
            zexp1 = 1
            DO i = 2, 50
              za = za*zss/i
              zexp1 = zexp1 + za
              IF (abs1(za)<epsilon1) EXIT
              IF (i==50) THEN
                info = 35
                RETURN
              END IF
            END DO
            z1 = 0
            IF (i5/=0) z1 = cmplx(0,i5*2*pi,kp)/zeps
            zexp1 = zexp1*(zalogzz-z1)
          ELSE
            zexp1 = z3/zeps
          END IF
          IF (real(zalogzz)*(nn-abs1(zeps))>alog_huge) THEN
            info = 20
            RETURN
          END IF
          zexp0 = (zz/2)**nn
          IF ((exponent(abs1(zz/2))*nn+exponent(abs1(zexp1)))*alog2>alog_huge) &
              THEN
            info = 10
            zbes10 = 0
            zbes11 = 0
            RETURN
          END IF
          zexp1 = zexp1*zexp0
          IF (nn==-1 .AND. abs1(zz)<1E-10 .AND. real(zeps*log(zz/2))<-5) THEN
            zexp = zexp0*z2
            CALL multi_z(za0,za1,zexp0,zexp1,zexp,zbes10,zbes11)
          ELSE
            CALL multi(zeps,za0,za1,zexp0,zexp1,zbes10,zbes11)
          END IF
        END SUBROUTINE def_bessel1

        SUBROUTINE multi(zeps,zin10,zin11,zin20,zin21,zout0,zout1)
          COMPLEX (kp), INTENT (IN) :: zeps, zin10, zin11, zin20, zin21
          COMPLEX (kp), INTENT (OUT) :: zout0, zout1
! This SUBROUTINE calculates zout0 and zout1.
! zout0+zeps*zout1=(zin10+zeps*zin11)*(zin20+zeps*zin21)
! Refer to Eqs. (6c) and (7c) of Ref. (3).
! The named constant kp is declared in MODULE mod_zbes.
! This is invoked by SUROUTINE def_bessel1.
          zout0 = zin10*zin20
          zout1 = zin10*zin21 + zin11*(zin20+zeps*zin21)
        END SUBROUTINE multi

        SUBROUTINE multi_z(zin10,zin11,zin20,zin21,zin22,zout0,zout1)
          COMPLEX (kp), INTENT (IN) :: zin10, zin11, zin20, zin21, zin22
          COMPLEX (kp), INTENT (OUT) :: zout0, zout1
! This SUBROUTINE calculates zout0 and zout1.
! zout0+zeps*zout1=(zin10+zeps*zin11)*zin22
! zin22=zin20+zeps*zin21
! Refer to Eqs. (6c) and (7c) of Ref. (3).
! The named constant kp is declared in MODULE mod_zbes.
! This is invoked by SUROUTINE def_bessel1.
          zout0 = zin10*zin20
          zout1 = zin10*zin21 + zin11*zin22
        END SUBROUTINE multi_z

        SUBROUTINE divis(zeps,ain1,zin1,ain2,zin2,aout,zout)
          COMPLEX (kp), INTENT (IN) :: zeps, zin1, zin2
          REAL (kp), INTENT (IN) :: ain1, ain2
          REAL (kp), INTENT (OUT) :: aout
          COMPLEX (kp), INTENT (OUT) :: zout
! This SUBROUTINE calculates aout and zout.
! aout+zeps*zout=(ain1+zeps*zin1)/(ain2+zeps*zin2)
! Refer to Eqs. (6d) and (7d) of Ref. (3).
! The named constant kp is declared in MODULE mod_zbes.
! This is invoked by SUROUTINE def_bessel1.
          aout = ain1/ain2
          zout = (zin1-(ain1/ain2)*zin2)/(ain2+zeps*zin2)
        END SUBROUTINE divis

        SUBROUTINE divis_z(zeps,zin0,zin1,aid0,zid1,zout0,zout1)
          COMPLEX (kp), INTENT (IN) :: zeps, zin0, zin1, zid1
          REAL (kp), INTENT (IN) :: aid0
          COMPLEX (kp), INTENT (OUT) :: zout0, zout1
! This SUBROUTINE calculates zout0 and zout1.
! zout0+zeps*zout1=(zin0+zeps*zin1)/(aid0+zeps*zid1)
! Refer to Eqs. (6d) and (7d) of Ref. (3).
! The named constant kp is declared in MODULE mod_zbes.
! This is invoked by SUROUTINE def_bessel1.
          zout0 = zin0/aid0
          zout1 = (zin1-(zin0/aid0)*zid1)/(aid0+zeps*zid1)
        END SUBROUTINE divis_z

      END SUBROUTINE bes2_series

      SUBROUTINE cylin_inte(znu,zz,zbes1_t,zhan1_t,zhan2_t,nfv,info)
        COMPLEX (kp), INTENT (IN) :: znu, zz
        COMPLEX (kp), INTENT (OUT) :: zbes1_t, zhan1_t, zhan2_t
        INTEGER, INTENT (OUT) :: nfv(3), info
! This SUBROUTINE calculates zbes1_t, zhan1_t zhan2_t and nfv by means of the
! numerical integration stated in Section 3.3 of Ref. (1) when nregion=3.
! zbessel1(znu,zz)=zbes1_t*zexd**nfv(1)
! zhankel1(znu,zz)=zhan1_t*zexd**nfv(2)
! zhankel2(znu,zz)=zhan2_t*zexd**nfv(3)
! zexd=EXP(pi*znu*zunit)
! This subroutine works under conditions of REAL(znu)>=0 and REAL(zz)>=0.
! The internal subprograms of this SUBROUTINE are integration, mposition,
! abs2, root1, zfun_stand, dfzero2, re_zfun2, newton, cal_phi, zdgaus8d and
! zfun.
! The named constants kp, tiny1, pi are declared in MODULE mod_zbes.
! This subroutine is invoked by SUBROUTINEs bessel1, bessel2, hankel1, hankel2.
! This subroutine invokes SUBROUTINEs integration, abs1.
! zgamma0 corresponds to gamma0 in Eq. (14) of Ref. (1).
! zinte corresponds to Ij in Eq. (19a) of Ref. (1).
! zt1 corresponds to t1 in Eqs. (34) of Ref. (1).
! radius0 corresponds to r0 in Eq. (21) of Ref. (1).
! jj corresponds to j in Eq. (15) of Ref. (1).
        COMPLEX (kp) :: zw0, zw0s, zw0t, zfun0, zfun00, znut, zam(3,4), zexd, &
          zgamma0, zinte, zt1, zinte1, zw1, z1, z2, zh, za
        REAL (kp) :: abs_re_zgamma, abs_re_zgamma_d, thetaz, phi_lim, a1, &
          pivot, radius, ai_zw0, am_area, tiny2
        REAL (kp), PARAMETER :: radius0 = 0.15, ab_er = 0.05
        INTEGER :: i, in, ip, is, j, jj, k, k1, l, mpsn1, mpsn2, nn, ind(2), &
          isd(3), id(3), kft(2,3,2,2)

! The integrals below mean the left-hand side of Eq. (19a) of Ref. (1).

! The starting point of the integral is given by fun(mpsn1), and
! the arriving point of the integral is given by fun(mpsn2),
! where mpsn1 and mpsn2 are integers, and function fun(mp) is given by TABLE A.

!                  TABLE A
!      mp                fun(mp)
!      11          inf+zunit*(-thetaz+2*pi)
!      10          inf+zunit*(-thetaz)
!       9          inf+zunit*(-thetaz-2*pi)
!       1         -inf+zunit*( thetaz+3*pi)
!       0         -inf+zunit*( thetaz+pi)
!      -1         -inf+zunit*( thetaz-pi)
!      -2         -inf+zunit*( thetaz-3*pi)
! inf=infinity.  zunit=(0,1).  thetaz=arg(zz).  mp is mpsn1 or mpsn2.

! When mpsn1 and mpsn2 are known, from Eqs. (12), (20) of Ref. (1) and TABLE A,
! the integral in Eq. (19a) of Ref. (1) is expressed as TABLE B.

!                  TABLE B
! mpsn1     mpsn2    Integral in Eq. (17a) of Ref. (1)
!   10       11              H1*zex2     +H2           =zinte
!    9       11              H1*(1+zex2) +H2*(1+1/zex2)=zinte
!    9       10              H1          +H2/zex2      =zinte
!    1       10     -2*J1*zex2           -H2           =zinte
!    0       11              H1*zex2                   =zinte
!    0       10                          -H2           =zinte
!    0        9             -H1          -H2*(1+1/zex2)=zinte
!    0        1      2*J1*zex2                         =zinte
!   -1       11              H1*(1+zex2) +H2           =zinte
!   -1       10              H1                        =zinte
!   -1        9                          -H2/zex2      =zinte
!   -1        1      2*J1*(1+zex2)                     =zinte
!   -1        0      2*J1                              =zinte
!   -2       11      2*J1/zex2+H1*(1+zex2)+H2          =zinte
!   -2       10      2*J1/zex2+H1                      =zinte
!   -2        0      2*J1*(1+1/zex2)                   =zinte
!   -2       -1      2*J1/zex2                         =zinte
! zex2=EXP(2*zunit*pi*znu).
! zinte=the value of the integral=Ij. Ij appears in Eq. (19a) of Ref. (1).
! J1=zbessel1(znu,zz),  H1=zhankel1(znu,zz),  H2=zhankel2(znu,zz)

! zam(jj,i)=kft(jj,i,1,1)*zex2**kft(jj,i,1,2)
!          +kft(jj,i,2,1)*zex2**kft(jj,i,2,2)
!   jj=1,2;  i=1,2,3.

! kft(jj,i,j,k):  jj=1,2= equation number. i=1: 2*J1;  i=2: H1;  i=3: H2.
!                 j=1,2= the term number.
!                 kft(jj,i,j,1)=-1,0,1;  kft(jj,i,j,2)=-2,-1,0,1,2.

!                TABLE C
! The 3-dimensional simultaneous linear equations for J1, H1 and H2 are:
! zam(1,1)*2*J1+zam(1,2)*H1+zam(1,3)*H2=zam(1,4)         (1)
! zam(2,1)*2*J1+zam(2,2)*H1+zam(2,3)*H2=zam(2,4)         (2)
!          2*J1         -H1         -H2=  0              (3)

! Eq. (1) corresponds to Eq. (19a) of Ref. (1) when j=1.
! Eq. (2) corresponds to Eq. (19a) of Ref. (1) when j=2.
! Eq. (3) corresponds to Eq. (19b) of Ref. (1).

! Correspondence between Eq. (19a) of Ref. (1) and Eqs. (1), (2) and (3).
! Eq. (19a) of Ref. (1)      Eqs. (1), (2), (3)
!           gj1                 zam(j,1)
!           gj2                 zam(j,2)
!           gj3                 zam(j,3)
!           Ij                  zam(j,4)
! Here, j=1, 2.

! The two paths c1 and c2 of numerical integration are defined in Section 3.3.3
! of Ref. (1).
! SUBROUTINE integration determines paths c1 and c2, and calculates zinte,
! mpsn1 and mpsn2.
! If zinte, mpsn1 and mpsn2 are known for paths c1 and c2,
! we can determine zam(jj,i) (jj=1,2; i=1,2,3,4) by using TABLEs B and C.
! Then, solving Eqs. (1), (2) and (3), we can determine J1, H1 and H2.

! zgamma0 below is gamma0 defined by the equation under Eq. (14) of Ref. (1).

        kft = 0
        zbes1_t = 0
        zhan1_t = 0
        zhan2_t = 0
        thetaz = atan2(aimag(zz),real(zz))
        znut = znu/zz
        za = sqrt(1-znut)*sqrt(znut+1)
        z1 = cmplx(-aimag(za),real(za),kp)
        zgamma0 = log(znut+z1)
        IF (real(zgamma0)<0) zgamma0 = -zgamma0
        abs_re_zgamma = abs(real(zgamma0))
        abs_re_zgamma_d = abs_re_zgamma + 1
        IF (abs(ai_znu)>ai_znu_m*0.565) THEN
          info = 10
          RETURN
        END IF
        IF (pi*abs(re_znu)>ai_arg_m) THEN
          info = 20
          RETURN
        END IF
        zexd = exp(cmplx(0,pi,kp)*znu)
        zam = 0
        tiny2 = tiny1/sqrt(epsilon1)
! To make 3-dimensional simultaneous linear equations.
! jj corresponds to j of Eq. (19a) of Ref. (1).
        DO jj = 1, 2
          CALL integration(zinte,info)
          IF (info>1) RETURN
          zam(jj,4) = zinte ! zinte: corresponds to Ij in Eq. (19a) of Ref. (1).
          SELECT CASE (mpsn1)
! From mpsn1 and mpsn2, zam(jj,k) (k=1,2,3) are determined by use of TABLE B.
          CASE (10)
            IF (mpsn2==11) THEN
              kft(jj,2,1,1) = 1
              kft(jj,2,1,2) = 1
              kft(jj,3,1,1) = 1
            ELSE
              info = 50
              RETURN
            END IF
          CASE (9)
            IF (mpsn2==11) THEN
              kft(jj,2,1,1) = 1
              kft(jj,2,1,2) = 0
              kft(jj,2,2,1) = 1
              kft(jj,2,2,2) = 1
              kft(jj,3,1,1) = 1
              kft(jj,3,1,2) = 0
              kft(jj,3,2,1) = 1
              kft(jj,3,2,2) = -1
            ELSE IF (mpsn2==10) THEN
              kft(jj,2,1,1) = 1
              kft(jj,3,1,1) = 1
              kft(jj,3,1,2) = -1
            ELSE
              info = 50
              RETURN
            END IF
          CASE (1)
            IF (mpsn2==10) THEN
              kft(jj,1,1,1) = -1
              kft(jj,1,1,2) = 1
              kft(jj,3,1,1) = -1
            ELSE
              info = 50
              RETURN
            END IF
          CASE (0)
            IF (mpsn2==11) THEN
              kft(jj,2,1,1) = 1
              kft(jj,2,1,2) = 1
            ELSE IF (mpsn2==10) THEN
              kft(jj,3,1,1) = -1
            ELSE IF (mpsn2==9) THEN
              kft(jj,2,1,1) = -1
              kft(jj,3,1,1) = -1
              kft(jj,3,1,2) = 0
              kft(jj,3,2,1) = -1
              kft(jj,3,2,2) = -1
            ELSE IF (mpsn2==1) THEN
              kft(jj,1,1,1) = 1
              kft(jj,1,1,2) = 1
            ELSE
              info = 50
              RETURN
            END IF
          CASE (-1)
            IF (mpsn2==11) THEN
              kft(jj,2,1,1) = 1
              kft(jj,2,1,2) = 0
              kft(jj,2,2,1) = 1
              kft(jj,2,2,2) = 1
              kft(jj,3,1,1) = 1
            ELSE IF (mpsn2==10) THEN
              kft(jj,2,1,1) = 1
            ELSE IF (mpsn2==9) THEN
              kft(jj,3,1,1) = -1
              kft(jj,3,1,2) = -1
            ELSE IF (mpsn2==1) THEN
              kft(jj,1,1,1) = 1
              kft(jj,1,1,2) = 0
              kft(jj,1,2,1) = 1
              kft(jj,1,2,2) = 1
            ELSE IF (mpsn2==0) THEN
              kft(jj,1,1,1) = 1
            ELSE
              info = 50
              RETURN
            END IF
          CASE (-2)
            IF (mpsn2==11) THEN
              kft(jj,1,1,1) = 1
              kft(jj,1,1,2) = -1
              kft(jj,2,1,1) = 1
              kft(jj,2,1,2) = 0
              kft(jj,2,2,1) = 1
              kft(jj,2,2,2) = 1
              kft(jj,3,1,1) = 1
            ELSE IF (mpsn2==10) THEN
              kft(jj,1,1,1) = 1
              kft(jj,1,1,2) = -1
              kft(jj,2,1,1) = 1
            ELSE IF (mpsn2==0) THEN
              kft(jj,1,1,1) = 1
              kft(jj,1,1,2) = 0
              kft(jj,1,2,1) = 1
              kft(jj,1,2,2) = -1
            ELSE IF (mpsn2==-1) THEN
              kft(jj,1,1,1) = 1
              kft(jj,1,1,2) = -1
            ELSE
              info = 50
              RETURN
            END IF
          CASE DEFAULT
            info = 51
            RETURN
          END SELECT
        END DO
! The elements zam(jj,k1) (jj,k1=1,2,3) are normalized by zexd**(2*nfv(k1))
! to avoid overflows of zam(jj,k1).
        IF (aimag(znu)>=0) THEN !  AIMAG(znu)>=0, |zexd|<<1
          nfv = 5
          DO k1 = 1, 3
            DO jj = 1, 2
              nfv(k1) = min(nfv(k1),kft(jj,k1,1,2),kft(jj,k1,2,2))
            END DO
          END DO
        ELSE ! AIMAG(znu)<0, |zexd|>>1
          nfv = -5
          DO k1 = 1, 3
            DO jj = 1, 2
              nfv(k1) = max(nfv(k1),kft(jj,k1,1,2),kft(jj,k1,2,2))
            END DO
          END DO
        END IF
        DO k1 = 1, 3
          DO jj = 1, 2
            kft(jj,k1,1,2) = kft(jj,k1,1,2) - nfv(k1)
            kft(jj,k1,2,2) = kft(jj,k1,2,2) - nfv(k1)
          END DO
        END DO
! AIMAG(znu)>=0, |zexd|<<1,    kft(jj,k1,1,2)>=0,   kft(jj,k1,2,2)>=0
! AIMAG(znu)<0,  |zexd|>>1,    kft(jj,k1,1,2)<=0,   kft(jj,k1,2,2)<=0
        z1 = 0
        z2 = 0
        k = 0
        DO jj = 1, 2
          DO k1 = 1, 3
            IF (kft(jj,k1,1,1)/=0) THEN
              IF (kft(jj,k1,1,2)>0) THEN
                z1 = zexd**(2*kft(jj,k1,1,2))
              ELSE IF (kft(jj,k1,1,2)==0) THEN
                z1 = 1
              ELSE
                IF ((2*kft(jj,k1,1,2))*pi*aimag(znu)>alog_huge) THEN
                  info = 10
                  RETURN
                END IF
                z1 = (1/zexd)**(-2*kft(jj,k1,1,2))
              END IF
              IF (abs1(z1)<tiny2) k = 1
            END IF
            IF (kft(jj,k1,2,1)/=0) THEN
              IF (kft(jj,k1,2,2)>0) THEN
                z2 = zexd**(2*kft(jj,k1,2,2))
              ELSE IF (kft(jj,k1,2,2)==0) THEN
                z2 = 1
              ELSE
                z2 = (1/zexd)**(-2*kft(jj,k1,2,2))
              END IF
              IF (abs1(z2)<tiny2) k = 1
            END IF
            IF (k==1) THEN
              info = 10
              RETURN
            END IF
            zam(jj,k1) = kft(jj,k1,1,1)*z1 + kft(jj,k1,2,1)*z2
          END DO
        END DO
        zam(3,1) = 1
        zam(3,2) = -1
        zam(3,3) = -1
        zam(3,4) = 0 ! Eq. (19b) of Ref. (1)
        DO k1 = 1, 3
          zam(3,k1) = zam(3,k1)*zexd**(-2*nfv(k1))
        END DO
        nfv = -2*nfv

! Re-arrangement of dimension zam(j,i)
        DO j = 1, 2
          is = 0
          DO i = 1, 3
            IF (abs1(zam(j,i))>=tiny1) THEN
              is = is + 1
              in = i
            END IF
          END DO
          IF (is==0) THEN
            info = 27
            RETURN
          END IF
          isd(j) = is
          ind(j) = in
        END DO
        isd(3) = 3
        IF (isd(1)>=2 .AND. isd(2)==1) THEN
          DO i = 1, 4
            zw1 = zam(1,i)
            zam(1,i) = zam(2,i)
            zam(2,i) = zw1
          END DO
          k = isd(1)
          isd(1) = isd(2)
          isd(2) = k
          k = ind(1)
          ind(1) = ind(2)
          ind(2) = k
        END IF
        id = (/ 1, 2, 3 /)
        DO i = 1, 2
          IF (isd(i)==1 .AND. ind(i)/=i) THEN
            DO j = 1, 3
              zw1 = zam(j,i)
              zam(j,i) = zam(j,ind(i))
              zam(j,ind(i)) = zw1
            END DO
            k = id(i)
            id(i) = id(ind(i))
            id(ind(i)) = k
            DO k = 1, 2
              IF (k==i) CYCLE
              IF (ind(k)==i) THEN
                ind(k) = ind(i)
                CYCLE
              END IF
              IF (ind(k)==ind(i)) ind(k) = i
            END DO
            ind(i) = i
          END IF
        END DO

! To solve the 3-dimensional simultaneous linear equations by the use of
! Gauss-Jordan's method.
        info = 0
        DO nn = 1, 3
! Pivoting
          IF (isd(nn)/=1) THEN
            pivot = 0
            ip = 1
            DO i = nn, 3
              a1 = abs1(zam(i,nn))
              IF (a1>pivot) THEN
                pivot = a1
                ip = i
              END IF
            END DO
            IF (pivot<=0) THEN
              info = 52
              RETURN
            END IF
            DO j = nn, 4
              zw1 = zam(nn,j)
              zam(nn,j) = zam(ip,j)
              zam(ip,j) = zw1
            END DO
          END IF
! Elimination
          DO j = nn + 1, 4
            zam(nn,j) = zam(nn,j)/zam(nn,nn)
          END DO
          DO i = 1, 3
            IF (i==nn) CYCLE
            DO j = nn + 1, 4
              zam(i,j) = zam(i,j) - zam(i,nn)*zam(nn,j)
            END DO
          END DO
        END DO
! Substitution of the solutions.
        DO i = 1, 5
          l = 0
          DO j = 1, 2
            IF (id(j)>id(j+1)) THEN
              l = l + 1
              zw1 = zam(j,4)
              zam(j,4) = zam(j+1,4)
              zam(j+1,4) = zw1
              k = id(j)
              id(j) = id(j+1)
              id(j+1) = k
            END IF
          END DO
          IF (l==0) EXIT
        END DO
        zbes1_t = zam(1,4)/2
        zhan1_t = zam(2,4)
        zhan2_t = zam(3,4)
      CONTAINS

        SUBROUTINE integration(zinte,info)
! This subroutine determines the paths cj, and calculates zinte, mpsn1
! and mpsn2.
! This subroutine invokes subprograms mposition, root1, dfzero2, newton,
! zdgaus8d, cal_phi and zfun.
! info: output condition
!   info=0: normal output.
!   info>1: failure in calculation
! Named constants kp, huge1, tiny1, ai_arg_m, alog2, pi and variables zex,
! alog_huge, alog_eps1, ai_znu_m, re_znu, ai_znu, re_zz, ai_zz are declared in
! MODULE mod_zbes.
! The variables used in common with the host program: radius, radius0, 
! zgamma0, jj, zfun00, zfun0, phi_lim, zw0, zw0s, zw0t, zinte1, am_area, 
! abs_re_zgamma, abs_re_zgamma_d, ai_zw0, thetaz, mpsn1, mpsn2.
! This subroutine is invoked by SUBROUTINE cylin_inte.
! This subroutine invokes SUBROUTINEs newton, root1, cal_phi, dfzero2, zdgaus8d,
! abs2, abs1.
          COMPLEX (kp), INTENT (OUT) :: zinte
          INTEGER, INTENT (OUT) :: info
          COMPLEX (kp) :: za, zgamma
          REAL (kp) :: phi0, phi1, phi2, phi_gam, radius_squ, dis1
          REAL (kp), PARAMETER :: phi5 = -1.571
          INTEGER :: i, kk, ll, lin

          radius_squ = (radius0*2.1)**2
          zgamma = zgamma0
          phi0 = 0
          IF (jj==2) THEN
            zgamma = -zgamma0
            phi0 = 3
          END IF
          zfun00 = zfun(zgamma)
          IF (real(zfun00)>alog_huge) THEN
            zinte = huge1
            info = 10
            RETURN
          END IF
          IF (real(zfun00)<alog_tiny) THEN
            zinte = tiny1
            info = 10
            RETURN
          END IF
          IF (abs1(zfun00)>ai_arg_m) THEN
            info = 20
            RETURN
          END IF

! zgamma corresponds to gammaj defined in Section 3.3.2 of Ref. (1).
! The following statements determine a descent curve from saddle point zgamma
! and calculate the integral zinte along this descent curve, which
! exists on the path cj, where suffix j is defined in Section 3.3.2 of Ref. (1)
! and corresponds to jj in this program.
! lin=0: the numerical integration is executed.
! lin=1: the numerical integration is not executed.
! ll=0: the path cj is determined by the steepest descent.
! ll=1: the path cj diverges from the steepest descent.
! ll corresponds to l defined in Step 4 of Section 3.3.3 of Ref. (1).
! kk corresponds to -k of Eqs. (21) of Ref. (1).
          zw0 = zgamma
          zfun0 = zfun00
          phi_lim = 0.785398 !=pi/4
          ai_zw0 = aimag(zw0)
          zinte = 0
          am_area = 0
          radius = radius0
          dis1 = 1E5
          DO i = -1, 1
            dis1 = min(abs(zgamma+cmplx(0,pi*i,kp)),dis1)
          END DO
          IF (radius0<5*dis1 .AND. dis1<1.1*radius0) radius = dis1/1.1
          CALL root1(phi0,phi1,info)
          IF (info>1) RETURN
          phi0 = phi1
          phi_gam = phi1
          radius = radius0
          ll = 0
          lin = 0
          DO kk = 1, 200
            IF (ll/=1 .AND. real(zw0)<abs_re_zgamma*0.8) THEN
              IF (ai_zw0>thetaz+2*pi) THEN
                zt1 = cmplx(-abs_re_zgamma_d-radius0,thetaz+3*pi)
                ll = 1
              ELSE IF (ai_zw0<thetaz-2*pi) THEN
                zt1 = cmplx(-abs_re_zgamma_d-radius0,thetaz-3*pi)
                ll = 1
              END IF
            END IF
            DO i = -1, 1
              za = zw0 + zgamma + cmplx(0,2*pi*i) ! za=zw0-(-zgamma-CMPLX(0,2*pi*i))
              IF (abs2(za)>radius_squ) CYCLE
              CALL root1(phi0+phi5,phi1,info)
              IF (info>1) RETURN
              phi0 = phi1
              GO TO 12
            END DO
            DO i = -1, 1, 2
              za = zw0 - zgamma - cmplx(0,2*pi*i) ! za=zw0-(zgamma+CMPLX(0,2*pi*i))
              IF (abs2(za)>radius_squ) CYCLE
              CALL root1(phi0+phi5,phi1,info)
              IF (info>1) RETURN
              phi0 = phi1
              GO TO 12
            END DO
            CALL newton(phi0,phi1,info)
            IF (info>1) RETURN
            phi0 = phi1
12          CONTINUE
            phi2 = phi0
            IF (ll==1) CALL cal_phi(zw0,phi0,phi2)
            zw0s = zw0
            zw0 = zw0 + radius0*cmplx(cos(phi2),sin(phi2),kp)
            ai_zw0 = aimag(zw0)
            zfun0 = zfun(zw0)
            zw0t = zw0
            IF (lin==0) THEN
              IF (real(zfun0-zfun00)<alog_eps1) THEN
                CALL dfzero2(info)
                IF (info>1) RETURN
                lin = 1
              END IF
              CALL zdgaus8d(info)
              IF (info>1) RETURN
            END IF
            IF (lin>0) THEN
              phi_lim = 1.4137 !=0.45*pi
              mpsn1 = mposition()
              IF (mpsn1<12) EXIT
            END IF
            IF (kk>150) THEN
              info = 80
              RETURN
            END IF
          END DO
          zinte1 = zinte

! The following statements determine another descent curve from saddle point
! zgamma and calculate the integral zinte along this descent curve, which
! exists on the curve cj, where suffix j is defined in Section 3.3.2 of Ref. (1)
! and corresponds to jj in this program.
! lin=0: the numerical integration is executed.
! lin=1: the numerical integration is not executed.
! ll=0: the path cj is determined by the steepest descent.
! ll=1: the path cj diverges from the steepest descent.
! ll corresponds to l defined in Step 4 of Section 3.3.3 of Ref. (1).
! kk corresponds to +k of Eqs. (21) of Ref. (1).
          phi_lim = 0.785398
          zw0 = zgamma
          zfun0 = zfun00
          ai_zw0 = aimag(zw0)
          zinte = 0
          phi0 = phi_gam + pi
          ll = 0
          lin = 0
          DO kk = 1, 200
            IF (ll/=1 .AND. real(zw0)<abs_re_zgamma*0.8) THEN
              IF (ai_zw0>thetaz+2*pi) THEN
                zt1 = cmplx(-abs_re_zgamma_d-radius0,thetaz+3*pi)
                ll = 1
              ELSE IF (ai_zw0<thetaz-2*pi) THEN
                zt1 = cmplx(-abs_re_zgamma_d-radius0,thetaz-3*pi)
                ll = 1
              END IF
            END IF
            DO i = -1, 1
              za = zw0 + zgamma + cmplx(0,2*pi*i)
              IF (abs2(za)>radius_squ) CYCLE
              CALL root1(phi0+phi5,phi1,info)
              IF (info>1) RETURN
              phi0 = phi1
              GO TO 22
            END DO
            DO i = -1, 1
              za = zw0 - zgamma - cmplx(0,2*pi*i)
              IF (abs2(za)>radius_squ) CYCLE
              CALL root1(phi0+phi5,phi1,info)
              IF (info>1) RETURN
              phi0 = phi1
              GO TO 22
            END DO
            CALL newton(phi0,phi1,info)
            IF (info>1) RETURN
            phi0 = phi1
22          CONTINUE
            phi2 = phi0
            IF (ll==1) CALL cal_phi(zw0,phi0,phi2)
            zw0s = zw0
            zw0 = zw0 + radius0*cmplx(cos(phi2),sin(phi2),kp)
            ai_zw0 = aimag(zw0)
            zfun0 = zfun(zw0)
            zw0t = zw0
            IF (lin==0) THEN
              IF (real(zfun0-zfun00)<alog_eps1) THEN
                CALL dfzero2(info)
                IF (info>1) RETURN
                lin = 1
              END IF
              CALL zdgaus8d(info)
              IF (info>1) RETURN
            END IF
            IF (lin>0) THEN
              phi_lim = 1.4137
              mpsn2 = mposition()
              IF (mpsn2<12) EXIT
            END IF
            IF (kk>150) THEN
              info = 80
              RETURN
            END IF
          END DO
          IF (mpsn1==mpsn2) THEN
            info = 88
            RETURN
          END IF
          za = zinte - zinte1
          zinte = cmplx(aimag(za),-real(za),kp)/pi
! Inversion. This corresponds to Eq. (20c) of Ref. (1).
          IF (mpsn1>mpsn2) THEN
            i = mpsn2
            mpsn2 = mpsn1
            mpsn1 = i
            zinte = -zinte
          END IF
        END SUBROUTINE integration

        FUNCTION mposition() RESULT (mpt)
! This function determines the position numbers mpsn1 and mpsn2.
! The named constants kp, pi are declared in MODULE mod_zbes.
! The variables used in common with the host program: ai_zw0, thetaz, zw0,
! abs_re_zgamma_d.
! This subroutine is invoked by SUBROUTINE integration.
          INTEGER :: i, mpt

          mpt = 15
          IF (real(zw0)<-abs_re_zgamma_d) THEN
            DO i = -2, 1
              IF (abs(ai_zw0-(thetaz+pi+i*2*pi))<pi) THEN
                mpt = i
                RETURN
              END IF
            END DO
          END IF
          IF (real(zw0)>abs_re_zgamma_d) THEN
            DO i = -1, 1
              IF (abs(ai_zw0-(-thetaz+i*2*pi))<pi) THEN
                mpt = i + 10
                RETURN
              END IF
            END DO
          END IF
          IF (ai_zw0>thetaz+3*pi) THEN
            mpt = 1
          ELSE IF (ai_zw0<thetaz-3*pi) THEN
            mpt = -2
          END IF
        END FUNCTION mposition

        FUNCTION abs2(za) RESULT (ab)
! FUNCTION abs2 = ABS(za)**2.
! The named constants kp, pi are declared in MODULE mod_zbes.
! This subroutine is invoked by SUBROUTINE integration.
          COMPLEX (kp), INTENT (IN) :: za
          REAL (kp) :: ab

          ab = real(za)**2 + aimag(za)**2
        END FUNCTION abs2

        SUBROUTINE root1(phi0,phi1,info)
! This subroutine calculates a root of Eqs. (22), (30), (31a), (36) of Ref. (1).
! This subroutine corresponds to Procedure A of Ref. (1).
! First this subroutine computes a rough value of a root of Eqs. (22), (30), 
! (31a), (36) of Ref. (1) by Eqs. (37) of Ref. (1), and then this computes the
! precise root from the rough root using a combination of bisection and the
! secant rule. This subroutine was written by modifying SUBROUTINE DFZERO, which
! was downloaded from the site http://gams.nist.gov/
! SUBROUTINE DFZERO was made by Shampine, L. F., (SNLA) and Watts, H. A., 
! (SNLA). SUBROUTINE DFZERO consists of a combination of bisection and the
! secant rule.
! zw=zw0+radius*CMPLX(COS(phi),SIN(phi))
! phi0: a first approximation of phi (input).
! phi1: the output value of phi (output).
! info: the output condition.
!   info=0 : normal output.
!   info=62: a root is not found.
! The named constants kp, pi are declared in MODULE mod_zbes.
! This subroutine is invoked by SUBROUTINE integration.
! This subroutine invokes SUBROUTINE zfun_stand.
! del_phi corresponds to delta phi in Eq. (37d) of Ref. (1).
! radius corresponds to r in Eq. (37c) of Ref. (1).
          REAL (kp), INTENT (IN) :: phi0
          REAL (kp), INTENT (OUT) :: phi1
          INTEGER, INTENT (OUT) :: info
          REAL (kp) :: pha, phb, re_fa, re_fb, ai_fa, ai_fb
          REAL (kp) :: phc, ai_fc
          REAL (kp) :: acbs, acmb, cmb, p, q
          REAL (kp), PARAMETER :: ab_er1 = 0.02, ab_er2 = ab_er1*0.3, &
            del_phi = 0.2
          INTEGER :: i, ic

          info = 0
          phi1 = 0
          pha = phi0
          phb = pha + del_phi
          CALL zfun_stand(pha,re_fa,ai_fa)
          CALL zfun_stand(phb,re_fb,ai_fb)
          DO
            IF (pha>phi0+2.0001*pi) THEN
              info = 62
              RETURN
            END IF
            IF (ai_fa*ai_fb<=0 .AND. (re_fa<=0 .OR. re_fb<=0)) THEN
              phi1 = 0
              phc = phb
              ai_fc = ai_fb
              ic = 0
              acbs = 0.5_kp*abs(pha-phb)
              DO i = 1, 500
                IF (abs(ai_fb)<abs(ai_fa)) THEN ! Perform interchange.
                  phc = pha
                  ai_fc = ai_fa
                  pha = phb
                  ai_fa = ai_fb
                  phb = phc
                  ai_fb = ai_fc
                END IF
                cmb = 0.5_kp*(phb-pha)
                acmb = abs(cmb)
! Test stopping criterion. Process results for proper setting of info.
                IF (acmb<ab_er1 .OR. abs(ai_fa)<ab_er2) THEN
                  info = 0
                  phi1 = pha
                  EXIT
                END IF
                IF (i>80) THEN
                  info = 63
                  RETURN
                END IF
! Calculate new iterate implicitly as pha+p/q, where we arrange p > 0.
! The implicit form is used to prevent overflow.
                p = (pha-phc)*ai_fa
                q = ai_fc - ai_fa
                IF (p<0) THEN
                  p = -p
                  q = -q
                END IF
! Update phc and check for satisfactory reduction in the size of the
! bracketing interval.  If not, perform bisection.
                phc = pha
                ai_fc = ai_fa
                ic = ic + 1
                IF (ic>=4) THEN
                  IF (4*acmb>acbs) GO TO 8
                  ic = 0
                  acbs = acmb
                END IF
! Root ought to be between pha and (phb+pha)/2.
                IF (p>cmb*q) GO TO 8
                pha = pha + p/q
                GO TO 9 ! Use secant rule.
8               pha = pha + cmb ! Use bisection (phb+pha)/2.
9               CALL zfun_stand(pha,re_fa,ai_fa)
! Decide whether next step is interpolation or extrapolation.
                IF (sign(1.0_kp,ai_fa)*sign(1.0_kp,ai_fb)>0) THEN
                  phb = phc
                  ai_fb = ai_fc
                END IF
              END DO
              phi1 = mod(phi1,2*pi)
              RETURN
            END IF
            pha = phb
            phb = phb + del_phi
            re_fa = re_fb
            ai_fa = ai_fb
            CALL zfun_stand(phb,re_fb,ai_fb)
          END DO
        END SUBROUTINE root1

        SUBROUTINE zfun_stand(phi,re_zfun_st,ai_zfun_st)
! The named constant kp is declared in MODULE mod_zbes.
! The variables used in common with the host program: radius, zfun0, zw0.
! This is invoked by SUBROUTINE root1.
! This subroutine invokes FUNCTION zfun.
! radius corresponds to r in Eq. (37c) of Ref. (1).
          REAL (kp), INTENT (IN) :: phi
          REAL (kp), INTENT (OUT) :: re_zfun_st, ai_zfun_st
          COMPLEX (kp) :: zfun_st, zw

          zw = zw0 + radius*cmplx(cos(phi),sin(phi),kp)
          zfun_st = zfun(zw) - zfun0
          re_zfun_st = real(zfun_st)
          ai_zfun_st = aimag(zfun_st)
        END SUBROUTINE zfun_stand

        SUBROUTINE dfzero2(info)
! This subroutine solves Eq. (26) of Ref. (1).
! This subroutine obtains root ss of the equation:
! REAL(zfun0-zfun00)=alog_eps1+ab_er, where zfun0=zfun(ss*zh+zw0s), zh=zw0-zw0s
! and 0<=ss<=1. zw0t=ss*zh+zw0s.
! Variable ss here corresponds to s of Eq. (26) of Ref. (1).
! This subroutine is made by modifying SUBROUTINE DFZERO, which is downloaded
! from the site http://gams.nist.gov/
! SUBROUTINE DFZERO was made by Shampine, L. F., (SNLA) and Watts, H. A., 
! (SNLA). SUBROUTINE DFZERO consists of a combination of bisection and the
! secant rule.
! The named constant kp is declared in MODULE mod_zbes.
! The variables used in common with the host program: ab_er, zh, zw0, zw0s 
! and zw0t.
! This subroutine is invoked by the internal SUBROUTINE integration.
! This subroutine invokes the internal SUBROUTINE re_zfun2.
          INTEGER, INTENT (OUT) :: info
          REAL (kp) :: a, ss, c, acbs, acmb, cmb, fa, fss, fc, p, q
          INTEGER i, ic

          zh = zw0 - zw0s
          ss = 0
          c = 1
          fss = re_zfun2(ss)
          fc = re_zfun2(c)
          a = c
          fa = fc
          ic = 0
          acbs = 0.5_kp*abs(ss-c)
          DO i = 1, 500
            IF (abs(fc)<abs(fss)) THEN
              a = ss
              fa = fss
              ss = c
              fss = fc
              c = a
              fc = fa ! Perform interchange.
            END IF
            cmb = 0.5_kp*(c-ss)
            acmb = abs(cmb)
! Test stopping criterion. Process results for proper setting of info.
            IF (fss<=ab_er) THEN
              info = 0
              zw0t = ss*zh + zw0s
              RETURN
            END IF
            IF (i>80) THEN
              info = 80
              RETURN
            END IF
! Calculate new iterate implicitly as ss+p/q, where we arrange p > 0.
! The implicit form is used to prevent overflow.
            p = (ss-a)*fss
            q = fa - fss
            IF (p<0) THEN
              p = -p
              q = -q
            END IF
! Update a and check for satisfactory reduction in the size of the
! bracketing interval.  If not, perform bisection.
            a = ss
            fa = fss
            ic = ic + 1
            IF (ic>=4) THEN
              IF (4*acmb>acbs) GO TO 8
              ic = 0
              acbs = acmb
            END IF
! Root ought to be between ss and (c+ss)/2.
            IF (p>cmb*q) GO TO 8
            ss = ss + p/q
            GO TO 9 ! Use secant rule.
8           ss = ss + cmb ! Use bisection (c+ss)/2.
9           fss = re_zfun2(ss) ! Have completed computation for new iterate ss.
! Decide whether next step is interpolation or extrapolation.
            IF (sign(1.0_kp,fss)*sign(1.0_kp,fc)>0) THEN
              c = a
              fc = fa
            END IF
          END DO
        END SUBROUTINE dfzero2

        FUNCTION re_zfun2(dhh) RESULT (ref2)
! The named constants kp, alog_eps1 are declared in MODULE mod_zbes.
! The variables used in common with the host program: ab_er, zfun00, zh, zw0s.
! This is invoked by SUBROUTINE dfzero2.
! This subroutine invokes the internal FUNCTION zfun.
          REAL (kp), INTENT (IN) :: dhh
          REAL (kp) :: ref2

          ref2 = real(zfun(dhh*zh+zw0s)-zfun00) - alog_eps1 - ab_er ! alog_eps1=-36.04
        END FUNCTION re_zfun2

        SUBROUTINE newton(phi0,phi1,info)
! This solves Eqs. (31) of Ref. (1) by the Newton-Raphson method.
! phi1 corresponds to phi(jk) in Eqs. (31) of Ref. (1).
! This subroutine invokes FUNCTION zfun.
! phi0: (in) the first approximate root
! phi1: (out) the outputted root
! info: (out) output information
!   info=0   normal output
!   info=65  phi1 diverges.
!   info=66  excess of trial number.
! eps1: the upper limit of the error of the outputted root phi1.
! mm: the upper limit of the trial number.
! ai_fun: the function whose roots are intended to be obtained.
!   ai_fun=AIMAG(znu*zw-zz*SINH(zw)-zfun0)
!   zw=zw0+radius*CMPLX(COS(phip),SIN(phip))
! dif_ai_fun: the derivative of function ai_fun.  dif_ai_fun=d(ai_fun)/d(phip)
! The named constant kp is declared in MODULE mod_zbes.
! The variables used in common with the host program: znu, zz, zw0, radius,
! zfun0.
! This subroutine is invoked by SUBROUTINE integration.
          REAL (kp), INTENT (IN) :: phi0
          REAL (kp), INTENT (OUT) :: phi1
          INTEGER, INTENT (OUT) :: info
          COMPLEX (kp) :: zw, zwd
          REAL (kp) :: ai_fun, dif, dif_ai_fun, sint, cost, phip, ar, ai, sih, &
            coh, sit, cot
          REAL (kp), PARAMETER :: eps1 = 1E-3
          INTEGER :: i
          INTEGER, PARAMETER :: mm = 8

          phip = phi0
          DO i = 1, 500
            sint = sin(phip)
            cost = cos(phip)
            zw = zw0 + radius*cmplx(cost,sint,kp)
! zfun(zw)=znu*zw-zz*SINH(zw),  zfun0=zfun(zw0)
            ar = real(zw)
            ai = aimag(zw)
            sih = sinh(ar)
            coh = cosh(ar)
            sit = sin(ai)
            cot = cos(ai)
            ai_fun = aimag(znu*zw-zz*cmplx(sih*cot,coh*sit,kp)-zfun0)
            zwd = radius*cmplx(-sint,cost,kp) ! =d(zw)/d(phip)
            dif_ai_fun = aimag((znu-zz*cmplx(coh*cot,sih*sit,kp))*zwd)
            dif = ai_fun/dif_ai_fun
            phip = phip - dif
            IF (abs(phip)>15) THEN
              phi1 = phip
              info = 65
              RETURN
            END IF
            IF (i>mm) THEN
              phi1 = phip
              info = 66
              RETURN
            END IF
            IF (abs(dif)<eps1) THEN
              phi1 = phip
              info = 0
              RETURN
            END IF
          END DO
        END SUBROUTINE newton

        SUBROUTINE cal_phi(zw0,phi0,phi1)
! Using Eqs. (33) in Ref. (1), this subroutine determines phi when ll=1.
! phi corresponds to hat-phi(jk) defined in the first equation of Eqs. (33b).
! zw0: zw0 corresponds to w(j,tau) in Eq. (33a) of Ref. (1). (input)
! zt1: the target showing the direction of progression of the path of the
!       integration. zt1 corresponds to t1 in Eqs. (34) of Ref. (1). (input)
! phi0: the direction of steepest descent. (input)
! phi1: the output of phi. (output)
! The named constants kp, pi are declared in MODULE mod_zbes.
! The variables used in common with the host program: radius, phi_lim, zt1.
! This subroutine is invoked by SUBROUTINE integration.
          COMPLEX (kp), INTENT (IN) :: zw0
          REAL (kp), INTENT (IN) :: phi0
          REAL (kp), INTENT (OUT) :: phi1
          COMPLEX (kp) :: za, zt12
          REAL (kp) :: phi_d, abs_ph

          zt12 = zt1
          IF (real(zw0-zt1)<0) zt12 = zt1 - 3
          za = zt12 - zw0
          phi_d = atan2(aimag(za),real(za)) - phi0
          phi_d = mod(phi_d,2*pi)
          IF (phi_d>pi) phi_d = phi_d - 2*pi
          IF (phi_d<-pi) phi_d = phi_d + 2*pi
          abs_ph = abs(phi_d)
          phi1 = phi0 + sign(min(abs_ph,phi_lim),phi_d)
        END SUBROUTINE cal_phi

        SUBROUTINE zdgaus8d(info)
! Subroutine zdgaus8d is written by modifying subroutine DGAUS8, which was 
! produced by Jones, R. E. (SNLA). SUBROUTINE DGAUS8 was downloaded from site 
! http://gams.nist.gov/. 
! Subroutine DGAUS8 uses an adaptive 8-point Legendre-Gauss algorithm.
! The modification is as follows.
! (1) DGAUS8 is written with Fortran 77, and this is translated into Fortran 90.
! (2) DGAUS8 integrates a real function. This is modified so that zdgaus8d can
!     integrate a complex function on a segment of the complex plane zw.
! (3) The numerical values of x1,...,x4,w1,...w4 in SUBROUTINE zg8 were 
!     rewritten to increase the precision.
! (4) This subroutine was modified to conform to this module.

! Explanation of variables used in this subroutine.
! The integrand in this subroutine is EXP(zfun(zw)), where zw is the 
! independent variable and zfun(zw) is the function defined in an internal
! function of subroutine cylin_inte.
! zw corresponds to w in Eqs. (12) of Ref. (1).
! zw0s: the starting point of the integration.
! zw0t: the arriving point of the integration.
! The path of integration is the segment from zw0s to zw0t.
! zinte: the output for the integration.
! info: the output condition.
!   info= 0: normal output.
!   info=60: zinte does not converge.
! The named constants kp, epsilon1 are declared in MODULE mod_zbes.
! The variables used in common with the host program: zinte, zw0s, zw0t,
! am_area.
! This subroutine is invoked by SUBROUTINE integration.
! This subroutine invokes SUBROUTINE zg8.
          INTEGER, INTENT (OUT) :: info
          COMPLEX (kp) :: zaa(0:20), zhh(0:20), zgr(0:20), zvl(0:20), zgl, &
            zglr, zest, zvr
          REAL (kp) :: ef, tol
          REAL (kp), PARAMETER :: fct_ef = 2.5
          INTEGER :: l, lr(0:20), l1
! Initialize
          l = 0
          zhh(l) = (zw0t-zw0s)/4
          zaa(l) = zw0s
          lr(l) = 1
          CALL zg8(zaa(l)+2*zhh(l),2*zhh(l),zest)
          am_area = max(abs1(zest),am_area)
          ef = 1
          tol = am_area*epsilon1*50
! Compute refined estimates, estimate the error, etc.
20        CONTINUE
          CALL zg8(zaa(l)+zhh(l),zhh(l),zgl)
          CALL zg8(zaa(l)+3*zhh(l),zhh(l),zgr(l))
          zglr = zgl + zgr(l)
          IF (abs1(zest-zglr)>tol*ef) THEN
            l1 = l + 1
            DO l = l1, 20
              ef = fct_ef*ef
              zhh(l) = zhh(l-1)/2 ! 4*zhh(l)=(zw0t-zw0s)/2**l
! Consider the left half of this level
              zest = zgl
              lr(l) = -1
              zaa(l) = zaa(l-1)
              CALL zg8(zaa(l)+zhh(l),zhh(l),zgl)
              CALL zg8(zaa(l)+3*zhh(l),zhh(l),zgr(l))
              zglr = zgl + zgr(l)
              IF (abs1(zest-zglr)>tol*ef) CYCLE
              zvl(l) = zglr
! Proceed to right half at this level
              zest = zgr(l-1)
              lr(l) = 1
              zaa(l) = zaa(l) + 4*zhh(l)
              CALL zg8(zaa(l)+zhh(l),zhh(l),zgl)
              CALL zg8(zaa(l)+3*zhh(l),zhh(l),zgr(l))
              zglr = zgl + zgr(l)
              IF (abs1(zest-zglr)<tol*ef) EXIT
              IF (l>12) EXIT
            END DO
          END IF
          zvr = zglr
          l1 = l - 1
! Return one level
          DO l = l1, 0, -1
            ef = ef/fct_ef
            IF (lr(l)<0) THEN
              zvl(l) = zvl(l+1) + zvr
              zest = zgr(l-1)
              lr(l) = 1
              zaa(l) = zaa(l) + 4*zhh(l)
              GO TO 20
            END IF
            zvr = zvl(l+1) + zvr
          END DO
          zinte = zinte + zvr
          info = 0
        END SUBROUTINE zdgaus8d

        SUBROUTINE zg8(zx,zh,zg8t)
! This subroutine integrates zfun(zw) from zw=zx to zw=zx+zh using the 8-point
! Legendre-Gauss algorithm.
! The named constant kp is declared in MODULE mod_zbes.
! This subroutine is invoked by subroutine zdgaus8d.
          COMPLEX (kp), INTENT (IN) :: zx, zh
          COMPLEX (kp), INTENT (OUT) :: zg8t
          REAL (kp), PARAMETER :: x1 = 0.1834346424956498049394761_kp, &
            w1 = 0.3626837833783619829651504_kp, x2 = &
            0.5255324099163289858177390_kp, w2 = &
            0.3137066458778872873379622_kp, x3 = &
            0.7966664774136267395915539_kp, w3 = &
            0.2223810344533744705443560_kp, x4 = &
            0.9602898564975362316835609_kp, w4 = &
            0.1012285362903762591525314_kp

          zg8t = w1*(exp(zfun(zx-x1*zh))+exp(zfun(zx+x1*zh))) + &
            w2*(exp(zfun(zx-x2*zh))+exp(zfun(zx+x2*zh))) + &
            w3*(exp(zfun(zx-x3*zh))+exp(zfun(zx+x3*zh))) + &
            w4*(exp(zfun(zx-x4*zh))+exp(zfun(zx+x4*zh)))
          zg8t = zg8t*zh
        END SUBROUTINE zg8

        FUNCTION zfun(zw) RESULT (zfunt)
! Function zfun(zw) corresponds to z*f(w) where f(w) is the first equation of
! Eq. (13) of Ref. (1).
! zfun(zw)=zz*zw*COSH(zgamma)-zz*SINH(zw)=znu*zw-zz*SINH(zw)
! COSH(zgamma)=znu/zz
! zgamma corresponds to gammaj defined in Section 3.3.2 of Ref. (1).
! The named constant kp is declared in MODULE mod_zbes.
! The variables used in common with the host program: znu, zz.
! This subroutine is invoked by subroutines zg8, newton, re_zfun2, zfun_stand,
! integration.
          COMPLEX (kp), INTENT (IN) :: zw
          COMPLEX (kp) :: zfunt
          REAL (kp) :: ar, ai

          ar = real(zw)
          ai = aimag(zw)
          zfunt = znu*zw - zz*cmplx(sinh(ar)*cos(ai),cosh(ar)*sin(ai),kp)
        END FUNCTION zfun

      END SUBROUTINE cylin_inte

      FUNCTION num_region(znu,zz) RESULT (nregion)
        COMPLEX (kp), INTENT (IN) :: znu, zz
        INTEGER :: nregion
! This FUNCTION determines nregion. nregion=num_region(znu,zz)
! This FUNCTION determines the values of epsilon1, alog_huge, alog_tiny, 
! ai_znu_m, alog_eps1, znua, zza, re_znu, ai_znu, znupi, ai_znupi, re_zz, ai_zz,
! znupi_i.
! Named constants kp, huge1, tiny1, ihuge, epsilon0, ai_arg_m, pi, alog2 and
! variables epsilon1, alog_huge, alog_tiny, ai_znu_m, alog_eps1, znua, zza,
! znupi_i, re_znu, ai_znu, re_zz, ai_zz, ai_znupi are declared in MODULE
! mod_zbes.
! This is invoked only by subroutines bessel1, bessel2, hankel1, hankel2.
! epsilon1 corresponds to epsilon1 in Eq. (25) in Ref. (1).
        REAL (kp) :: abs_reznu

        IF (alog_huge<0) THEN
          epsilon1 = max(epsilon0,1E-20_kp) ! =2.22E-16
          alog_eps1 = log(epsilon1) ! =-36.04
          alog_huge = log(huge1/2) ! =707.99
          ai_znu_m = alog_huge/pi ! =224.
          alog_tiny = log(12*tiny1/epsilon1) ! =-669.86
          alog2 = LOG(2.)
        END IF
        znua = -znu
        zza = -zz
        re_znu = real(znu)
        ai_znu = aimag(znu)
        znupi = znu*pi
        ai_znupi = aimag(znupi)
        re_zz = real(zz)
        ai_zz = aimag(zz)
        znupi_i = cmplx(-aimag(znupi),real(znupi),kp)
        abs_reznu = abs(re_znu)
! Out of the range. nregion=0
        IF (abs1(zz)>ai_arg_m) THEN
          nregion = 0
          RETURN
        END IF
! By series expansions. nregion=1,2
        IF (abs(zz)<1) THEN
          IF (abs_reznu>ihuge) THEN
            nregion = 0
            RETURN
          END IF
          IF (abs(znu-nint(re_znu))>0.3) THEN
            nregion = 1
            RETURN
          END IF
          nregion = 2
          RETURN
        END IF
        IF (abs_reznu>ai_arg_m) THEN
          nregion = 0
          RETURN
        END IF
        nregion = 3 ! By numerical integration.
      END FUNCTION num_region

      FUNCTION abs1(za) RESULT (ab)
        COMPLEX (kp), INTENT (IN) :: za
        REAL (kp) :: ab
! FUNCTION abs1 calculates a rough absolute value of complex argument za.
! The named constant kp is declared in MODULE mod_zbes.
! This is invoked by bessel1, bessel2, hankel1, hankel2, num_region,
! bes1_series, bes2_series, bes2_srs_init, def_bessel1, cylin_inte, cal_phi and
! zdgaus8d.
        ab = abs(real(za)) + abs(aimag(za))
      END FUNCTION abs1

    END MODULE mod_zbes

