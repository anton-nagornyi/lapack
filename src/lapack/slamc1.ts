import {Pointer} from "../pointer";
import {slamc3} from "./slamc3";

let first = true;
const lt = new Pointer(0);
const lrnd = new Pointer(false);
const lbeta = new Pointer(0);
const lieee1 = new Pointer(false);

export const slamc1 = (beta: Pointer<number>, t: Pointer<number>, rnd: Pointer<boolean>, ieee1: Pointer<boolean>): number =>
{
    /* Initialized data */

    /* System generated locals */
    const r__1 = new Pointer(0);
    const r__2 = new Pointer(0);

    /*  -- LAPACK auxiliary routine (version 3.2) -- */
    /*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
    /*     November 2006 */

    /*     .. Scalar Arguments .. */
    /*     .. */

    /*  Purpose */
    /*  ======= */

    /*  SLAMC1 determines the machine parameters given by BETA, T, RND, and */
    /*  IEEE1. */

    /*  Arguments */
    /*  ========= */

    /*  BETA    (output) INTEGER */
    /*          The base of the machine. */

    /*  T       (output) INTEGER */
    /*          The number of ( BETA ) digits in the mantissa. */

    /*  RND     (output) LOGICAL */
    /*          Specifies whether proper rounding  ( RND = .TRUE. )  or */
    /*          chopping  ( RND = .FALSE. )  occurs in addition. This may not */
    /*          be a reliable guide to the way in which the machine performs */
    /*          its arithmetic. */

    /*  IEEE1   (output) LOGICAL */
    /*          Specifies whether rounding appears to be done in the IEEE */
    /*          'round to nearest' style. */

    /*  Further Details */
    /*  =============== */

    /*  The routine is based on the routine  ENVRON  by Malcolm and */
    /*  incorporates suggestions by Gentleman and Marovich. See */

    /*     Malcolm M. A. (1972) Algorithms to reveal properties of */
    /*        floating-point arithmetic. Comms. of the ACM, 15, 949-951. */

    /*     Gentleman W. M. and Marovich S. B. (1974) More on algorithms */
    /*        that reveal properties of floating point arithmetic units. */
    /*        Comms. of the ACM, 17, 276-277. */

    /* ===================================================================== */

    /*     .. Local Scalars .. */
    /*     .. */
    /*     .. External Functions .. */
    /*     .. */
    /*     .. Save statement .. */
    /*     .. */
    /*     .. Data statements .. */
    /*     .. */
    /*     .. Executable Statements .. */

    if (first) {
        const one = new Pointer(1.);

        /*        LBETA,  LIEEE1,  LT and  LRND  are the  local values  of  BETA, */
        /*        IEEE1, T and RND. */

        /*        Throughout this routine  we use the function  SLAMC3  to ensure */
        /*        that relevant values are  stored and not held in registers,  or */
        /*        are not affected by optimizers. */

        /*        Compute  a = 2.0**m  with the  smallest positive integer m such */
        /*        that */

        /*           fl( a + 1.0 ) = a. */

        const a = new Pointer(1.);
        const c__ = new Pointer(1.);
        /* +       WHILE( C.EQ.ONE )LOOP */
            while (c__.value == one.value)
            {
                a.value *= 2;
                a.value = Math.fround(a.value);
                c__.value = slamc3(a, one);
                r__1.value = Math.fround(-a.value);
                c__.value = slamc3(c__, r__1);
            }

        /* +       END WHILE */

        /*        Now compute  b = 2.0**m  with the smallest positive integer m */
        /*        such that */

        /*           fl( a + b ) .gt. a. */

        const b = new Pointer(1.);
        c__.value = slamc3(a, b);

        /* +       WHILE( C.EQ.A )LOOP */
            while (c__.value == a.value)
            {
                b.value *= 2;
                b.value = Math.fround(b.value);
                c__.value = slamc3(a, b);
            }
        /* +       END WHILE */

        /*        Now compute the base.  a and c  are neighbouring floating point */
        /*        numbers  in the  interval  ( beta**t, beta**( t + 1 ) )  and so */
        /*        their difference is beta. Adding 0.25 to c is to ensure that it */
        /*        is truncated to beta and not ( beta - 1 ). */

        const qtr = one.value / 4;
        const savec = new Pointer(c__.value);
        r__1.value = Math.fround(-a.value);
        c__.value = slamc3(c__, r__1);
        lbeta.value = Math.trunc(c__.value + qtr);

        /*        Now determine whether rounding or chopping occurs,  by adding a */
        /*        bit  less  than  beta/2  and a  bit  more  than  beta/2  to  a. */

        b.assign(lbeta);
        r__1.value = b.value / 2;
        r__2.value = -b.value / 100;
        const f = new Pointer(slamc3(r__1, r__2));
        c__.value = slamc3(f, a);
        if (c__.value == a.value)
        {
            lrnd.value = true;
        } else {
            lrnd.value = false;
        }
        r__1.value = b.value / 2;
        r__2.value = b.value / 100;
        f.value = slamc3(r__1, r__2);
        c__.value = slamc3(f, a);
        if (lrnd.value && c__.value == a.value)
        {
            lrnd.value = false;
        }

        /*        Try and decide whether rounding is done in the  IEEE  'round to */
        /*        nearest' style. B/2 is half a unit in the last place of the two */
        /*        numbers A and SAVEC. Furthermore, A is even, i.e. has last  bit */
        /*        zero, and SAVEC is odd. Thus adding B/2 to A should not  change */
        /*        A, but adding B/2 to SAVEC should change SAVEC. */

        r__1.value = b.value / 2;
        const t1 = slamc3(r__1, a);
        r__1.value = b.value / 2;
        const t2 = slamc3(r__1, savec);
        lieee1.value = t1 == a.value && t2 > savec.value && lrnd.value;

        /*        Now find  the  mantissa, t.  It should  be the  integer part of */
        /*        log to the base beta of a,  however it is safer to determine  t */
        /*        by powering.  So we find t as the smallest positive integer for */
        /*        which */

        /*           fl( beta**t + 1.0 ) = 1.0. */

        lt.value = 0;
        a.value = 1.;
        c__.value = 1.;

        /* +       WHILE( C.EQ.ONE )LOOP */
            while (c__.value == one.value) {
                ++lt.value;
                a.value *= lbeta.value;
                a.value = Math.fround(a.value);
                c__.value = slamc3(a, one);
                r__1.value = Math.fround(-a.value);
                c__.value = slamc3(c__, r__1);
            }
        /* +       END WHILE */

    }


beta.value = Math.trunc(lbeta.value);
t.assign(lt);
rnd.assign(lrnd);
ieee1.assign(lieee1);
first = false;

    return 0;

    /*     End of SLAMC1 */

}; /* slamc1_ */