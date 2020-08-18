import {Pointer} from "../pointer";
import {slamc3} from "./slamc3";
import {slamc4} from "./slamc4";
import {slamc5} from "./slamc5";
import {pow_ri} from "./pow_ri";
import {slamc1} from "./slamc1";

let first = true;
let iwarn = false;
const fmt_9999 = (emin: number) => `WARNING. The value EMIN may be incorrect: ${emin}
EMIN = ${emin}
If, after inspection, the value EMIN looks {v1} acceptable please comment out {v1} the IF block as marked within the code of 
routine {v1} SLAMC2 {v2} otherwise supply EMIN explicitly.`;

const lt = new Pointer(0);
const leps = new Pointer(0);
const lbeta = new Pointer(0);
const lemin = new Pointer(0);
const lemax = new Pointer(0);
const lrmin = new Pointer(0);
const lrmax = new Pointer(0);

const {abs, min, max} = Math;

export const slamc2 = (beta: Pointer<number>, t: Pointer<number>, rnd: Pointer<boolean>, eps: Pointer<number>,
                       emin: Pointer<number>, rmin: Pointer<number>, emax: Pointer<number>, rmax: Pointer<number>) =>
{
    /* Initialized data */

    /* Format strings */


    /* System generated locals */
    const i__1 = new Pointer(0);
    const r__1 = new Pointer(0);
    const r__2 = new Pointer(0);
    const r__3 = new Pointer(0);
    const r__4 = new Pointer(0);
    const r__5 = new Pointer(0);

    /* Local variables */
    const a = new Pointer(0);
    const b = new Pointer(0);
    const c__ = new Pointer(0);

    let ieee = new Pointer(false);
    const half = new Pointer(0);
    const lrnd = new Pointer(false);

    const gnmin = new Pointer(0);
    const small = new Pointer(0);
    const gpmin = new Pointer(0);
    const third = new Pointer(0);

    const sixth = new Pointer(0);
    const lieee1 = new Pointer(false);
    const ngnmin = new Pointer(0);
    const ngpmin = new Pointer(0);

    /* Fortran I/O blocks */
// static cilist io___58 = { 0, 6, 0, fmt_9999, 0 };


    /*  -- LAPACK auxiliary routine (version 3.2) -- */
    /*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
    /*     November 2006 */

    /*     .. Scalar Arguments .. */
    /*     .. */

    /*  Purpose */
    /*  ======= */

    /*  SLAMC2 determines the machine parameters specified in its argument */
    /*  list. */

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

    /*  EPS     (output) REAL */
    /*          The smallest positive number such that */

    /*             fl( 1.0 - EPS ) .LT. 1.0, */

    /*          where fl denotes the computed value. */

    /*  EMIN    (output) INTEGER */
    /*          The minimum exponent before (gradual) underflow occurs. */

    /*  RMIN    (output) REAL */
    /*          The smallest normalized number for the machine, given by */
    /*          BASE**( EMIN - 1 ), where  BASE  is the floating point value */
    /*          of BETA. */

    /*  EMAX    (output) INTEGER */
    /*          The maximum exponent before overflow occurs. */

    /*  RMAX    (output) REAL */
    /*          The largest positive number for the machine, given by */
    /*          BASE**EMAX * ( 1 - EPS ), where  BASE  is the floating point */
    /*          value of BETA. */

    /*  Further Details */
    /*  =============== */

    /*  The computation of  EPS  is based on a routine PARANOIA by */
    /*  W. Kahan of the University of California at Berkeley. */

    /* ===================================================================== */

    /*     .. Local Scalars .. */
    /*     .. */
    /*     .. External Functions .. */
    /*     .. */
    /*     .. External Subroutines .. */
    /*     .. */
    /*     .. Intrinsic Functions .. */
    /*     .. */
    /*     .. Save statement .. */
    /*     .. */
    /*     .. Data statements .. */
    /*     .. */
    /*     .. Executable Statements .. */

    if (first)
    {
        const zero = new Pointer(0.);
        const one = new Pointer(1.);
        const two = 2.;

        /*        LBETA, LT, LRND, LEPS, LEMIN and LRMIN  are the local values of */
        /*        BETA, T, RND, EPS, EMIN and RMIN. */

        /*        Throughout this routine  we use the function  SLAMC3  to ensure */
        /*        that relevant values are stored  and not held in registers,  or */
        /*        are not affected by optimizers. */

        /*        SLAMC1 returns the parameters  LBETA, LT, LRND and LIEEE1. */

        lbeta.value = Math.trunc(lbeta.value);
        lt.value = Math.trunc(lt.value);
        slamc1(lbeta, lt, lrnd, lieee1);

        /*        Start to find EPS. */

        b.assign(lbeta);
        i__1.value = -lt.value;
        a.value = pow_ri(b, i__1);
        leps.assign(a);

        /*        Try some tricks to see whether or not this is the correct  EPS. */

        b.value = two / 3;
        half.value = one.value / 2;
        r__1.value = -half.value;
        sixth.value = slamc3(b, r__1);
        third.value = slamc3(sixth, sixth);
        r__1.value = -half.value;
        b.value = slamc3(third, r__1);
        b.value = slamc3(b, sixth);
        b.value = abs(b.value);

        if (b.value < leps.value)
        {
            b.assign(leps);
        }

        leps.value = 1.;
        /* +       WHILE( ( LEPS.GT.B ).AND.( B.GT.ZERO ) )LOOP */
        while (leps.value > b.value && b.value > zero.value)
        {
            leps.assign(b);
            r__1.value = half.value * leps.value;
            /* Computing 5th power */
            r__3.value = two;
            r__4.assign(r__3);
            r__3.value *= r__3.value;
            /* Computing 2nd power */
            r__5.assign(leps);
            r__2.value = r__4.value * (r__3.value * r__3.value) * (r__5.value * r__5.value);
            c__.value = slamc3(r__1, r__2);
            r__1.value = -c__.value;
            c__.value = slamc3(half, r__1);
            b.value = slamc3(half, c__);
            r__1.value = -b.value;
            c__.value = slamc3(half, r__1);
            b.value = slamc3(half, c__);
        }
        /* +       END WHILE */

        if (a.value < leps.value)
        {
            leps.assign(a);
        }

        /*        Computation of EPS complete. */

        /*        Now find  EMIN.  Let A = + or - 1, and + or - (1 + BASE**(-3)). */
        /*        Keep dividing  A by BETA until (gradual) underflow occurs. This */
        /*        is detected when we cannot recover the previous A. */

        const rbase = Math.fround(one.value / lbeta.value);
        small.assign(one);
        for (let i = 1; i <= 3; ++i)
        {
            r__1.value = Math.fround(small.value * rbase);
            small.value = slamc3(r__1, zero);
        }
        a.value = slamc3(one, small);
        slamc4(ngpmin, one, lbeta);
        r__1.value = -one.value;
        slamc4(ngnmin, r__1, lbeta);
        slamc4(gpmin, a, lbeta);
        r__1.value = -a.value;
        slamc4(gnmin, r__1, lbeta);
        ieee.value = false;

        if (ngpmin.value == ngnmin.value && gpmin.value == gnmin.value)
        {
            if (ngpmin.value == gpmin.value)
            {
                lemin.assign(ngpmin);
                /*            ( Non twos-complement machines, no gradual underflow; */
                /*              e.g.,  VAX ) */
            } else if (gpmin.value - ngpmin.value == 3)
            {
                lemin.value = ngpmin.value - 1 + lt.value;
                ieee.value = true;
                /*            ( Non twos-complement machines, with gradual underflow; */
                /*              e.g., IEEE standard followers ) */
            } else
            {
                lemin.value = min(ngpmin.value, gpmin.value);
                /*            ( A guess; no known machine ) */
                iwarn = true;
            }

        } else if (ngpmin.value == gpmin.value && ngnmin.value == gnmin.value)
        {
            if ((i__1.value = ngpmin.value - ngnmin.value, abs(i__1.value)) == 1)
            {
                lemin.value = max(ngpmin.value, ngnmin.value);
                /*            ( Twos-complement machines, no gradual underflow; */
                /*              e.g., CYBER 205 ) */
            } else
            {
                lemin.value = min(ngpmin.value, ngnmin.value);
                /*            ( A guess; no known machine ) */
                iwarn = true;
            }

        } else if ((i__1.value = ngpmin.value - ngnmin.value, abs(i__1.value)) == 1 && gpmin.value == gnmin.value)
        {
            if (gpmin.value - min(ngpmin.value, ngnmin.value) == 3)
            {
                lemin.value = max(ngpmin.value, ngnmin.value) - 1 + lt.value;
                /*            ( Twos-complement machines with gradual underflow; */
                /*              no known machine ) */
            } else
            {
                lemin.value = min(ngpmin.value, ngnmin.value);
                /*            ( A guess; no known machine ) */
                iwarn = true;
            }

        } else
        {
            /* Computing MIN */
            i__1.value = min(ngpmin.value, ngnmin.value);
            i__1.value = min(i__1.value, gpmin.value);
            lemin.value = min(i__1.value, gnmin.value);
            /*         ( A guess; no known machine ) */
            iwarn = true;
        }
        first = false;
        /* ** */
        /* Comment out this if block if EMIN is ok */
        if (iwarn)
        {
            console.warn(fmt_9999(lemin.value));
            first = true;
        }
        /* ** */

        /*        Assume IEEE arithmetic if we found denormalised  numbers above, */
        /*        or if arithmetic seems to round in the  IEEE style,  determined */
        /*        in routine SLAMC1. A true IEEE machine should have both  things */
        /*        true; however, faulty machines may have one or the other. */

        ieee = ieee || lieee1.value;

        /*        Compute  RMIN by successive division by  BETA. We could compute */
        /*        RMIN as BASE**( EMIN - 1 ),  but some machines underflow during */
        /*        this computation. */

        lrmin.value = 1.;
        i__1.value = 1 - lemin.value;
        for (let i = 1; i <= i__1.value; ++i)
        {
            r__1.value = Math.fround(lrmin.value * rbase);
            lrmin.value = slamc3(r__1, zero);
            /* L30: */
        }

        /*        Finally, call SLAMC5 to compute EMAX and RMAX. */

        slamc5(lbeta, lt, lemin, ieee, lemax, lrmax);
    }

    beta.assign(lbeta);
    t.assign(lt);
    rnd.assign(lrnd);
    eps.assign(leps);
    emin.assign(lemin);
    rmin.assign(lrmin);
    emax.assign(lemax);
    rmax.assign(lrmax);

    return 0;


    /*     End of SLAMC2 */

}; /* slamc2_ */