import {Pointer} from "../pointer";
import {slamc3} from "./slamc3";

const c_b32 = new Pointer(0.);

export const slamc5 = (beta: Pointer<number>, p: Pointer<number>, emin: Pointer<number>,
                       ieee: Pointer<boolean>, emax: Pointer<number>, rmax: Pointer<number>) =>
{
    beta.value = Math.trunc(beta.value);
    p.value = Math.trunc(p.value);
    emin.value = Math.trunc(emin.value);
    emax.value = Math.trunc(emax.value);
    /*  -- LAPACK auxiliary routine (version 3.2) -- */
    /*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
    /*     November 2006 */

    /*     .. Scalar Arguments .. */
    /*     .. */

    /*  Purpose */
    /*  ======= */

    /*  SLAMC5 attempts to compute RMAX, the largest machine floating-point */
    /*  number, without overflow.  It assumes that EMAX + abs(EMIN) sum */
    /*  approximately to a power of 2.  It will fail on machines where this */
    /*  assumption does not hold, for example, the Cyber 205 (EMIN = -28625, */
    /*  EMAX = 28718).  It will also fail if the value supplied for EMIN is */
    /*  too large (i.e. too close to zero), probably with overflow. */

    /*  Arguments */
    /*  ========= */

    /*  BETA    (input) INTEGER */
    /*          The base of floating-point arithmetic. */

    /*  P       (input) INTEGER */
    /*          The number of base BETA digits in the mantissa of a */
    /*          floating-point value. */

    /*  EMIN    (input) INTEGER */
    /*          The minimum exponent before (gradual) underflow. */

    /*  IEEE    (input) LOGICAL */
    /*          A logical flag specifying whether or not the arithmetic */
    /*          system is thought to comply with the IEEE standard. */

    /*  EMAX    (output) INTEGER */
    /*          The largest exponent before overflow */

    /*  RMAX    (output) REAL */
    /*          The largest machine floating-point number. */

    /* ===================================================================== */

    /*     .. Parameters .. */
    /*     .. */
    /*     .. Local Scalars .. */
    /*     .. */
    /*     .. External Functions .. */
    /*     .. */
    /*     .. Intrinsic Functions .. */
    /*     .. */
    /*     .. Executable Statements .. */

    /*     First compute LEXP and UEXP, two powers of 2 that bound */
    /*     abs(EMIN). We then assume that EMAX + abs(EMIN) will sum */
    /*     approximately to the bound that is closest to abs(EMIN). */
    /*     (EMAX is the exponent of the required number RMAX). */

    let lexp = 1;
    let exbits = 1;
    let try__;
    for (; ;)
    {
        try__ = lexp << 1;
        if (try__ <= -emin.value)
        {
            lexp = try__;
            ++exbits;
            continue
        }
        break;
    }

    let uexp;
    if (lexp == -emin.value)
    {
        uexp = lexp;
    } else
    {
        uexp = try__;
        ++exbits;
    }

    /*     Now -LEXP is less than or equal to EMIN, and -UEXP is greater */
    /*     than or equal to EMIN. EXBITS is the number of bits needed to */
    /*     store the exponent. */

    const expsum = (uexp + emin.value > -lexp - emin.value) ? lexp << 1 : uexp << 1;

    /*     EXPSUM is the exponent range, approximately equal to */
    /*     EMAX - EMIN + 1 . */

    emax.value = expsum + emin.value - 1;
    const nbits = exbits + 1 + p.value;

    /*     NBITS is the total number of bits needed to store a */
    /*     floating-point number. */

    if (nbits % 2 == 1 && beta.value == 2)
    {

        /*        Either there are an odd number of bits used to store a */
        /*        floating-point number, which is unlikely, or some bits are */
        /*        not used in the representation of numbers, which is possible, */
        /*        (e.g. Cray machines) or the mantissa has an implicit bit, */
        /*        (e.g. IEEE machines, Dec Vax machines), which is perhaps the */
        /*        most likely. We have to assume the last alternative. */
        /*        If this is true, then we need to reduce EMAX by one because */
        /*        there must be some way of representing zero in an implicit-bit */
        /*        system. On machines like Cray, we are reducing EMAX by one */
        /*        unnecessarily. */

        --emax.value;
    }

    if (ieee.value)
    {

        /*        Assume we are on an IEEE machine which reserves one exponent */
        /*        for infinity and NaN. */

        --emax.value;
    }

    /*     Now create RMAX, the largest machine number, which should */
    /*     be equal to (1.0 - BETA**(-P)) * BETA**EMAX . */

    /*     First compute 1.0 - BETA**(-P), being careful that the */
    /*     result is less than 1.0 . */

    const recbas = Math.fround(1. / beta.value);
    const z__ = new Pointer(Math.fround(beta.value - 1.));
    const y = new Pointer(0.);
    let i__1 = p.value;
    let oldy = 0.;
    for (let i = 1; i <= i__1; ++i)
    {
        z__.value *= recbas;
        z__.value = Math.fround(z__.value);
        if (y.value < 1.)
        {
            oldy = y.value;
        }
        y.value = slamc3(y, z__);
        /* L20: */
    }
    if (y.value >= 1.)
    {
        y.value = oldy;
    }

    /*     Now multiply by BETA**EMAX to get RMAX. */

    i__1 = emax.value;
    const r__1 = new Pointer(0);
    for (let i = 1; i <= i__1; ++i)
    {
        r__1.value = Math.fround(y.value * beta.value);
        y.value = slamc3(r__1, c_b32);
        /* L30: */
    }

    rmax.assign(y);
    return 0;

    /*     End of SLAMC5 */

}; /* slamc5_ */