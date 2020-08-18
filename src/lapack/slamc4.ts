import {Pointer} from "../pointer";
import {slamc3} from "./slamc3";

export const slamc4 = (emin: Pointer<number>, start: Pointer<number>, base: Pointer<number>): number =>
{
    emin.value = Math.trunc(emin.value);
    base.value = Math.trunc(base.value);
    start.value = Math.fround(start.value);

    /* System generated locals */
    const r__1 = new Pointer(0);

    /*  -- LAPACK auxiliary routine (version 3.2) -- */
    /*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
    /*     November 2006 */

    /*     .. Scalar Arguments .. */
    /*     .. */

    /*  Purpose */
    /*  ======= */

    /*  SLAMC4 is a service routine for SLAMC2. */

    /*  Arguments */
    /*  ========= */

    /*  EMIN    (output) INTEGER */
    /*          The minimum exponent before (gradual) underflow, computed by */
    /*          setting A = START and dividing by BASE until the previous A */
    /*          can not be recovered. */

    /*  START   (input) REAL */
    /*          The starting point for determining EMIN. */

    /*  BASE    (input) INTEGER */
    /*          The base of the machine. */

    /* ===================================================================== */

    /*     .. Local Scalars .. */
    /*     .. */
    /*     .. External Functions .. */
    /*     .. */
    /*     .. Executable Statements .. */

    let a = start.value;
    const one = 1.;
    const rbase = Math.fround(one / base.value);
    const zero = new Pointer(0.);
    emin.value = 1;
    r__1.value = Math.fround(a * rbase);
    let b1 = slamc3(r__1, zero);
    let c1 = a;
    let c2 = a;
    let d1 = a;
    let d2 = a;
    /* +    WHILE( ( C1.EQ.A ).AND.( C2.EQ.A ).AND. */
    /*    $       ( D1.EQ.A ).AND.( D2.EQ.A )      )LOOP */
    while (c1 == a && c2 == a && d1 == a && d2 == a)
    {
        --emin.value;
        a = b1;
        r__1.value = Math.fround(a / base.value);
        b1 = slamc3(r__1, zero);
        r__1.value = Math.fround(b1 * base.value);
        c1 = slamc3(r__1, zero);
        d1 = zero.value;
        let i__1 = base.value;
        for (let i = 1; i <= i__1; ++i)
        {
            d1 += b1;
        }
        r__1.value = Math.fround(a * rbase);
        const b2 = slamc3(r__1, zero);
        r__1.value = Math.fround(b2 / rbase);
        c2 = slamc3(r__1, zero);
        d2 = zero.value;
        i__1 = Math.fround(base.value);
        for (let i = 1; i <= i__1; ++i)
        {
            d2 += b2;
            d2 = Math.fround(d2);
        }
    }
    /* +    END WHILE */

    return 0;

    /*     End of SLAMC4 */

}; /* slamc4_ */