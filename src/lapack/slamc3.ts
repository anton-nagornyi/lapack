import {Pointer} from "../pointer";

export const slamc3 = (a: Pointer<number>, b: Pointer<number>): number =>
{
    /*  -- LAPACK auxiliary routine (version 3.2) -- */
    /*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
    /*     November 2006 */

    /*     .. Scalar Arguments .. */
    /*     .. */

    /*  Purpose */
    /*  ======= */

    /*  SLAMC3  is intended to force  A  and  B  to be stored prior to doing */
    /*  the addition of  A  and  B ,  for use in situations where optimizers */
    /*  might hold one of these in a register. */

    /*  Arguments */
    /*  ========= */

    /*  A       (input) REAL */
    /*  B       (input) REAL */
    /*          The values A and B. */

    /* ===================================================================== */

    /*     .. Executable Statements .. */
    return Math.fround(Math.fround(a.value) + Math.fround(b.value));

    /*     End of SLAMC3 */

}; /* slamc3_ */