import {slamch} from "./slamch";
import {sscal} from "../blas/sscal";
import {slabad} from "./slabad";
import {Pointer} from "../pointer";

const {abs} = Math;

export const srscl = (n: number, sa: number, sx: Float32Array, incx: number): number =>
{
    let mul: number, cden: number;
    let done = false;
    let cnum: number, cden1: number, cnum1: number;
    let bignum = new Pointer(0);
    let smlnum = new Pointer(0);


    /*  -- LAPACK auxiliary routine (version 3.2) -- */
    /*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
    /*     November 2006 */

    /*     .. Scalar Arguments .. */
    /*     .. */
    /*     .. Array Arguments .. */
    /*     .. */

    /*  Purpose */
    /*  ======= */

    /*  SRSCL multiplies an n-element real vector x by the real scalar 1/a. */
    /*  This is done without overflow or underflow as long as */
    /*  the final result x/a does not overflow or underflow. */

    /*  Arguments */
    /*  ========= */

    /*  N       (input) INTEGER */
    /*          The number of components of the vector x. */

    /*  SA      (input) REAL */
    /*          The scalar a which is used to divide each component of x. */
    /*          SA must be >= 0, or the subroutine will divide by zero. */

    /*  SX      (input/output) REAL array, dimension */
    /*                         (1+(N-1)*abs(INCX)) */
    /*          The n-element vector x. */

    /*  INCX    (input) INTEGER */
    /*          The increment between successive values of the vector SX. */
    /*          > 0:  SX(1) = X(1) and SX(1+(i-1)*INCX) = x(i),     1< i<= n */

    /* ===================================================================== */

    /*     .. Parameters .. */
    /*     .. */
    /*     .. Local Scalars .. */
    /*     .. */
    /*     .. External Functions .. */
    /*     .. */
    /*     .. External Subroutines .. */
    /*     .. */
    /*     .. Intrinsic Functions .. */
    /*     .. */
    /*     .. Executable Statements .. */

    /*     Quick return if possible */

    /* Function Body */
    if (n <= 0)
    {
        return 0;
    }

    /*     Get machine parameters */

    smlnum.value = slamch("S");
    bignum.value = 1. / smlnum.value;
    slabad(smlnum, bignum);

    /*     Initialize the denominator to SA and the numerator to 1. */

    cden = sa;
    cnum = 1.;

    while (!done)
    {
        cden1 = cden * smlnum.value;
        cnum1 = cnum / bignum.value;
        if (abs(cden1) > abs(cnum) && cnum != 0.)
        {

            /*        Pre-multiply X by SMLNUM if CDEN is large compared to CNUM. */

            mul = smlnum.value;
            done = false;
            cden = cden1;
        } else if (abs(cnum1) > abs(cden))
        {

            /*        Pre-multiply X by BIGNUM if CDEN is small compared to CNUM. */

            mul = bignum.value;
            done = false;
            cnum = cnum1;
        } else
        {

            /*        Multiply X by CNUM / CDEN and return. */

            mul = cnum / cden;
            done = true;
        }

        /*     Scale the vector X by MUL */

        sscal(n, mul, sx, incx);
    }
    return 0;

    /*     End of SRSCL */

}; /* srscl_ */