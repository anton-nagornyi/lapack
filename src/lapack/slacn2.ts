import {Pointer} from "../pointer";
import {isamax} from "../blas/isamax";
import {scopy} from "../blas/scopy";
import {offset} from "../blas/offset";
import {r_sign} from "./r_sign";
import {i_nint} from "./i_nint";
import {sasum} from "../blas/sasum";

const c__1 = 1;
const c_b11 = 1.;
const a_offset = -1;
const {abs} = Math;

export const slacn2 = (n: number, v: Float32Array, x: Float32Array, isgn: Int32Array,
                       est: Pointer<number>, kase: Pointer<number>, isave: Int32Array): number =>
{

    let altsgn: number, estold: number;


    /*  -- LAPACK auxiliary routine (version 3.2) -- */
    /*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
    /*     November 2006 */

    /*     .. Scalar Arguments .. */
    /*     .. */
    /*     .. Array Arguments .. */
    /*     .. */

    /*  Purpose */
    /*  ======= */

    /*  SLACN2 estimates the 1-norm of a square, real matrix A. */
    /*  Reverse communication is used for evaluating matrix-vector products. */

    /*  Arguments */
    /*  ========= */

    /*  N      (input) INTEGER */
    /*         The order of the matrix.  N >= 1. */

    /*  V      (workspace) REAL array, dimension (N) */
    /*         On the final return, V = A*W,  where  EST = norm(V)/norm(W) */
    /*         (W is not returned). */

    /*  X      (input/output) REAL array, dimension (N) */
    /*         On an intermediate return, X should be overwritten by */
    /*               A * X,   if KASE=1, */
    /*               A' * X,  if KASE=2, */
    /*         and SLACN2 must be re-called with all the other parameters */
    /*         unchanged. */

    /*  ISGN   (workspace) INTEGER array, dimension (N) */

    /*  EST    (input/output) REAL */
    /*         On entry with KASE = 1 or 2 and ISAVE(1) = 3, EST should be */
    /*         unchanged from the previous call to SLACN2. */
    /*         On exit, EST is an estimate (a lower bound) for norm(A). */

    /*  KASE   (input/output) INTEGER */
    /*         On the initial call to SLACN2, KASE should be 0. */
    /*         On an intermediate return, KASE will be 1 or 2, indicating */
    /*         whether X should be overwritten by A * X  or A' * X. */
    /*         On the final return from SLACN2, KASE will again be 0. */

    /*  ISAVE  (input/output) INTEGER array, dimension (3) */
    /*         ISAVE is used to save variables between calls to SLACN2 */

    /*  Further Details */
    /*  ======= ======= */

    /*  Contributed by Nick Higham, University of Manchester. */
    /*  Originally named SONEST, dated March 16, 1988. */

    /*  Reference: N.J. Higham, "FORTRAN codes for estimating the one-norm of */
    /*  a real or complex matrix, with applications to condition estimation", */
    /*  ACM Trans. Math. Soft., vol. 14, no. 4, pp. 381-396, December 1988. */

    /*  This is a thread safe version of SLACON, which uses the array ISAVE */
    /*  in place of a SAVE statement, as follows: */

    /*     SLACON     SLACN2 */
    /*      JUMP     ISAVE(1) */
    /*      J        ISAVE(2) */
    /*      ITER     ISAVE(3) */

    /*  ===================================================================== */

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

    /* Function Body */
    if (kase.value == 0)
    {
        for (let i = 1; i <= n; ++i)
        {
            x[a_offset + i] = 1. / n;
            /* L10: */
        }
        kase.value = 1;
        isave[a_offset + 1] = 1;
        return 0;
    }

    switch (isave[a_offset + 1])
    {
        case 1:
            return L20(n, v, x, est, kase, isave, isgn);
        case 2:
            return L40(n, x, isave, kase);
        case 3:
            return L70(n, x, v, isave, isgn, est, kase);
        case 4:
            return L110(n, x, isave, kase);
        case 5:
            return L140(n, x, v, est, kase);
        default:
            return L70(n, x, v, isave, isgn, est, kase);
    }

    /*     End of SLACN2 */

}; /* slacn2_ */

const L20 = (n: number, v: Float32Array, x: Float32Array, est: Pointer<number>, kase: Pointer<number>, isave: Int32Array, isgn: Int32Array,): number =>
{
    /*     ................ ENTRY   (ISAVE( 1 ) = 1) */
    /*     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY A*X. */

    if (n == 1)
    {
        v[a_offset + 1] = x[a_offset + 1];
        est.value = abs(v[a_offset + 1]);
        /*        ... QUIT */
        return L150(kase);
    }
    est.value = sasum(n, offset(x, a_offset + 1), c__1);

    for (let i = 1; i <= n; ++i)
    {
        x[a_offset + i] = r_sign(c_b11, x[a_offset + i]);
        isgn[a_offset + i] = i_nint(x[a_offset + i]);
        /* L30: */
    }
    kase.value = 2;
    isave[a_offset + 1] = 2;
    return 0;
};

const L40 = (n: number, x: Float32Array, isave: Int32Array, kase: Pointer<number>): number =>
{
    /*     ................ ENTRY   (ISAVE( 1 ) = 2) */
    /*     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X. */
    isave[a_offset + 2] = isamax(n, offset(x, a_offset + 1), c__1);
    isave[a_offset + 3] = 2;
    return L50(n, x, isave, kase);
};
const L50 = (n: number, x: Float32Array, isave: Int32Array, kase: Pointer<number>): number =>
{
    /*     MAIN LOOP - ITERATIONS 2,3,...,ITMAX. */

    for (let i = 1; i <= n; ++i)
    {
        x[a_offset + i] = 0.;
        /* L60: */
    }
    x[a_offset + isave[a_offset + 2]] = 1.;
    kase.value = 1;
    isave[a_offset + 1] = 3;
    return 0;
};
const L70 = (n: number, x: Float32Array, v: Float32Array, isave: Int32Array, isgn: Int32Array, est: Pointer<number>, kase: Pointer<number>): number =>
{
    /*     ................ ENTRY   (ISAVE( 1 ) = 3) */
    /*     X HAS BEEN OVERWRITTEN BY A*X. */

    scopy(n, offset(x, a_offset + 1), c__1, offset(v, a_offset + 1), c__1);
    const estold = est.value;
    est.value = sasum(n, offset(v, a_offset + 1), c__1);

    for (let i = 1; i <= n; ++i)
    {
        const r__1 = r_sign(c_b11, x[a_offset + i]);
        if (i_nint(r__1) != isgn[a_offset + i])
        {
            return L90(n, x, isave, isgn, estold, est, kase);
        }
        /* L80: */
    }
    /*     REPEATED SIGN VECTOR DETECTED, HENCE ALGORITHM HAS CONVERGED. */
    return L120(n, x, isave, kase);
};
const L90 = (n: number, x: Float32Array, isave: Int32Array, isgn: Int32Array, estold: number, est: Pointer<number>, kase: Pointer<number>): number =>
{
    /*     TEST FOR CYCLING. */
    if (est.value <= estold)
    {
        return L120(n, x, isave, kase);
    }
    for (let i = 1; i <= n; ++i)
    {
        x[a_offset + i] = r_sign(c_b11, x[a_offset + i]);
        isgn[a_offset + i] = i_nint(x[a_offset + i]);
        /* L100: */
    }
    kase.value = 2;
    isave[a_offset + 1] = 4;
    return 0;
};
const L110 = (n: number, x: Float32Array, isave: Int32Array, kase: Pointer<number>): number =>
{
    /*     ................ ENTRY   (ISAVE( 1 ) = 4) */
    /*     X HAS BEEN OVERWRITTEN BY TRANSPOSE(A)*X. */

    const jlast = isave[a_offset + 2];
    let r__1: number;
    isave[a_offset + 2] = isamax(n, offset(x, a_offset + 1), c__1);
    if (x[a_offset + jlast] != (r__1 = x[a_offset + isave[a_offset + 2]], abs(r__1)) && isave[a_offset + 3] < 5)
    {
        ++isave[a_offset + 3];
        return L50(n, x, isave, kase);
    }
    return L120(n, x, isave, kase);
};
const L120 = (n: number, x: Float32Array, isave: Int32Array, kase: Pointer<number>): number =>
{
    /*     ITERATION COMPLETE.  FINAL STAGE. */

    let altsgn = 1.;
    for (let i = 1; i <= n; ++i)
    {
        x[a_offset + i] = altsgn * ((i - 1) / (n - 1) + 1.);
        altsgn = -altsgn;
        /* L130: */
    }
    kase.value = 1;
    isave[a_offset + 1] = 5;
    return 0;
};
const L140 = (n: number, x: Float32Array, v: Float32Array, est: Pointer<number>, kase: Pointer<number>): number =>
{
    /*     ................ ENTRY   (ISAVE( 1 ) = 5) */
    /*     X HAS BEEN OVERWRITTEN BY A*X. */

    const temp = sasum(n, offset(x, a_offset + 1), c__1) / (n * 3) * 2.;
    if (temp > est.value)
    {
        scopy(n, offset(x, a_offset + 1), c__1, offset(v, a_offset + 1), c__1);
        est.value = temp;
    }
    return L150(kase);
};
const L150 = (kase: Pointer<number>): number =>
{
    kase.value = 0;
    return 0;
};