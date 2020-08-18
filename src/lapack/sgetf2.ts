import {LapackInfo} from "../info";
import {slamch} from "./slamch";
import {isamax} from "../blas/isamax";
import {sswap} from "../blas/sswap";
import {sscal} from "../blas/sscal";
import {sger} from "../blas/sger";
import {offset} from "../blas/offset";

const c__1 = 1;
const c_b8 = -1.;

const {max, min, abs} = Math;

export const sgetf2 = (m: number, n: number, a: Float32Array, lda: number,
                       ipiv: Int32Array, info: LapackInfo): number =>
{
    /* System generated locals */
    let a_dim1: number, a_offset: number, i__1: number, i__2: number, i__3: number;
    let r__1: number;

    /* Local variables */
    let i__: number, j: number, jp: number;
    let sfmin: number;

    /*  -- LAPACK routine (version 3.2) -- */
    /*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
    /*     November 2006 */

    /*     .. Scalar Arguments .. */
    /*     .. */
    /*     .. Array Arguments .. */
    /*     .. */

    /*  Purpose */
    /*  ======= */

    /*  SGETF2 computes an LU factorization of a general m-by-n matrix A */
    /*  using partial pivoting with row interchanges. */

    /*  The factorization has the form */
    /*     A = P * L * U */
    /*  where P is a permutation matrix, L is lower triangular with unit */
    /*  diagonal elements (lower trapezoidal if m > n), and U is upper */
    /*  triangular (upper trapezoidal if m < n). */

    /*  This is the right-looking Level 2 BLAS version of the algorithm. */

    /*  Arguments */
    /*  ========= */

    /*  M       (input) INTEGER */
    /*          The number of rows of the matrix A.  M >= 0. */

    /*  N       (input) INTEGER */
    /*          The number of columns of the matrix A.  N >= 0. */

    /*  A       (input/output) REAL array, dimension (LDA,N) */
    /*          On entry, the m by n matrix to be factored. */
    /*          On exit, the factors L and U from the factorization */
    /*          A = P*L*U; the unit diagonal elements of L are not stored. */

    /*  LDA     (input) INTEGER */
    /*          The leading dimension of the array A.  LDA >= max(1,M). */

    /*  IPIV    (output) INTEGER array, dimension (min(M,N)) */
    /*          The pivot indices; for 1 <= i <= min(M,N), row i of the */
    /*          matrix was interchanged with row IPIV(i). */

    /*  INFO    (output) INTEGER */
    /*          = 0: successful exit */
    /*          < 0: if INFO = -k, the k-th argument had an illegal value */
    /*          > 0: if INFO = k, U(k,k) is exactly zero. The factorization */
    /*               has been completed, but the factor U is exactly */
    /*               singular, and division by zero will occur if it is used */
    /*               to solve a system of equations. */

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

    /*     Test the input parameters. */

    /* Parameter adjustments */
    a_dim1 = lda;
    a_offset = -(1 + a_dim1);
    const ipiv_offset = -1;

    /* Function Body */
    info.value = 0;
    if (m < 0)
    {
        info.value = -1;
    } else if (n < 0)
    {
        info.value = -2;
    } else if (lda < max(1, m))
    {
        info.value = -4;
    }
    if (info.value != 0)
    {
        i__1 = -(info.value);
        console.error(`SGETF2 parameter ${i__1}`,);
        return 0;
    }

    /*     Quick return if possible */

    if (m == 0 || n == 0)
    {
        return 0;
    }

    /*     Compute machine safe minimum */

    sfmin = slamch("S");

    i__1 = min(m, n);
    for (j = 1; j <= i__1; ++j)
    {

        /*        Find pivot and test for singularity. */

        i__2 = m - j + 1;
        jp = j - 1 + isamax(i__2, offset(a, a_offset + j + j * a_dim1), c__1);
        ipiv[ipiv_offset + j] = jp;
        if (a[a_offset + jp + j * a_dim1] != 0.)
        {

            /*           Apply the interchange to columns 1:N. */

            if (jp != j)
            {
                sswap(n, offset(a, a_offset + j + a_dim1), lda, offset(a, a_offset + jp + a_dim1), lda);
            }

            /*           Compute elements J+1:M of J-th column. */

            if (j < m)
            {
                if ((r__1 = a[a_offset + j + j * a_dim1], abs(r__1)) >= sfmin)
                {
                    i__2 = m - j;
                    r__1 = 1. / a[a_offset + j + j * a_dim1];
                    sscal(i__2, r__1, offset(a, a_offset + j + 1 + j * a_dim1), c__1);
                } else
                {
                    i__2 = m - j;
                    for (i__ = 1; i__ <= i__2; ++i__)
                    {
                        a[a_offset + j + i__ + j * a_dim1] /= a[a_offset + j + j * a_dim1];
                        /* L20: */
                    }
                }
            }

        } else if (info.value == 0)
        {

            info.value = j;
        }

        if (j < min(m, n))
        {

            /*           Update trailing submatrix. */

            i__2 = m - j;
            i__3 = n - j;
            sger(i__2, i__3, c_b8, offset(a, a_offset + j + 1 + j * a_dim1), c__1, offset(a, a_offset + j + (
                j + 1) * a_dim1), lda, offset(a, a_offset + j + 1 + (j + 1) * a_dim1), lda);
        }
        /* L10: */
    }
    return 0;

    /*     End of SGETF2 */

}; /* sgetf2_ */