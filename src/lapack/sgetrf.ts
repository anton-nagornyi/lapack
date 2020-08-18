import {LapackInfo} from "../info";
import {ilaenv} from "./ilaenv";
import {slaswp} from "../blas/slaswp";
import {strsm} from "../blas/strsm";
import {sgemm} from "../blas/sgemm";
import {sgetf2} from "./sgetf2";
import {offset} from "../blas/offset";

const {max, min} = Math;
const c__1 = 1;
const c_n1 = -1;
const c_b16 = 1.;
const c_b19 = -1.;

export const sgetrf = (m: number, n: number, a: Float32Array, lda: number,
                       ipiv: Int32Array, info: LapackInfo): number =>
{
    /* System generated locals */
    let a_dim1: number, a_offset: number, i__1: number, i__2: number, i__3: number, i__4: number, i__5: number;

    /* Local variables */
    let i__: number, j: number, jb: number, nb: number;
    let iinfo = new LapackInfo();

    /*  -- LAPACK routine (version 3.2) -- */
    /*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
    /*     November 2006 */

    /*     .. Scalar Arguments .. */
    /*     .. */
    /*     .. Array Arguments .. */
    /*     .. */

    /*  Purpose */
    /*  ======= */

    /*  SGETRF computes an LU factorization of a general M-by-N matrix A */
    /*  using partial pivoting with row interchanges. */

    /*  The factorization has the form */
    /*     A = P * L * U */
    /*  where P is a permutation matrix, L is lower triangular with unit */
    /*  diagonal elements (lower trapezoidal if m > n), and U is upper */
    /*  triangular (upper trapezoidal if m < n). */

    /*  This is the right-looking Level 3 BLAS version of the algorithm. */

    /*  Arguments */
    /*  ========= */

    /*  M       (input) INTEGER */
    /*          The number of rows of the matrix A.  M >= 0. */

    /*  N       (input) INTEGER */
    /*          The number of columns of the matrix A.  N >= 0. */

    /*  A       (input/output) REAL array, dimension (LDA,N) */
    /*          On entry, the M-by-N matrix to be factored. */
    /*          On exit, the factors L and U from the factorization */
    /*          A = P*L*U; the unit diagonal elements of L are not stored. */

    /*  LDA     (input) INTEGER */
    /*          The leading dimension of the array A.  LDA >= max(1,M). */

    /*  IPIV    (output) INTEGER array, dimension (min(M,N)) */
    /*          The pivot indices; for 1 <= i <= min(M,N), row i of the */
    /*          matrix was interchanged with row IPIV(i). */

    /*  INFO    (output) INTEGER */
    /*          = 0:  successful exit */
    /*          < 0:  if INFO = -i, the i-th argument had an illegal value */
    /*          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization */
    /*                has been completed, but the factor U is exactly */
    /*                singular, and division by zero will occur if it is used */
    /*                to solve a system of equations. */

    /*  ===================================================================== */

    /*     .. Parameters .. */
    /*     .. */
    /*     .. Local Scalars .. */
    /*     .. */
    /*     .. External Subroutines .. */
    /*     .. */
    /*     .. External Functions .. */
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
        console.error(`SGETRF parameter ${i__1}`,);
        return 0;
    }

    /*     Quick return if possible */

    if (m == 0 || n == 0)
    {
        return 0;
    }

    /*     Determine the block size for this environment. */

    nb = ilaenv(c__1, "SGETRF", " ", m, n, c_n1, c_n1);
    if (nb <= 1 || nb >= min(m, n))
    {

        /*        Use unblocked code. */

        sgetf2(m, n, a, lda, ipiv, info);
    } else
    {

        /*        Use blocked code. */

        i__1 = min(m, n);
        i__2 = nb;
        for (j = 1; i__2 < 0 ? j >= i__1 : j <= i__1; j += i__2)
        {
            /* Computing MIN */
            i__3 = min(m, n) - j + 1;
            jb = min(i__3, nb);

            /*           Factor diagonal and subdiagonal blocks and test for exact */
            /*           singularity. */

            i__3 = m - j + 1;
            sgetf2(i__3, jb, offset(a, a_offset + j + j * a_dim1), lda, offset(ipiv, ipiv_offset + j), iinfo);

            /*           Adjust INFO and the pivot indices. */

            if (info.value == 0 && iinfo.value > 0)
            {
                info.value = iinfo.value + j - 1;
            }
            /* Computing MIN */
            i__4 = m;
            i__5 = j + jb - 1;
            i__3 = min(i__4, i__5);
            for (i__ = j; i__ <= i__3; ++i__)
            {
                ipiv[ipiv_offset + i__] = j - 1 + ipiv[ipiv_offset + i__];
                /* L10: */
            }

            /*           Apply interchanges to columns 1:J-1. */

            i__3 = j - 1;
            i__4 = j + jb - 1;
            slaswp(i__3, a, lda, j, i__4, ipiv, c__1);

            if (j + jb <= n)
            {

                /*              Apply interchanges to columns J+JB:N. */

                i__3 = n - j - jb + 1;
                i__4 = j + jb - 1;
                slaswp(i__3, offset(a, a_offset + (j + jb) * a_dim1 + 1), lda, j, i__4,
                    offset(ipiv, ipiv_offset + 1), c__1);

                /*              Compute block row of U. */

                i__3 = n - j - jb + 1;
                strsm("Left", "Lower", "No transpose", "Unit", jb, i__3,
                    c_b16, offset(a, a_offset + j + j * a_dim1), lda, offset(a, a_offset + j + (j + jb) *
                        a_dim1), lda);
                if (j + jb <= m)
                {

                    /*                 Update trailing submatrix. */

                    i__3 = m - j - jb + 1;
                    i__4 = n - j - jb + 1;
                    sgemm("No transpose", "No transpose", i__3, i__4, jb,
                        c_b19, offset(a, a_offset + j + jb + j * a_dim1), lda, offset(a, a_offset + j + (j +
                            jb) * a_dim1), lda, c_b16, offset(a, a_offset + j + jb + (j + jb) *
                            a_dim1), lda);
                }
            }
            /* L20: */
        }
    }
    return 0;

    /*     End of SGETRF */

}; /* sgetrf_ */
