import {LapackInfo} from "../info";
import {slaswp} from "../blas/slaswp";
import {strsm} from "../blas/strsm";
import {offset} from "../blas/offset";

const c__1 = 1;
const c_b12 = 1.;
const c_n1 = -1;

const {max} = Math;

export const sgetrs = (trans: string, n: number, nrhs: number, a: Float32Array,
                       lda: number, ipiv: Int32Array, b: Float32Array, ldb: number, info: LapackInfo): number =>
{
    /* System generated locals */
    let a_dim1: number, a_offset: number, b_dim1: number, b_offset: number, i__1: number;

    /* Local variables */
    let notran: boolean;

    /*  -- LAPACK routine (version 3.2) -- */
    /*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
    /*     November 2006 */

    /*     .. Scalar Arguments .. */
    /*     .. */
    /*     .. Array Arguments .. */
    /*     .. */

    /*  Purpose */
    /*  ======= */

    /*  SGETRS solves a system of linear equations */
    /*     A * X = B  or  A' * X = B */
    /*  with a general N-by-N matrix A using the LU factorization computed */
    /*  by SGETRF. */

    /*  Arguments */
    /*  ========= */

    /*  TRANS   (input) CHARACTER*1 */
    /*          Specifies the form of the system of equations: */
    /*          = 'N':  A * X = B  (No transpose) */
    /*          = 'T':  A'* X = B  (Transpose) */
    /*          = 'C':  A'* X = B  (Conjugate transpose = Transpose) */

    /*  N       (input) INTEGER */
    /*          The order of the matrix A.  N >= 0. */

    /*  NRHS    (input) INTEGER */
    /*          The number of right hand sides, i.e., the number of columns */
    /*          of the matrix B.  NRHS >= 0. */

    /*  A       (input) REAL array, dimension (LDA,N) */
    /*          The factors L and U from the factorization A = P*L*U */
    /*          as computed by SGETRF. */

    /*  LDA     (input) INTEGER */
    /*          The leading dimension of the array A.  LDA >= max(1,N). */

    /*  IPIV    (input) INTEGER array, dimension (N) */
    /*          The pivot indices from SGETRF; for 1<=i<=N, row i of the */
    /*          matrix was interchanged with row IPIV(i). */

    /*  B       (input/output) REAL array, dimension (LDB,NRHS) */
    /*          On entry, the right hand side matrix B. */
    /*          On exit, the solution matrix X. */

    /*  LDB     (input) INTEGER */
    /*          The leading dimension of the array B.  LDB >= max(1,N). */

    /*  INFO    (output) INTEGER */
    /*          = 0:  successful exit */
    /*          < 0:  if INFO = -i, the i-th argument had an illegal value */

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

    b_dim1 = ldb;
    b_offset = -(1 + b_dim1);

    /* Function Body */
    info.value = 0;
    notran = trans === "N";
    if (!notran && trans !== "T" && trans !== "C")
    {
        info.value = -1;
    } else if (n < 0)
    {
        info.value = -2;
    } else if (nrhs < 0)
    {
        info.value = -3;
    } else if (lda < max(1, n))
    {
        info.value = -5;
    } else if (ldb < max(1, n))
    {
        info.value = -8;
    }
    if (info.value != 0)
    {
        i__1 = -(info.value);

        console.error(`SGETRS parameter ${i__1}`,);
        return 0;
    }

    /*     Quick return if possible */

    if (n == 0 || nrhs == 0)
    {
        return 0;
    }

    if (notran)
    {

        /*        Solve A * X = B. */

        /*        Apply row interchanges to the right hand sides. */

        slaswp(nrhs, b, ldb, c__1, n, offset(ipiv, ipiv_offset + 1), c__1);

        /*        Solve L*X = B, overwriting B with X. */

        strsm("Left", "Lower", "No transpose", "Unit", n, nrhs, c_b12, a, lda, b, ldb);

        /*        Solve U*X = B, overwriting B with X. */

        strsm("Left", "Upper", "No transpose", "Non-unit", n, nrhs, c_b12,
            a, lda, b, ldb);
    } else
    {

        /*        Solve A' * X = B. */

        /*        Solve U'*X = B, overwriting B with X. */

        strsm("Left", "Upper", "Transpose", "Non-unit", n, nrhs, c_b12, a, lda, b, ldb);

        /*        Solve L'*X = B, overwriting B with X. */

        strsm("Left", "Lower", "Transpose", "Unit", n, nrhs, c_b12, a, lda, b, ldb);

        /*        Apply row interchanges to the solution vectors. */

        slaswp(nrhs, b, ldb, c__1, n, offset(ipiv, ipiv_offset + 1), c_n1);
    }

    return 0;

    /*     End of SGETRS */

}; /* sgetrs_ */
