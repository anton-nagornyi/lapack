import {LapackInfo} from "../info";
import {sswap} from "../blas/sswap";
import {sger} from "../blas/sger";
import {stbsv} from "../blas/stbsv";
import {sgemv} from "../blas/sgemv";
import {offset} from "../blas/offset";

const c_b7 = -1;
const c__1 = 1;
const c_b23 = 1;

const {min, max} = Math;

export const sgbtrs = (trans: string, n: number, kl: number,
                       ku: number, nrhs: number, ab: Float32Array, ldab: number, ipiv: Int32Array, b: Float32Array,
                       ldb: number, info: LapackInfo) =>
{
    /* System generated locals */
    let ab_dim1: number, ab_offset: number, b_dim1: number, b_offset: number, i__1: number, i__2: number, i__3: number;

    /* Local variables */
    let i__: number, j: number, l: number, kd: number, lm: number;
    let lnoti: boolean, notran: boolean;

    /*  -- LAPACK routine (version 3.2) -- */
    /*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
    /*     November 2006 */

    /*     .. Scalar Arguments .. */
    /*     .. */
    /*     .. Array Arguments .. */
    /*     .. */

    /*  Purpose */
    /*  ======= */

    /*  SGBTRS solves a system of linear equations */
    /*     A * X = B  or  A' * X = B */
    /*  with a general band matrix A using the LU factorization computed */
    /*  by SGBTRF. */

    /*  Arguments */
    /*  ========= */

    /*  TRANS   (input) CHARACTER*1 */
    /*          Specifies the form of the system of equations. */
    /*          = 'N':  A * X = B  (No transpose) */
    /*          = 'T':  A'* X = B  (Transpose) */
    /*          = 'C':  A'* X = B  (Conjugate transpose = Transpose) */

    /*  N       (input) INTEGER */
    /*          The order of the matrix A.  N >= 0. */

    /*  KL      (input) INTEGER */
    /*          The number of subdiagonals within the band of A.  KL >= 0. */

    /*  KU      (input) INTEGER */
    /*          The number of superdiagonals within the band of A.  KU >= 0. */

    /*  NRHS    (input) INTEGER */
    /*          The number of right hand sides, i.e., the number of columns */
    /*          of the matrix B.  NRHS >= 0. */

    /*  AB      (input) REAL array, dimension (LDAB,N) */
    /*          Details of the LU factorization of the band matrix A, as */
    /*          computed by SGBTRF.  U is stored as an upper triangular band */
    /*          matrix with KL+KU superdiagonals in rows 1 to KL+KU+1, and */
    /*          the multipliers used during the factorization are stored in */
    /*          rows KL+KU+2 to 2*KL+KU+1. */

    /*  LDAB    (input) INTEGER */
    /*          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1. */

    /*  IPIV    (input) INTEGER array, dimension (N) */
    /*          The pivot indices; for 1 <= i <= N, row i of the matrix was */
    /*          interchanged with row IPIV(i). */

    /*  B       (input/output) REAL array, dimension (LDB,NRHS) */
    /*          On entry, the right hand side matrix B. */
    /*          On exit, the solution matrix X. */

    /*  LDB     (input) INTEGER */
    /*          The leading dimension of the array B.  LDB >= max(1,N). */

    /*  INFO    (output) INTEGER */
    /*          = 0:  successful exit */
    /*          < 0: if INFO = -i, the i-th argument had an illegal value */

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
    ab_dim1 = ldab;
    ab_offset = 1 + ab_dim1;
    const ipiv_offset = 1;
    b_dim1 = ldb;
    b_offset = 1 + b_dim1;

    /* Function Body */
    info.value = 0;
    notran = trans.startsWith("N");
    if (!notran && !trans.startsWith("T") && !trans.startsWith("C"))
    {
        info.value = -1;
    } else if (n < 0)
    {
        info.value = -2;
    } else if (kl < 0)
    {
        info.value = -3;
    } else if (ku < 0)
    {
        info.value = -4;
    } else if (nrhs < 0)
    {
        info.value = -5;
    } else if (ldab < (kl << 1) + ku + 1)
    {
        info.value = -7;
    } else if (ldb < max(1, n))
    {
        info.value = -10;
    }
    if (info.value != 0)
    {
        i__1 = -(info.value);
        console.error(`SGBTRS parameter ${i__1}`,);
        return 0;
    }

    /*     Quick return if possible */

    if (n == 0 || nrhs == 0)
    {
        return 0;
    }

    kd = ku + kl + 1;
    lnoti = kl > 0;

    if (notran)
    {

        /*        Solve  A*X = B. */

        /*        Solve L*X = B, overwriting B with X. */

        /*        L is represented as a product of permutations and unit lower */
        /*        triangular matrices L = P(1) * L(1) * ... * P(n-1) * L(n-1), */
        /*        where each transformation L(i) is a rank-one modification of */
        /*        the identity matrix. */

        if (lnoti)
        {
            i__1 = n - 1;
            for (j = 1; j <= i__1; ++j)
            {
                /* Computing MIN */
                i__2 = kl;
                i__3 = n - j;
                lm = min(i__2, i__3);
                l = ipiv[j - ipiv_offset];

                if (l != j)
                {
                    sswap(nrhs, offset(b, l + b_dim1 - b_offset), ldb, offset(b, j + b_dim1 - b_offset), ldb);
                }
                sger(lm, nrhs, c_b7, offset(ab, kd + 1 + j * ab_dim1 - ab_offset), c__1,
                    offset(b, j + b_dim1 - b_offset), ldb, offset(b, j + 1 + b_dim1 -b_offset), ldb);
                /* L10: */
            }
        }

        i__1 = nrhs;
        for (i__ = 1; i__ <= i__1; ++i__)
        {

            /*           Solve U*X = B, overwriting B with X. */

            i__2 = kl + ku;
            stbsv("Upper", "No transpose", "Non-unit", n, i__2, ab, ldab,
                offset(b, i__ * b_dim1 + 1 - b_offset), c__1);
            /* L20: */
        }

    } else
    {

        /*        Solve A'*X = B. */

        i__1 = nrhs;
        for (i__ = 1; i__ <= i__1; ++i__)
        {

            /*           Solve U'*X = B, overwriting B with X. */

            i__2 = kl + ku;
            stbsv("Upper", "Transpose", "Non-unit", n, i__2, ab,
                ldab, offset(b, i__ * b_dim1 + 1 - b_offset), c__1);
            /* L30: */
        }

        /*        Solve L'*X = B, overwriting B with X. */

        if (lnoti)
        {
            for (j = n - 1; j >= 1; --j)
            {
                /* Computing MIN */
                i__1 = kl;
                i__2 = n - j;
                lm = min(i__1, i__2);
                sgemv("Transpose", lm, nrhs, c_b7, offset(b, j + 1 + b_dim1 - b_offset), ldb,
                    offset(ab, kd + 1 + j * ab_dim1 - ab_offset), c__1, c_b23, offset(b, j +
                    b_dim1 - b_offset), ldb);
                l = ipiv[j- ipiv_offset];
                if (l != j)
                {
                    sswap(nrhs, offset(b, l + b_dim1 - b_offset), ldb, offset(b, j + b_dim1 - b_offset), ldb);
                }
                /* L40: */
            }
        }
    }
    return 0;

    /*     End of SGBTRS */

}; /* sgbtrs_ */
