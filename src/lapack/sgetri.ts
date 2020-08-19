import {LapackInfo} from "../info";
import {ilaenv} from "./ilaenv";
import {sgemv} from "../blas/sgemv";
import {sgemm} from "../blas/sgemm";
import {strsm} from "../blas/strsm";
import {sswap} from "../blas/sswap";
import {offset} from "../blas/offset";
import {strtri} from "./strtri";

const c__1 = 1;
const c_n1 = -1;
const c__2 = 2;
const c_b20 = -1.;
const c_b22 = 1.;

const {min, max} = Math;

export const sgetri = (n: number, a: Float32Array, lda: number, ipiv: Int32Array,
                       work: Float32Array, lwork: number, info: LapackInfo) =>
{
    /* System generated locals */
    let a_dim1: number, a_offset: number, i__1: number, i__2: number, i__3: number;

    /* Local variables */
    let i__: number, j: number, jb: number, nb: number, jj: number, jp: number, nn: number, iws: number, nbmin: number;

    let ldwork: number, lwkopt: number;
    let lquery: boolean;

    /*  -- LAPACK routine (version 3.2) -- */
    /*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
    /*     November 2006 */

    /*     .. Scalar Arguments .. */
    /*     .. */
    /*     .. Array Arguments .. */
    /*     .. */

    /*  Purpose */
    /*  ======= */

    /*  SGETRI computes the inverse of a matrix using the LU factorization */
    /*  computed by SGETRF. */

    /*  This method inverts U and then computes inv(A) by solving the system */
    /*  inv(A)*L = inv(U) for inv(A). */

    /*  Arguments */
    /*  ========= */

    /*  N       (input) INTEGER */
    /*          The order of the matrix A.  N >= 0. */

    /*  A       (input/output) REAL array, dimension (LDA,N) */
    /*          On entry, the factors L and U from the factorization */
    /*          A = P*L*U as computed by SGETRF. */
    /*          On exit, if INFO = 0, the inverse of the original matrix A. */

    /*  LDA     (input) INTEGER */
    /*          The leading dimension of the array A.  LDA >= max(1,N). */

    /*  IPIV    (input) INTEGER array, dimension (N) */
    /*          The pivot indices from SGETRF; for 1<=i<=N, row i of the */
    /*          matrix was interchanged with row IPIV(i). */

    /*  WORK    (workspace/output) REAL array, dimension (MAX(1,LWORK)) */
    /*          On exit, if INFO=0, then WORK(1) returns the optimal LWORK. */

    /*  LWORK   (input) INTEGER */
    /*          The dimension of the array WORK.  LWORK >= max(1,N). */
    /*          For optimal performance LWORK >= N*NB, where NB is */
    /*          the optimal blocksize returned by ILAENV. */

    /*          If LWORK = -1, then a workspace query is assumed; the routine */
    /*          only calculates the optimal size of the WORK array, returns */
    /*          this value as the first entry of the WORK array, and no error */
    /*          message related to LWORK is issued by XERBLA. */

    /*  INFO    (output) INTEGER */
    /*          = 0:  successful exit */
    /*          < 0:  if INFO = -i, the i-th argument had an illegal value */
    /*          > 0:  if INFO = i, U(i,i) is exactly zero; the matrix is */
    /*                singular and its inverse could not be computed. */

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
    const work_offset = -1;

    /* Function Body */
    info.value = 0;
    nb = ilaenv(c__1, "SGETRI", " ", n, c_n1, c_n1, c_n1);
    lwkopt = n * nb;
    work[work_offset + 1] = lwkopt;
    lquery = lwork == -1;
    if (n < 0)
    {
        info.value = -1;
    } else if (lda < max(1, n))
    {
        info.value = -3;
    } else if (lwork < max(1, n) && !lquery)
    {
        info.value = -6;
    }
    if (info.value != 0)
    {
        i__1 = -(info.value);

        console.error(`SGETRI parameter ${i__1}`,);
        return 0;
    } else if (lquery)
    {
        return 0;
    }

    /*     Quick return if possible */

    if (n == 0)
    {
        return 0;
    }

    /*     Form inv(U).  If INFO > 0 from STRTRI, then U is singular, */
    /*     and the inverse is not computed. */

    strtri("Upper", "Non-unit", n, a, lda, info);
    if (info.value > 0)
    {
        return 0;
    }

    nbmin = 2;
    ldwork = n;
    if (nb > 1 && nb < n)
    {
        /* Computing MAX */
        i__1 = ldwork * nb;
        iws = max(i__1, 1);
        if (lwork < iws)
        {
            nb = lwork / ldwork;
            /* Computing MAX */
            i__1 = 2, i__2 = ilaenv(c__2, "SGETRI", " ", n, c_n1, c_n1, c_n1);
            nbmin = max(i__1, i__2);
        }
    } else
    {
        iws = n;
    }

    /*     Solve the equation inv(A)*L = inv(U) for inv(A). */

    if (nb < nbmin || nb >= n)
    {

        /*        Use unblocked code. */

        for (j = n; j >= 1; --j)
        {

            /*           Copy current column of L to WORK and replace with zeros. */

            i__1 = n;
            for (i__ = j + 1; i__ <= i__1; ++i__)
            {
                work[work_offset + i__] = a[a_offset + i__ + j * a_dim1];
                a[a_offset + i__ + j * a_dim1] = 0.;
                /* L10: */
            }

            /*           Compute current column of inv(A). */

            if (j < n)
            {
                i__1 = n - j;
                sgemv("No transpose", n, i__1, c_b20, offset(a, a_offset + (j + 1) * a_dim1
                    + 1), lda, offset(work, work_offset + j + 1), c__1, c_b22, offset(a, a_offset + j * a_dim1
                    + 1), c__1);
            }
            /* L20: */
        }
    } else
    {

        /*        Use blocked code. */

        nn = (n - 1) / nb * nb + 1;
        i__1 = -nb;
        for (j = nn; i__1 < 0 ? j >= 1 : j <= 1; j += i__1)
        {
            /* Computing MIN */
            i__2 = nb, i__3 = n - j + 1;
            jb = min(i__2, i__3);

            /*           Copy current block column of L to WORK and replace with */
            /*           zeros. */

            i__2 = j + jb - 1;
            for (jj = j; jj <= i__2; ++jj)
            {
                i__3 = n;
                for (i__ = jj + 1; i__ <= i__3; ++i__)
                {
                    work[work_offset + i__ + (jj - j) * ldwork] = a[a_offset + i__ + jj * a_dim1];
                    a[a_offset + i__ + jj * a_dim1] = 0.;
                    /* L30: */
                }
                /* L40: */
            }

            /*           Compute current block column of inv(A). */

            if (j + jb <= n)
            {
                i__2 = n - j - jb + 1;
                sgemm("No transpose", "No transpose", n, jb, i__2, c_b20,
                    offset(a, a_offset + (j + jb) * a_dim1 + 1), lda, offset(work, work_offset + j + jb),
                    ldwork, c_b22, offset(a, a_offset + j * a_dim1 + 1), lda);
            }
            strsm("Right", "Lower", "No transpose", "Unit", n, jb, c_b22,
                offset(work, work_offset + j), ldwork, offset(a, a_offset + j * a_dim1 + 1), lda);
            /* L50: */
        }
    }

    /*     Apply column interchanges. */

    for (j = n - 1; j >= 1; --j)
    {
        jp = ipiv[ipiv_offset + j];
        if (jp != j)
        {
            sswap(n, offset(a, a_offset + j * a_dim1 + 1), c__1, offset(a, a_offset + jp * a_dim1 + 1), c__1);
        }
        /* L60: */
    }

    work[work_offset + 1] = iws;
    return 0;

    /*     End of SGETRI */

}; /* sgetri_ */
