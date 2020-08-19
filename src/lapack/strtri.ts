/* Table of constant values */

import {LapackInfo} from "../info";
import {lsame_} from "./lsame";
import {ilaenv} from "./ilaenv";
import {strsm} from "../blas/strsm";
import {strmm} from "../blas/strmm";
import {strti2} from "./strti2";
import {offset} from "../blas/offset";

const c__1 = 1;
const c_n1 = -1;
const c_b18 = 1.;
const c_b22 = -1.;

const {min, max} = Math;
export const strtri = (uplo: string, diag: string, n: number, a: Float32Array,
                       lda: number, info: LapackInfo): number =>
{
    /* System generated locals */
    let a__1: string;
    let a_dim1: number, a_offset: number, i__1: number, i__3: number, i__4: number, i__5: number;
    let ch__1: string;

    /* Builtin functions */
    /* Local variables */
    let j: number, jb: number, nb: number, nn: number;
    let upper: boolean;
    let nounit: boolean;

    /*  -- LAPACK routine (version 3.2) -- */
    /*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
    /*     November 2006 */

    /*     .. Scalar Arguments .. */
    /*     .. */
    /*     .. Array Arguments .. */
    /*     .. */

    /*  Purpose */
    /*  ======= */

    /*  STRTRI computes the inverse of a real upper or lower triangular */
    /*  matrix A. */

    /*  This is the Level 3 BLAS version of the algorithm. */

    /*  Arguments */
    /*  ========= */

    /*  UPLO    (input) CHARACTER*1 */
    /*          = 'U':  A is upper triangular; */
    /*          = 'L':  A is lower triangular. */

    /*  DIAG    (input) CHARACTER*1 */
    /*          = 'N':  A is non-unit triangular; */
    /*          = 'U':  A is unit triangular. */

    /*  N       (input) INTEGER */
    /*          The order of the matrix A.  N >= 0. */

    /*  A       (input/output) REAL array, dimension (LDA,N) */
    /*          On entry, the triangular matrix A.  If UPLO = 'U', the */
    /*          leading N-by-N upper triangular part of the array A contains */
    /*          the upper triangular matrix, and the strictly lower */
    /*          triangular part of A is not referenced.  If UPLO = 'L', the */
    /*          leading N-by-N lower triangular part of the array A contains */
    /*          the lower triangular matrix, and the strictly upper */
    /*          triangular part of A is not referenced.  If DIAG = 'U', the */
    /*          diagonal elements of A are also not referenced and are */
    /*          assumed to be 1. */
    /*          On exit, the (triangular) inverse of the original matrix, in */
    /*          the same storage format. */

    /*  LDA     (input) INTEGER */
    /*          The leading dimension of the array A.  LDA >= max(1,N). */

    /*  INFO    (output) INTEGER */
    /*          = 0: successful exit */
    /*          < 0: if INFO = -i, the i-th argument had an illegal value */
    /*          > 0: if INFO = i, A(i,i) is exactly zero.  The triangular */
    /*               matrix is singular and its inverse can not be computed. */

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

    /* Function Body */
    info.value = 0;
    upper = lsame_(uplo, "U");
    nounit = lsame_(diag, "N");
    if (!upper && !lsame_(uplo, "L"))
    {
        info.value = -1;
    } else if (!nounit && !lsame_(diag, "U"))
    {
        info.value = -2;
    } else if (n < 0)
    {
        info.value = -3;
    } else if (lda < max(1, n))
    {
        info.value = -5;
    }
    if (info.value != 0)
    {
        i__1 = -(info.value);
        console.error(`STRTRI parameter ${i__1}`,);
        return 0;
    }

    /*     Quick return if possible */

    if (n == 0)
    {
        return 0;
    }

    /*     Check for singularity if non-unit. */

    if (nounit)
    {
        i__1 = n;
        for (info.value = 1; info.value <= i__1; ++(info.value))
        {
            if (a[a_offset + info.value + info.value * a_dim1] == 0.)
            {
                return 0;
            }
            /* L10: */
        }
        info.value = 0;
    }

    /*     Determine the block size for this environment. */

    /* Writing concatenation */
    nb = ilaenv(c__1, "STRTRI", uplo[0] + diag[0], n, c_n1, c_n1, c_n1);
    if (nb <= 1 || nb >= n)
    {

        /*        Use unblocked code */
        strti2(uplo, diag, n, a, lda, info);
    } else
    {
        /*        Use blocked code */

        if (upper)
        {

            /*           Compute inverse of upper triangular matrix */

            i__1 = n;
            i__3 = nb;
            for (j = 1; i__3 < 0 ? j >= i__1 : j <= i__1; j += i__3)
            {
                /* Computing MIN */
                i__4 = nb, i__5 = n - j + 1;
                jb = min(i__4, i__5);

                /*              Compute rows 1:j-1 of current block column */

                i__4 = j - 1;
                strmm("Left", "Upper", "No transpose", diag, i__4, jb,
                    c_b18, a, lda, offset(a, j * a_dim1 + 1), lda);
                i__4 = j - 1;
                strsm("Right", "Upper", "No transpose", diag, i__4, jb,
                    c_b22, offset(a, j + j * a_dim1), lda, offset(a, j * a_dim1 + 1),
                    lda);

                /*              Compute inverse of current diagonal block */

                strti2("Upper", diag, jb, offset(a, j + j * a_dim1), lda, info);
                /* L20: */
            }
        } else
        {

            /*           Compute inverse of lower triangular matrix */

            nn = (n - 1) / nb * nb + 1;
            i__3 = -nb;
            for (j = nn; i__3 < 0 ? j >= 1 : j <= 1; j += i__3)
            {
                /* Computing MIN */
                i__1 = nb, i__4 = n - j + 1;
                jb = min(i__1, i__4);
                if (j + jb <= n)
                {

                    /*                 Compute rows j+jb:n of current block column */

                    i__1 = n - j - jb + 1;
                    strmm("Left", "Lower", "No transpose", diag, i__1, jb,
                        c_b18, offset(a, j + jb + (j + jb) * a_dim1), lda, offset(a, j
                            + jb + j * a_dim1), lda);
                    i__1 = n - j - jb + 1;
                    strsm("Right", "Lower", "No transpose", diag, i__1, jb,
                        c_b22, offset(a, j + j * a_dim1), lda, offset(a, j + jb + j *
                            a_dim1), lda);
                }

                /*              Compute inverse of current diagonal block */

                strti2("Lower", diag, jb, offset(a, j + j * a_dim1), lda, info);
                /* L30: */
            }
        }
    }

    return 0;

    /*     End of STRTRI */

}; /* strtri_ */