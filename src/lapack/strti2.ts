/* Table of constant values */

import {LapackInfo} from "../info";
import {lsame_} from "./lsame";
import {sscal} from "../blas/sscal";
import {offset} from "../blas/offset";
import {strmv} from "../blas/strmv";

const c__1 = 1;
const {max} = Math;

export const strti2 = (uplo: string, diag: string, n: number, a: Float32Array,
                lda: number, info: LapackInfo): number =>
{
    /* System generated locals */
    let a_dim1: number, a_offset: number, i__1: number, i__2: number;

    /* Local variables */
    let j: number;
    let ajj: number;
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

    /*  STRTI2 computes the inverse of a real upper or lower triangular */
    /*  matrix. */

    /*  This is the Level 2 BLAS version of the algorithm. */

    /*  Arguments */
    /*  ========= */

    /*  UPLO    (input) CHARACTER*1 */
    /*          Specifies whether the matrix A is upper or lower triangular. */
    /*          = 'U':  Upper triangular */
    /*          = 'L':  Lower triangular */

    /*  DIAG    (input) CHARACTER*1 */
    /*          Specifies whether or not the matrix A is unit triangular. */
    /*          = 'N':  Non-unit triangular */
    /*          = 'U':  Unit triangular */

    /*  N       (input) INTEGER */
    /*          The order of the matrix A.  N >= 0. */

    /*  A       (input/output) REAL array, dimension (LDA,N) */
    /*          On entry, the triangular matrix A.  If UPLO = 'U', the */
    /*          leading n by n upper triangular part of the array A contains */
    /*          the upper triangular matrix, and the strictly lower */
    /*          triangular part of A is not referenced.  If UPLO = 'L', the */
    /*          leading n by n lower triangular part of the array A contains */
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
    /*          < 0: if INFO = -k, the k-th argument had an illegal value */

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
        console.error(`STRTI2 parameter ${i__1}`,);
        return 0;
    }

    if (upper)
    {

        /*        Compute inverse of upper triangular matrix. */

        i__1 = n;
        for (j = 1; j <= i__1; ++j)
        {
            if (nounit)
            {
                a[a_offset + j + j * a_dim1] = 1. / a[a_offset + j + j * a_dim1];
                ajj = -a[a_offset + j + j * a_dim1];
            } else
            {
                ajj = -1.;
            }

            /*           Compute elements 1:j-1 of j-th column. */

            i__2 = j - 1;
            strmv("Upper", "No transpose", diag, i__2, a, lda,
                offset(a,a_offset + j * a_dim1 + 1), c__1);
            i__2 = j - 1;
            sscal(i__2, ajj, offset(a, a_offset + j * a_dim1 + 1), c__1);
            /* L10: */
        }
    } else
    {

        /*        Compute inverse of lower triangular matrix. */

        for (j = n; j >= 1; --j)
        {
            if (nounit)
            {
                a[a_offset + j + j * a_dim1] = 1. / a[a_offset + j + j * a_dim1];
                ajj = -a[a_offset + j + j * a_dim1];
            } else
            {
                ajj = -1.;
            }
            if (j < n)
            {

                /*              Compute elements j+1:n of j-th column. */

                i__1 = n - j;
                strmv("Lower", "No transpose", diag, i__1, offset(a, a_offset + j + 1 + (j +
                    1) * a_dim1), lda, offset(a, a_offset + j + 1 + j * a_dim1), c__1);
                i__1 = n - j;
                sscal(i__1, ajj, offset(a, a_offset + j + 1 + j * a_dim1), c__1);
            }
            /* L20: */
        }
    }

    return 0;

    /*     End of STRTI2 */

}; /* strti2_ */