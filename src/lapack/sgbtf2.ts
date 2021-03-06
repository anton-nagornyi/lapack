/* Table of constant values */

import {LapackInfo} from "../info";
import {isamax} from "../blas/isamax";
import {sswap} from "../blas/sswap";
import {sscal} from "../blas/sscal";
import {sger} from "../blas/sger";
import {offset} from "../blas/offset";

const c_b9 = -1.;

const {min, max} = Math;

export const sgbtf2 = (m: number, n: number, kl: number, ku: number,
                       ab: Float32Array, ldab: number, ipiv: Int32Array, info: LapackInfo) =>
{
    /* System generated locals */
    let ab_dim1: number, ab_offset: number, i__1: number, i__2: number, i__3: number, i__4: number;
    let r__1: number;

    /* Local variables */
    let i__: number, j: number, km: number, jp: number, ju: number, kv: number;

    /*  -- LAPACK routine (version 3.2) -- */
    /*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
    /*     November 2006 */

    /*     .. Scalar Arguments .. */
    /*     .. */
    /*     .. Array Arguments .. */
    /*     .. */

    /*  Purpose */
    /*  ======= */

    /*  SGBTF2 computes an LU factorization of a real m-by-n band matrix A */
    /*  using partial pivoting with row interchanges. */

    /*  This is the unblocked version of the algorithm, calling Level 2 BLAS. */

    /*  Arguments */
    /*  ========= */

    /*  M       (input) INTEGER */
    /*          The number of rows of the matrix A.  M >= 0. */

    /*  N       (input) INTEGER */
    /*          The number of columns of the matrix A.  N >= 0. */

    /*  KL      (input) INTEGER */
    /*          The number of subdiagonals within the band of A.  KL >= 0. */

    /*  KU      (input) INTEGER */
    /*          The number of superdiagonals within the band of A.  KU >= 0. */

    /*  AB      (input/output) REAL array, dimension (LDAB,N) */
    /*          On entry, the matrix A in band storage, in rows KL+1 to */
    /*          2*KL+KU+1; rows 1 to KL of the array need not be set. */
    /*          The j-th column of A is stored in the j-th column of the */
    /*          array AB as follows: */
    /*          AB(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl) */

    /*          On exit, details of the factorization: U is stored as an */
    /*          upper triangular band matrix with KL+KU superdiagonals in */
    /*          rows 1 to KL+KU+1, and the multipliers used during the */
    /*          factorization are stored in rows KL+KU+2 to 2*KL+KU+1. */
    /*          See below for further details. */

    /*  LDAB    (input) INTEGER */
    /*          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1. */

    /*  IPIV    (output) INTEGER array, dimension (min(M,N)) */
    /*          The pivot indices; for 1 <= i <= min(M,N), row i of the */
    /*          matrix was interchanged with row IPIV(i). */

    /*  INFO    (output) INTEGER */
    /*          = 0: successful exit */
    /*          < 0: if INFO = -i, the i-th argument had an illegal value */
    /*          > 0: if INFO = +i, U(i,i) is exactly zero. The factorization */
    /*               has been completed, but the factor U is exactly */
    /*               singular, and division by zero will occur if it is used */
    /*               to solve a system of equations. */

    /*  Further Details */
    /*  =============== */

    /*  The band storage scheme is illustrated by the following example, when */
    /*  M = N = 6, KL = 2, KU = 1: */

    /*  On entry:                       On exit: */

    /*      *    *    *    +    +    +       *    *    *   u14  u25  u36 */
    /*      *    *    +    +    +    +       *    *   u13  u24  u35  u46 */
    /*      *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56 */
    /*     a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66 */
    /*     a21  a32  a43  a54  a65   *      m21  m32  m43  m54  m65   * */
    /*     a31  a42  a53  a64   *    *      m31  m42  m53  m64   *    * */

    /*  Array elements marked * are not used by the routine; elements marked */
    /*  + need not be set on entry, but are required by the routine to store */
    /*  elements of U, because of fill-in resulting from the row */
    /*  interchanges. */

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

    /*     KV is the number of superdiagonals in the factor U, allowing for */
    /*     fill-in. */

    /* Parameter adjustments */
    ab_dim1 = ldab;
    ab_offset = -(1 + ab_dim1);
    const ipiv_offset = -1;
    /* Function Body */
    kv = ku + kl;

    /*     Test the input parameters. */

    info.value = 0;
    if (m < 0)
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
    } else if (ldab < kl + kv + 1)
    {
        info.value = -6;
    }
    if (info.value != 0)
    {
        i__1 = -(info.value);
        console.error(`SGBTF2 parameter ${i__1}`,);
        return 0;
    }

    /*     Quick return if possible */

    if (m == 0 || n == 0)
    {
        return 0;
    }

    /*     Gaussian elimination with partial pivoting */

    /*     Set fill-in elements in columns KU+2 to KV to zero. */

    for (j = ku + 2; j <= min(kv, n); ++j)
    {
        for (let i = kv - j + 2; i <= kl; ++i)
        {
            ab[ab_offset + i + j * ab_dim1] = 0.;
            /* L10: */
        }
        /* L20: */
    }

    /*     JU is the index of the last column affected by the current stage */
    /*     of the factorization. */

    ju = 1;

    for (j = 1; j <= min(m, n); ++j)
    {

        // console.log(`MM: ${j}`);
        // for (let i = 0; i < ldab; ++i)
        // {
        //     let res = "";
        //     for (let j = 0; j < n; ++j)
        //     {
        //         res += `${ab[j * ldab + i]} `;
        //     }
        //     console.log(res);
        // }
        // console.log(`MM: ${j} ----`);
        /*        Set fill-in elements in column J+KV to zero. */

        if (j + kv <= n)
        {
            for (let i = 1; i <= kl; ++i)
            {
                ab[ab_offset + i + (j + kv) * ab_dim1] = 0.;
                /* L30: */
            }
        }

        /*        Find pivot and test for singularity. KM is the number of */
        /*        subdiagonal elements in the current column. */

        /* Computing MIN */
        km = min(kl, m - j);
        jp = isamax(km + 1, offset(ab, ab_offset + kv + 1 + j * ab_dim1), 1);
        ipiv[ipiv_offset + j] = jp + j - 1;
        // console.log("JP:", jp);
        if (ab[ab_offset + kv + jp + j * ab_dim1] != 0.)
        {
            /* Computing MAX */
            /* Computing MIN */
            i__4 = j + ku + jp - 1;
            i__2 = ju;
            i__3 = min(i__4, n);
            ju = max(i__2, i__3);

            /*           Apply interchange to columns J to JU. */

            if (jp != 1)
            {
                i__2 = ju - j + 1;
                const incx = ldab - 1;
                const incy = ldab - 1;

                if (j === 5)
                {
                    console.log(kv + jp + j * ab_dim1, kv, jp, j * ab_dim1);
                    // console.log(`MM!!: ${j}`);
                    // for (let i = 0; i < ldab; ++i)
                    // {
                    //     let res = "";
                    //     for (let j = 0; j < n; ++j)
                    //     {
                    //         res += `${(ab[j * ldab + i]).toFixed(6)} `;
                    //     }
                    //     console.log(res);
                    // }
                    // console.log(`MM!!: ${j} ----`);
                }
                sswap(i__2, offset(ab, ab_offset + kv + jp + j * ab_dim1), incx, offset(ab, ab_offset + kv + 1 +
                j * ab_dim1), incy);

                // if (j === 5)
                // {
                //     console.log(`MM: ${j}`);
                //     for (let i = 0; i < ldab; ++i)
                //     {
                //         let res = "";
                //         for (let j = 0; j < n; ++j)
                //         {
                //             res += `${(ab[j * ldab + i]).toFixed(6)} `;
                //         }
                //         console.log(res);
                //     }
                //     console.log(`MM: ${j} ----`);
                // }
            }

            if (km > 0)
            {
                /*              Compute multipliers. */

                r__1 = 1. / ab[ab_offset + kv + 1 + j * ab_dim1];
                sscal(km, r__1, offset(ab, ab_offset + kv + 2 + j * ab_dim1), 1);
                /*              Update trailing submatrix within the band. */

                if (ju > j)
                {
                    i__2 = ju - j;
                    i__3 = ldab - 1;
                    i__4 = ldab - 1;
                    sger(km, i__2, c_b9, offset(ab, ab_offset + kv + 2 + j * ab_dim1), 1,
                        offset(ab, ab_offset + kv + (j + 1) * ab_dim1), i__3, offset(ab, ab_offset + kv + 1 +
                        (j + 1) * ab_dim1), i__4);
                }
            }
        } else
        {

            /*           If pivot is zero, set INFO to the index of the pivot */
            /*           unless a zero pivot has already been found. */

            if (info.value == 0)
            {
                info.value = j;
            }
        }
        /* L40: */
    }
    return 0;

    /*     End of SGBTF2 */

}; /* sgbtf2_ */
