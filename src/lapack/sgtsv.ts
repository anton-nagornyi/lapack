import {LapackInfo} from "../info";

const {abs, max} = Math;

export const sgtsv = (n: number, nrhs: number, dl: Float32Array, d__: Float32Array,
                      du: Float32Array, b: Float32Array, ldb: number, info: LapackInfo): number =>
{
    /* System generated locals */
    let b_dim1: number, b_offset: number, i__1: number, i__2: number;
    let r__1: number, r__2: number;

    /* Local variables */
    let i__: number, j: number;
    let fact: number, temp: number;

    /*  -- LAPACK routine (version 3.2) -- */
    /*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
    /*     November 2006 */

    /*     .. Scalar Arguments .. */
    /*     .. */
    /*     .. Array Arguments .. */
    /*     .. */

    /*  Purpose */
    /*  ======= */

    /*  SGTSV  solves the equation */

    /*     A*X = B, */

    /*  where A is an n by n tridiagonal matrix, by Gaussian elimination with */
    /*  partial pivoting. */

    /*  Note that the equation  A'*X = B  may be solved by interchanging the */
    /*  order of the arguments DU and DL. */

    /*  Arguments */
    /*  ========= */

    /*  N       (input) INTEGER */
    /*          The order of the matrix A.  N >= 0. */

    /*  NRHS    (input) INTEGER */
    /*          The number of right hand sides, i.e., the number of columns */
    /*          of the matrix B.  NRHS >= 0. */

    /*  DL      (input/output) REAL array, dimension (N-1) */
    /*          On entry, DL must contain the (n-1) sub-diagonal elements of */
    /*          A. */

    /*          On exit, DL is overwritten by the (n-2) elements of the */
    /*          second super-diagonal of the upper triangular matrix U from */
    /*          the LU factorization of A, in DL(1), ..., DL(n-2). */

    /*  D       (input/output) REAL array, dimension (N) */
    /*          On entry, D must contain the diagonal elements of A. */

    /*          On exit, D is overwritten by the n diagonal elements of U. */

    /*  DU      (input/output) REAL array, dimension (N-1) */
    /*          On entry, DU must contain the (n-1) super-diagonal elements */
    /*          of A. */

    /*          On exit, DU is overwritten by the (n-1) elements of the first */
    /*          super-diagonal of U. */

    /*  B       (input/output) REAL array, dimension (LDB,NRHS) */
    /*          On entry, the N by NRHS matrix of right hand side matrix B. */
    /*          On exit, if INFO = 0, the N by NRHS solution matrix X. */

    /*  LDB     (input) INTEGER */
    /*          The leading dimension of the array B.  LDB >= max(1,N). */

    /*  INFO    (output) INTEGER */
    /*          = 0: successful exit */
    /*          < 0: if INFO = -i, the i-th argument had an illegal value */
    /*          > 0: if INFO = i, U(i,i) is exactly zero, and the solution */
    /*               has not been computed.  The factorization has not been */
    /*               completed unless i = N. */

    /*  ===================================================================== */

    /*     .. Parameters .. */
    /*     .. */
    /*     .. Local Scalars .. */
    /*     .. */
    /*     .. Intrinsic Functions .. */
    /*     .. */
    /*     .. External Subroutines .. */
    /*     .. */
    /*     .. Executable Statements .. */

    /* Parameter adjustments */
    // --dl;
    // --d__;
    // --du;
    b_dim1 = ldb;
    b_offset = -(1 + b_dim1);
    const dl_offset = -1;
    const d_offset = -1;
    const du_offset = -1;
    
    // b -= b_offset;

    /* Function Body */
    info.value = 0;
    if (n < 0)
    {
        info.value = -1;
    } else if (nrhs < 0)
    {
        info.value = -2;
    } else if (ldb < max(1, n))
    {
        info.value = -7;
    }
    if (info.value != 0)
    {
        i__1 = -(info.value);
        console.error(`SGTSV parameter ${i__1}`,);
        return 0;
    }

    if (n == 0)
    {
        return 0;
    }

    if (nrhs == 1)
    {
        i__1 = n - 2;
        for (i__ = 1; i__ <= i__1; ++i__)
        {
            if ((r__1 = d__[d_offset + i__], abs(r__1)) >= (r__2 = dl[dl_offset + i__], abs(r__2)))
            {

                /*              No row interchange required */

                if (d__[d_offset + i__] != 0)
                {
                    fact = dl[dl_offset + i__] / d__[d_offset + i__];
                    d__[d_offset + i__ + 1] -= fact * du[du_offset + i__];
                    b[b_offset + i__ + 1 + b_dim1] -= fact * b[b_offset + i__ + b_dim1];
                } else
                {
                    info.value = i__;
                    return 0;
                }
                dl[dl_offset + i__] = 0;
            } else
            {

                /*              Interchange rows I and I+1 */

                fact = d__[d_offset + i__] / dl[dl_offset + i__];
                d__[d_offset + i__] = dl[dl_offset + i__];
                temp = d__[d_offset + i__ + 1];
                d__[d_offset + i__ + 1] = du[du_offset + i__] - fact * temp;
                dl[dl_offset + i__] = du[du_offset + i__ + 1];
                du[du_offset + i__ + 1] = -fact * dl[dl_offset + i__];
                du[du_offset + i__] = temp;
                temp = b[b_offset + i__ + b_dim1];
                b[b_offset + i__ + b_dim1] = b[b_offset + i__ + 1 + b_dim1];
                b[b_offset + i__ + 1 + b_dim1] = temp - fact * b[b_offset + i__ + 1 + b_dim1];
            }
            /* L10: */
        }
        if (n > 1)
        {
            i__ = n - 1;
            if ((r__1 = d__[d_offset + i__], abs(r__1)) >= (r__2 = dl[dl_offset + i__], abs(r__2)))
            {
                if (d__[d_offset + i__] != 0)
                {
                    fact = dl[dl_offset + i__] / d__[d_offset + i__];
                    d__[d_offset + i__ + 1] -= fact * du[du_offset + i__];
                    b[b_offset + i__ + 1 + b_dim1] -= fact * b[b_offset + i__ + b_dim1];
                } else
                {
                    info.value = i__;
                    return 0;
                }
            } else
            {
                fact = d__[d_offset + i__] / dl[dl_offset + i__];
                d__[d_offset + i__] = dl[dl_offset + i__];
                temp = d__[d_offset + i__ + 1];
                d__[d_offset + i__ + 1] = du[du_offset + i__] - fact * temp;
                du[du_offset + i__] = temp;
                temp = b[b_offset + i__ + b_dim1];
                b[b_offset + i__ + b_dim1] = b[b_offset + i__ + 1 + b_dim1];
                b[b_offset + i__ + 1 + b_dim1] = temp - fact * b[b_offset + i__ + 1 + b_dim1];
            }
        }
        if (d__[d_offset + n] == 0)
        {
            info.value = n;
            return 0;
        }
    } else
    {
        i__1 = n - 2;
        for (i__ = 1; i__ <= i__1; ++i__)
        {
            if ((r__1 = d__[d_offset + i__], abs(r__1)) >= (r__2 = dl[dl_offset + i__], abs(r__2)))
            {

                /*              No row interchange required */

                if (d__[d_offset + i__] != 0)
                {
                    fact = dl[dl_offset + i__] / d__[d_offset + i__];
                    d__[d_offset + i__ + 1] -= fact * du[du_offset + i__];
                    i__2 = nrhs;
                    for (j = 1; j <= i__2; ++j)
                    {
                        b[b_offset + i__ + 1 + j * b_dim1] -= fact * b[b_offset + i__ + j * b_dim1];
                        /* L20: */
                    }
                } else
                {
                    info.value = i__;
                    return 0;
                }
                dl[dl_offset + i__] = 0;
            } else
            {

                /*              Interchange rows I and I+1 */

                fact = d__[d_offset + i__] / dl[dl_offset + i__];
                d__[d_offset + i__] = dl[dl_offset + i__];
                temp = d__[d_offset + i__ + 1];
                d__[d_offset + i__ + 1] = du[du_offset + i__] - fact * temp;
                dl[dl_offset + i__] = du[du_offset + i__ + 1];
                du[du_offset + i__ + 1] = -fact * dl[dl_offset + i__];
                du[du_offset + i__] = temp;
                i__2 = nrhs;
                for (j = 1; j <= i__2; ++j)
                {
                    temp = b[b_offset + i__ + j * b_dim1];
                    b[b_offset + i__ + j * b_dim1] = b[b_offset + i__ + 1 + j * b_dim1];
                    b[b_offset + i__ + 1 + j * b_dim1] = temp - fact * b[b_offset + i__ + 1 + j *
                    b_dim1];
                    /* L30: */
                }
            }
            /* L40: */
        }
        if (n > 1)
        {
            i__ = n - 1;
            if ((r__1 = d__[d_offset + i__], abs(r__1)) >= (r__2 = dl[dl_offset + i__], abs(r__2)))
            {
                if (d__[d_offset + i__] != 0)
                {
                    fact = dl[dl_offset + i__] / d__[d_offset + i__];
                    d__[d_offset + i__ + 1] -= fact * du[du_offset + i__];
                    i__1 = nrhs;
                    for (j = 1; j <= i__1; ++j)
                    {
                        b[b_offset + i__ + 1 + j * b_dim1] -= fact * b[b_offset + i__ + j * b_dim1];
                        /* L50: */
                    }
                } else
                {
                    info.value = i__;
                    return 0;
                }
            } else
            {
                fact = d__[d_offset + i__] / dl[dl_offset + i__];
                d__[d_offset + i__] = dl[dl_offset + i__];
                temp = d__[d_offset + i__ + 1];
                d__[d_offset + i__ + 1] = du[du_offset + i__] - fact * temp;
                du[du_offset + i__] = temp;
                i__1 = nrhs;
                for (j = 1; j <= i__1; ++j)
                {
                    temp = b[b_offset + i__ + j * b_dim1];
                    b[b_offset + i__ + j * b_dim1] = b[b_offset + i__ + 1 + j * b_dim1];
                    b[b_offset + i__ + 1 + j * b_dim1] = temp - fact * b[b_offset + i__ + 1 + j *
                    b_dim1];
                    /* L60: */
                }
            }
        }
        if (d__[d_offset + n] == 0)
        {
            info.value = n;
            return 0;
        }
    }

    /*     Back solve with the matrix U from the factorization. */

    if (nrhs <= 2)
    {
        for (j = 1; j <= nrhs; ++j)
        {
            b[b_offset + n + j * b_dim1] /= d__[d_offset + n];
            if (n > 1)
            {
                b[b_offset + n - 1 + j * b_dim1] = (b[b_offset + n - 1 + j * b_dim1] - du[du_offset + n - 1] * b[b_offset + 
                n + j * b_dim1]) / d__[d_offset + n - 1];
            }
            for (i__ = n - 2; i__ >= 1; --i__)
            {
                b[b_offset + i__ + j * b_dim1] = (b[b_offset + i__ + j * b_dim1] - du[du_offset + i__] * b[b_offset + i__ + 1
                + j * b_dim1] - dl[dl_offset + i__] * b[b_offset + i__ + 2 + j * b_dim1]) / d__[d_offset + 
                    i__];
                /* L80: */
            }
        }
    } else
    {
        i__1 = nrhs;
        for (j = 1; j <= i__1; ++j)
        {
            b[b_offset + n + j * b_dim1] /= d__[d_offset + n];
            if (n > 1)
            {
                b[b_offset + n - 1 + j * b_dim1] = (b[b_offset + n - 1 + j * b_dim1] - du[du_offset + n - 1]
                    * b[b_offset + n + j * b_dim1]) / d__[d_offset + n - 1];
            }
            for (i__ = n - 2; i__ >= 1; --i__)
            {
                b[b_offset + i__ + j * b_dim1] = (b[b_offset + i__ + j * b_dim1] - du[du_offset + i__] * b[b_offset + i__
                    + 1 + j * b_dim1] - dl[dl_offset + i__] * b[b_offset + i__ + 2 + j * b_dim1])
                    / d__[d_offset + i__];
                /* L90: */
            }
            /* L100: */
        }
    }

    return 0;

    /*     End of SGTSV */

}; /* sgtsv_ */
