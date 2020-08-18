const {max} = Math;

export const sger = (m: number, n: number, alpha: number, x: Float32Array,
                     incx: number, y: Float32Array, incy: number, a: Float32Array, lda: number): number =>
{
    /* System generated locals */
    let a_dim1: number, a_offset: number, i__1: number, i__2: number;

    /* Local variables */
    let i__: number, j: number, ix: number, jy: number, kx: number, info: number;
    let temp: number;

    /*     .. Scalar Arguments .. */
    /*     .. */
    /*     .. Array Arguments .. */
    /*     .. */

    /*  Purpose */
    /*  ======= */

    /*  SGER   performs the rank 1 operation */

    /*     A := alpha*x*y' + A, */

    /*  where alpha is a scalar, x is an m element vector, y is an n element */
    /*  vector and A is an m by n matrix. */

    /*  Arguments */
    /*  ========== */

    /*  M      - INTEGER. */
    /*           On entry, M specifies the number of rows of the matrix A. */
    /*           M must be at least zero. */
    /*           Unchanged on exit. */

    /*  N      - INTEGER. */
    /*           On entry, N specifies the number of columns of the matrix A. */
    /*           N must be at least zero. */
    /*           Unchanged on exit. */

    /*  ALPHA  - REAL            . */
    /*           On entry, ALPHA specifies the scalar alpha. */
    /*           Unchanged on exit. */

    /*  X      - REAL             array of dimension at least */
    /*           ( 1 + ( m - 1 )*abs( INCX ) ). */
    /*           Before entry, the incremented array X must contain the m */
    /*           element vector x. */
    /*           Unchanged on exit. */

    /*  INCX   - INTEGER. */
    /*           On entry, INCX specifies the increment for the elements of */
    /*           X. INCX must not be zero. */
    /*           Unchanged on exit. */

    /*  Y      - REAL             array of dimension at least */
    /*           ( 1 + ( n - 1 )*abs( INCY ) ). */
    /*           Before entry, the incremented array Y must contain the n */
    /*           element vector y. */
    /*           Unchanged on exit. */

    /*  INCY   - INTEGER. */
    /*           On entry, INCY specifies the increment for the elements of */
    /*           Y. INCY must not be zero. */
    /*           Unchanged on exit. */

    /*  A      - REAL             array of DIMENSION ( LDA, n ). */
    /*           Before entry, the leading m by n part of the array A must */
    /*           contain the matrix of coefficients. On exit, A is */
    /*           overwritten by the updated matrix. */

    /*  LDA    - INTEGER. */
    /*           On entry, LDA specifies the first dimension of A as declared */
    /*           in the calling (sub) program. LDA must be at least */
    /*           max( 1, m ). */
    /*           Unchanged on exit. */


    /*  Level 2 Blas routine. */

    /*  -- Written on 22-October-1986. */
    /*     Jack Dongarra, Argonne National Lab. */
    /*     Jeremy Du Croz, Nag Central Office. */
    /*     Sven Hammarling, Nag Central Office. */
    /*     Richard Hanson, Sandia National Labs. */


    /*     .. Parameters .. */
    /*     .. */
    /*     .. Local Scalars .. */
    /*     .. */
    /*     .. External Subroutines .. */
    /*     .. */
    /*     .. Intrinsic Functions .. */
    /*     .. */

    /*     Test the input parameters. */

    /* Parameter adjustments */
    a_dim1 = lda;
    a_offset = -(1 + a_dim1);
    const x_offset = -1;
    const y_offset = -1;

    /* Function Body */
    info = 0;
    if (m < 0)
    {
        info = 1;
    } else if (n < 0)
    {
        info = 2;
    } else if (incx == 0)
    {
        info = 5;
    } else if (incy == 0)
    {
        info = 7;
    } else if (lda < max(1, m))
    {
        info = 9;
    }
    if (info != 0)
    {
        console.error(`SGER parameter ${info}`);
        return 0;
    }

    // console.log("SGER", n);
    // const res = new Array<number>();
    // for (let i = 0; i < m; i += incx)
    // {
    //     res.push(x[i]);
    // }
    // console.log(res.join(", "));
    // res.length = 0;
    // for (let i = 0; i < n; i += incy)
    // {
    //     res.push(y[i]);
    // }
    // console.log(res.join(", "));
    // console.log("SGER ---");
    /*     Quick return if possible. */

    if (m == 0 || n == 0 || alpha == 0)
    {
        return 0;
    }

    /*     Start the operations. In this version the elements of A are */
    /*     accessed sequentially with one pass through A. */

    if (incy > 0)
    {
        jy = 1;
    } else
    {
        jy = 1 - (n - 1) * incy;
    }
    if (incx == 1)
    {
        i__1 = n;
        for (j = 1; j <= i__1; ++j)
        {
            if (y[y_offset + jy] != 0)
            {
                temp = alpha * y[y_offset + jy];
                i__2 = m;
                for (i__ = 1; i__ <= i__2; ++i__)
                {
                    a[a_offset + i__ + j * a_dim1] += x[x_offset + i__] * temp;
                    /* L10: */
                }
            }
            jy += incy;
            /* L20: */
        }
    } else
    {
        if (incx > 0)
        {
            kx = 1;
        } else
        {
            kx = 1 - (m - 1) * incx;
        }
        i__1 = n;
        for (j = 1; j <= i__1; ++j)
        {
            if (y[y_offset + jy] != 0)
            {
                temp = alpha * y[y_offset + jy];
                ix = kx;
                i__2 = m;
                for (i__ = 1; i__ <= i__2; ++i__)
                {
                    a[a_offset + i__ + j * a_dim1] += x[x_offset + ix] * temp;
                    ix += incx;
                    /* L30: */
                }
            }
            jy += incy;
            /* L40: */
        }
    }

    return 0;

    /*     End of SGER  . */

}; /* sger_ */
