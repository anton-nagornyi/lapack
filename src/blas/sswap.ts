export const sswap = (n: number, sx: Float32Array, incx: number, sy: Float32Array,
                      incy: number) =>
{
    /* System generated locals */
    let i__1: number;

    /* Local variables */
    let i__: number, m: number, ix: number, iy: number, mp1: number;
    let stemp: number;

    /*     .. Scalar Arguments .. */
    /*     .. */
    /*     .. Array Arguments .. */
    /*     .. */

    /*  Purpose */
    /*  ======= */

    /*     interchanges two vectors. */
    /*     uses unrolled loops for increments equal to 1. */
    /*     jack dongarra, linpack, 3/11/78. */
    /*     modified 12/3/93, array(1) declarations changed to array(*) */


    /*     .. Local Scalars .. */
    /*     .. */
    /*     .. Intrinsic Functions .. */
    /*     .. */
    /* Parameter adjustments */
    const x_offset = -1;
    const y_offset = -1;
    

    /* Function Body */
    if (n <= 0)
    {
        return 0;
    }
    if (incx == 1 && incy == 1)
    {
        return L20(n, x_offset, y_offset, sx, sy);
    } else
    {
        /*       code for unequal increments or equal increments not equal */
        /*         to 1 */
        ix = 1;
        iy = 1;
        if (incx < 0)
        {
            ix = (-(n) + 1) * incx + 1;
        }
        if (incy < 0)
        {
            iy = (-(n) + 1) * incy + 1;
        }
        i__1 = n;
        for (i__ = 1; i__ <= i__1; ++i__)
        {
            stemp = sx[x_offset + ix];
            sx[x_offset + ix] = sy[y_offset + iy];
            sy[y_offset + iy] = stemp;
            ix += incx;
            iy += incy;
            /* L10: */
        }
        return 0;
    }
}; /* sswap_ */

const L20 = (n: number, x_offset: number, y_offset: number, sx: Float32Array, sy: Float32Array): number =>
{
    const m = n % 3;
    if (m == 0)
    {
        return L40(m, n, x_offset, y_offset, sx, sy);
    }
    const i__1 = m;
    for (let i = 1; i <= i__1; ++i)
    {
        const stemp = sx[x_offset + i];
        sx[x_offset + i] = sy[y_offset + i];
        sy[y_offset + i] = stemp;
        /* L30: */
    }
    if (n < 3)
    {
        return 0;
    }
    return L40(m, n, x_offset, y_offset, sx, sy);
};
const L40 = (m: number, n: number, x_offset: number, y_offset: number, sx: Float32Array, sy: Float32Array) =>
{
    const mp1 = m + 1;
    const i__1 = n;
    for (let i = mp1; i <= i__1; i += 3)
    {
        let stemp = sx[i];
        sx[x_offset + i] = sy[y_offset + i];
        sy[y_offset + i] = stemp;
        stemp = sx[x_offset + i + 1];
        sx[x_offset + i + 1] = sy[y_offset + i + 1];
        sy[y_offset + i + 1] = stemp;
        stemp = sx[x_offset + i + 2];
        sx[x_offset + i + 2] = sy[y_offset + i + 2];
        sy[y_offset + i + 2] = stemp;
        /* L50: */
    }
    return 0;
};