const arr_offset = -1;

export const saxpy = (n: number, sa: number, sx: Float32Array, incx: number,
                      sy: Float32Array, incy: number): number =>
{
    /* System generated locals */
    let i__1: number;

    /* Local variables */
    let i__: number, m: number, ix: number, iy: number, mp1: number;

    /*     .. Scalar Arguments .. */
    /*     .. */
    /*     .. Array Arguments .. */
    /*     .. */

    /*  Purpose */
    /*  ======= */

    /*     SAXPY constant times a vector plus a vector. */
    /*     uses unrolled loop for increments equal to one. */
    /*     jack dongarra, linpack, 3/11/78. */
    /*     modified 12/3/93, array(1) declarations changed to array(*) */


    /*     .. Local Scalars .. */
    /*     .. */
    /*     .. Intrinsic Functions .. */
    /*     .. */
    /* Parameter adjustments */

    /* Function Body */
    if (n <= 0)
    {
        return 0;
    }
    if (sa == 0.)
    {
        return 0;
    }
    if (incx == 1 && incy == 1)
    {
        /*        code for both increments equal to 1 */


        /*        clean-up loop */
        m = n % 4;
        if (m == 0)
        {
            mp1 = m + 1;
            i__1 = n;
            for (i__ = mp1; i__ <= i__1; i__ += 4)
            {
                sy[arr_offset + i__] += sa * sx[arr_offset + i__];
                sy[arr_offset + i__ + 1] += sa * sx[arr_offset + i__ + 1];
                sy[arr_offset + i__ + 2] += sa * sx[arr_offset + i__ + 2];
                sy[arr_offset + i__ + 3] += sa * sx[arr_offset + i__ + 3];
                /* L50: */
            }
        }
        i__1 = m;
        for (i__ = 1; i__ <= i__1; ++i__)
        {
            sy[arr_offset + i__] += sa * sx[arr_offset + i__];
            /* L30: */
        }
        if (n < 4)
        {
            return 0;
        }
        return 0;
    }

    /*        code for unequal increments or equal increments */
    /*          not equal to 1 */
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
        sy[arr_offset + iy] += sa * sx[arr_offset + ix];
        ix += incx;
        iy += incy;
        /* L10: */
    }
    return 0;
}; /* saxpy_ */