export const sdot = (n: number, sx: Float32Array, incx: number, sy: Float32Array, incy: number): number =>
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

    /*     forms the dot product of two vectors. */
    /*     uses unrolled loops for increments equal to one. */
    /*     jack dongarra, linpack, 3/11/78. */
    /*     modified 12/3/93, array(1) declarations changed to array(*) */


    /*     .. Local Scalars .. */
    /*     .. */
    /*     .. Intrinsic Functions .. */
    /*     .. */
    /* Parameter adjustments */
    const a_offset = -1;

    /* Function Body */
    stemp = 0.;
    if (n <= 0)
    {
        return 0.;
    }
    if (incx == 1 && incy == 1)
    {
        /*        code for both increments equal to 1 */


        /*        clean-up loop */

        m = n % 5;
        if (m != 0)
        {
            i__1 = m;
            for (i__ = 1; i__ <= i__1; ++i__)
            {
                stemp += sx[a_offset + i__] * sy[a_offset + i__];
                /* L30: */
            }
            if (n < 5)
            {
                return stemp;
            }
        }
        mp1 = m + 1;
        i__1 = n;
        for (i__ = mp1; i__ <= i__1; i__ += 5)
        {
            stemp = stemp + sx[a_offset + i__] * sy[a_offset + i__] + sx[a_offset + i__ + 1] * sy[a_offset + i__ + 1] + sx[a_offset + 
            i__ + 2] * sy[a_offset + i__ + 2] + sx[a_offset + i__ + 3] * sy[a_offset + i__ + 3] + sx[a_offset + i__ +
            4] * sy[a_offset + i__ + 4];
            /* L50: */
        }
        return stemp;
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
        stemp += sx[a_offset + ix] * sy[a_offset + iy];
        ix += incx;
        iy += incy;
        /* L10: */
    }
    return stemp;


}; /* sdot_ */