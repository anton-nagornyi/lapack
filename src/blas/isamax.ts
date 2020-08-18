const {abs} = Math;

export const isamax = (n: number, sx: Float32Array, incx: number) =>
{
    /* System generated locals */
    let ret_val: number, i__1: number;
    let r__1: number;

    /* Local variables */
    let i__: number, ix: number;
    let smax: number;

    /*     .. Scalar Arguments .. */
    /*     .. */
    /*     .. Array Arguments .. */
    /*     .. */

    /*  Purpose */
    /*  ======= */

    /*     finds the index of element having max. absolute value. */
    /*     jack dongarra, linpack, 3/11/78. */
    /*     modified 3/93 to return if incx .le. 0. */
    /*     modified 12/3/93, array(1) declarations changed to array(*) */


    /*     .. Local Scalars .. */
    /*     .. */
    /*     .. Intrinsic Functions .. */
    /*     .. */
    /* Parameter adjustments */
    const sx_offset = -1;

    /* Function Body */
    ret_val = 0;
    if (n < 1 || incx <= 0)
    {
        return ret_val;
    }
    ret_val = 1;
    if (n == 1)
    {
        return ret_val;
    }
    if (incx == 1)
    {
        /*        code for increment equal to 1 */
        smax = abs(sx[sx_offset + 1]);
        i__1 = n;
        for (i__ = 2; i__ <= i__1; ++i__)
        {
            if ((r__1 = sx[sx_offset + i__], abs(r__1)) <= smax)
            {
                continue;
            }
            ret_val = i__;
            smax = (r__1 = sx[sx_offset + i__], abs(r__1));
        }
        return ret_val;
    }

    /*        code for increment not equal to 1 */

    ix = 1;
    smax = abs(sx[sx_offset + 1]);
    ix += incx;
    i__1 = n;
    for (i__ = 2; i__ <= i__1; ++i__)
    {
        if ((r__1 = sx[sx_offset + ix], abs(r__1)) <= smax)
        {
            ix += incx;
        } else
        {
            ret_val = i__;
            smax = (r__1 = sx[sx_offset + ix], abs(r__1));
        }
        /* L10: */
    }
    return ret_val;
    
}; /* isamax_ */
