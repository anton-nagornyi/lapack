const {abs} = Math;
const sx_offset = -1;

export const sasum = (n: number, sx: Float32Array, incx: number): number =>
{
    /* System generated locals */
    let i__1: number, i__2: number;
    let ret_val: number, r__1: number, r__2: number, r__3: number, r__4: number, r__5: number, r__6: number;

    /* Local variables */
    let i__: number, m: number, mp1: number, nincx: number;
    let stemp: number;

    /*     .. Scalar Arguments .. */
    /*     .. */
    /*     .. Array Arguments .. */
    /*     .. */

    /*  Purpose */
    /*  ======= */

    /*     takes the sum of the absolute values. */
    /*     uses unrolled loops for increment equal to one. */
    /*     jack dongarra, linpack, 3/11/78. */
    /*     modified 3/93 to return if incx .le. 0. */
    /*     modified 12/3/93, array(1) declarations changed to array(*) */


    /*     .. Local Scalars .. */
    /*     .. */
    /*     .. Intrinsic Functions .. */
    /*     .. */
    /* Parameter adjustments */


    /* Function Body */
    ret_val = 0.;
    stemp = 0.;
    if (n <= 0 || incx <= 0)
    {
        return ret_val;
    }
    if (incx == 1)
    {
        return L20(n, stemp, sx);
    }

    /*        code for increment not equal to 1 */

    nincx = n * incx;
    i__1 = nincx;
    i__2 = incx;
    for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2)
    {
        stemp += (r__1 = sx[sx_offset + i__], abs(r__1));
        /* L10: */
    }
    return stemp;

}; /* sasum_ */

const L20 = (n: number, stemp: number, sx: Float32Array): number =>
{
    /*        code for increment equal to 1 */


    /*        clean-up loop */
    const m = n % 6;
    if (m == 0)
    {
        return L40(m, n, stemp, sx);
    }

    let r__1: number;
    for (let i = 1; i <= m; ++i)
    {
        stemp += (r__1 = sx[sx_offset + i], abs(r__1));
        /* L30: */
    }
    if (n < 6)
    {
        return L60(stemp);
    }
    return L40(m, n, stemp, sx);
};
const L40 = (m: number, n: number, stemp: number, sx: Float32Array): number =>
{
    const mp1 = m + 1;
    let r__1: number, r__2: number, r__3: number, r__4: number, r__5: number, r__6: number;
    for (let i = mp1; i <= n; i += 6)
    {
        stemp = stemp + (r__1 = sx[sx_offset + i], abs(r__1)) + (r__2 = sx[sx_offset + i + 1],
            abs(r__2)) + (r__3 = sx[sx_offset + i + 2], abs(r__3)) + (r__4 = sx[sx_offset +
        i + 3], abs(r__4)) + (r__5 = sx[sx_offset + i + 4], abs(r__5)) + (
            r__6 = sx[sx_offset + i + 5], abs(r__6));
        /* L50: */
    }
    return L60(stemp);
};
const L60 = (stemp: number): number =>
{
    return stemp;
};