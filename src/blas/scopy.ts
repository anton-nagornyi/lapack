export function scopy(
    n: number,
    sx: Float32Array,
    incx: number,
    sy: Float32Array,
    incy: number): void {

    if (n <= 0) return;
    const xb = 1;
    const yb = 1;

    if (incx === 1 && incy === 1) {
        const m = n % 7;
        if (m !== 0) {
            for (let i = 1; i <= m; i++) {
                sy[i - yb] = sx[i - xb];
            }
            if (n < 7) {
                return;
            }
        }
        const mp1 = m + 1;
        for (let i = mp1; i <= n; i++) {
            let kx = i - xb;
            let ky = i - yb;
            // prolly this helped the compiler(fortran) unwind for loops
            sy[ky++] = sx[kx++]; //1
            sy[ky++] = sx[kx++];
            sy[ky++] = sx[kx++];
            sy[ky++] = sx[kx++];
            sy[ky++] = sx[kx++];
            sy[ky++] = sx[kx++];
            sy[ky++] = sx[kx++]; //7
        }
    }
    else {
        let ix = 1;
        let iy = 1;
        if (incx < 0) ix = (-n + 1) * incx + 1;
        if (incy < 0) iy = (-n + 1) * incy + 1;
        for (let i = 1; i <= n; i++) {
            let kx = ix - xb;
            let ky = iy - yb;
            sy[ky] = sx[kx];
            ix += incx;
            iy += incy;
        }
    }
}