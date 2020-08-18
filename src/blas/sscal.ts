export function sscal(
    n: number,
    sa: number,
    sx: Float32Array,
    incx: number): void {


    //alias
    const sb = 1;


    if (n <= 0 || incx <= 0) return;
    if (sa === 1) return;
    if (incx === 1) {
        /*code for increment equal to 1
        *
        *
        *        clean-up loop*/
        // Munch in batches of 5
        let m = n % 5;
        if (m !== 0) {
            for (let i = 1; i <= m; i++) {
                sx[i - sb] = sa * sx[i - sb];
            }
            if (n < 5) return;
        }
        let mp1 = m + 1;
        for (let i = mp1; i <= n; i += 5) {
            sx[i - sb] = sa * sx[i - sb];
            sx[i + 1 - sb] *= sa;
            sx[i + 2 - sb] *= sa;
            sx[i + 3 - sb] *= sa;
            sx[i + 4 - sb] *= sa;
        }
    }
    else {

        let NINCX = n * incx;
        for (let i = 1; i <= NINCX; i += incx) {
            sx[i - sb] *= sa;
        }
    }
}