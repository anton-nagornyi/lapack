import {offset} from "./offset";

const { max } = Math;

export function sger1(
    m: number,
    n: number,
    alpha: number,
    x: Float32Array,
    incx: number,
    y: Float32Array,
    incy: number,
    a: Float32Array,
    lda: number): void {

    let err = 0;
    switch (true) {
        case (m < 0): err = 1; break;
        case (n < 0): err = 2; break;
        case (incx === 0): err = 5; break;
        case (incy === 0): err = 7; break;
        case (lda < max(1, m)): err = 9; break;
        default:
            err = 0;
    }

    if (err) {
        throw new Error('sger');
    }
    const b = 1;
    if (m === 0 || n === 0 || alpha === 0) return;

    let jy = incy < 0 ? 1 - (n - 1) * incy : 1;
    let kx = incx < 0 ? 1 - (m - 1) * incx : 1;
    for (let j = 1; j <= n; j++) {
        if (y[jy - b] !== 0) {
            let temp = alpha * y[jy - b];
            let ix = kx;
            const coords = offset(a, (j - b) * n);//a.colOfEx(j);
            for (let i = 1; i <= m; i++) {
                coords[i - b] += x[ix - b] * temp;
                ix += incx;
            }
        }
        jy += incy;
    }
}