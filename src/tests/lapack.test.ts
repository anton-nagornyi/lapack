import {LapackInfo} from "../info";
import {sgtsv} from "../lapack/sgtsv";
import {sgbsv} from "../lapack/sgbsv";
import {sgbtrf} from "../lapack/sgbtrf";
import {sger} from "../blas/sger";
import {sscal} from "../blas/sscal";
import {sgbtf2} from "../lapack/sgbtf2";
import {slamch} from "../lapack/slamch";
import {slange} from "../lapack/slange";
import {sgetrf} from "../lapack/sgetrf";
import {sgetrs} from "../lapack/sgetrs";
import {saxpy} from "../blas/saxpy";
import {expectToBeCloseToArray, extractMNArray, fillArray, flatten, roundTo} from "./utils";
import {strtri} from "../lapack/strtri";
import {sgetri} from "../lapack/sgetri";

describe("sgtsv", () =>
{
    it("B is 1d", () =>
    {
        const dl = [-1.0, 2.0];
        const du = [1.0, 4.0];
        const d = [2.0, 7.0, -3.0];
        const b = [4.0, 25.0, -5.0];
        const b1 = new Float32Array(b);
        const nrhs = 1;
        const info = new LapackInfo();
        sgtsv(b.length, nrhs, new Float32Array(dl), new Float32Array(d), new Float32Array(du), b1, b.length, info);
        const res = new Array<number>(b.length);
        for (let i = 0; i < b.length; ++i)
        {
            res[i] = b1[i];
        }
        expect(Array.from(b1)).toEqual([1.0, 2.0, 3.0]);
    });
    it("B is 3d", () =>
    {
        const dl = [3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0];
        const du = [4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0];
        const d = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0];
        const b = [
            [1.0, 1.0, 1.0],
            [1.0, -1.0, 2.0],
            [1.0, 1.0, 3.0],
            [1.0, -1.0, 4.0],
            [1.0, 1.0, 5.0],
            [1.0, -1.0, 6.0],
            [1.0, 1.0, 7.0],
            [1.0, -1.0, 8.0],
            [1.0, 1.0, 9.0],
        ];
        const b1 = new Float32Array(b.length * b[0].length);
        const bb = b as number[][];
        for (let i = 0; i < bb.length; ++i)
        {
            for (let j = 0; j < bb[i].length; ++j)
            {
                b1[i * bb[i].length + j] = bb[i][j];
            }
        }
        const nrhs = b[0].length;
        const info = new LapackInfo();

        sgtsv(b.length, nrhs, new Float32Array(dl), new Float32Array(d), new Float32Array(du), b1, b.length, info);

        const dd = b[0].length;
        const res = new Array<Array<number>>(b.length);
        for (let i = 0; i < b.length; ++i)
        {
            res[i] = new Array<number>(dd);
            for (let j = 0; j < dd; ++j)
            {
                res[i][j] = b1[i * dd + j];
            }
        }

        expectToBeCloseToArray(res, [
            [1.98011, -0.245028, -1.17383],
            [0.727227, 0.948563, -1.03256],
            [0.0467182, 1.01274, -0.038224],
            [5.68859, -1.17215, -4.2234],
            [2.93496, 2.68381, -2.62217],
            [-0.107316, 2.24346, -0.730378],
            [8.16923, -1.79231, -5.42885],
            [4.45144, 3.20877, -4.39078],
            [0.691113, 3.3703, -1.11091]
        ]);
    });
});
it("sger", () =>
{
    const alpha = -1;
    const x = [2., 3., 4., 5., 7., 8., 9.];
    const y = [1., 2., 3., 6., 7., 8., 9.];
    const incx = 1;
    const incy = 1;
    const a = new Float32Array(x.length * y.length);
    sger(x.length, y.length, alpha, new Float32Array(x), incx, new Float32Array(y), incy, a, x.length);

    const a1 = new Array<Array<number>>(y.length);
    for (let j = 0; j < x.length; ++j)
    {
        a1[j] = new Array<number>(x.length);
        for (let i = 0; i < a1.length; ++i)
        {
            a1[j][i] = a[j * x.length + i];
        }
    }
    // console.log(a1);
    expectToBeCloseToArray(a1, [ [ -2, -3, -4, -5, -7, -8, -9 ],
        [ -4, -6, -8, -10, -14, -16, -18 ],
        [ -6, -9, -12, -15, -21, -24, -27 ],
        [ -12, -18, -24, -30, -42, -48, -54 ],
        [ -14, -21, -28, -35, -49, -56, -63 ],
        [ -16, -24, -32, -40, -56, -64, -72 ],
        [ -18, -27, -36, -45, -63, -72, -81 ] ]);
});
it("sger2", () =>
{
    const m = 2;
    const n = 3;
    let lda = 5;
    const incx = 2;
    const incy = 1;
    const alpha = 0.5;

    const len_x = 10;
    const len_y = 10;
    const rmaxa = m + 1;
    const cmaxa = n;

    lda = rmaxa;
    const x = new Array<number>(len_x);
    const y = new Array<number>(len_y);
    const a = new Float32Array(rmaxa * cmaxa);

    for (let i = 0; i < 10; i++)
    {
        x[i] = 1.0;
        y[i] = 1.0;
    }

    for (let i = 0; i < m; i++)
    {
        for (let j = 0; j < n; j++)
        {
            a[i + j * lda] = j + 1;
        }
    }

    sger(m, n, alpha, new Float32Array(x), incx, new Float32Array(y), incy, a, lda);

    // console.log(a.join(","));
    const a1 = new Array<Array<number>>(m);
    for (let i = 0; i < m; ++i)
    {
        a1[i] = new Array<number>(n);
        for (let j = 0; j < n; ++j)
        {
            a1[i][j] = a[j * n + i];
        }
    }
    expectToBeCloseToArray(a1, [ [ 1.5, 2.5, 3.5 ], [ 1.5, 2.5, 3.5 ] ]);
    // console.log(a1);
});
describe("sgbtrf", () =>
{
    it("1", () =>
    {
        const a1 = [
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,],
            [0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0,],
            [0.0, 0.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0,],
            [0.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0,],

            [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,],

            [2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 0.0],
            [3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 0.0, 0.0]
        ];
        const a = new Array<Array<number>>(a1[0].length);
        for (let j = 0; j < a1[0].length; ++j)
        {
            a[j] = new Array<number>(a1.length);
            for (let i = 0; i < a1.length; ++i)
            {
                a[j][i] = a1[i][j];
            }
        }
        const m = 9;
        const n = 9;
        const l = 2;
        const u = 3;
        const lda = 8;

        const ipiv = new Int32Array(Math.min(m, n));
        const info = new LapackInfo();
        const ab = new Float32Array(flatten(a));
        // console.log(ab);
        sgbtf2(m, n, l, u, ab, lda, ipiv, info);
        // console.log(info.value, ipiv.join(", "), ab.join(", "));

        const aa = new Array<Array<number>>(lda);
        for (let i = 0; i < lda; ++i)
        {
            aa[i] = new Array<number>(m);
            for (let j = 0; j < n; ++j)
            {
                aa[i][j] = roundTo(ab[j * lda + i], 3);
            }
        }
        expectToBeCloseToArray(aa, [
            [0, 0, 0, 0, 0, 4, 4, 4, 4],
            [0, 0, 0, 0, 3, 3, 3, 3, -2.272],
            [0, 0, 0, 2, 2, 2, 2, -4.074, -1.747],
            [0, 0, 1, 1, 1, 1, -4.691, -4.177, 1],
            [0, 2, 2, 2, 2, -4.42, -3.174, 2, 0.927],
            [3, 3, 3, 3, -3.617, -5.096, 3, -1.547, 3.038],
            [0.667, 0.444, 0.519, 0.568, -0.334, 0.327, 0.044, -0.79, 0],
            [0.333, -0.111, 0.593, 0.247, -0.829, -0.589, -0.618, -0, 0]
        ]);
        expect(ipiv.join(",")).toBe("3,4,5,6,5,6,9,8,9");
    })
});
it("sscal", () =>
{
    const a = 2;
    const x = new Float32Array([1., 2., 3., 4.]);
    sscal(x.length, a, x, 1);
    expect(x.join(",")).toBe("2,4,6,8");
    // console.log(x.join(","));
});
describe("sgbsv", () =>
{
    it("B is 1d", () =>
    {
        const a1 = [
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,],
            [0.0, 0.0, 0.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0,],
            [0.0, 0.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0,],
            [0.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0,],

            [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,],

            [2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 0.0],
            [3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 0.0, 0.0]
        ];
        const ab = fillArray(a1);

        const b = [
            [1.0, 1.0, 1.0],
            [1.0, -1.0, 2.0],
            [1.0, 1.0, 3.0],
            [1.0, -1.0, 4.0],
            [1.0, 1.0, 5.0],
            [1.0, -1.0, 6.0],
            [1.0, 1.0, 7.0],
            [1.0, -1.0, 8.0],
            [1.0, 1.0, 9.0],
        ];
        const b1 = fillArray(b);
        const nrhs = 3;
        const info = new LapackInfo();
        const l = 2;
        const u = 3;

        const ipiv = new Int32Array(b.length);

        sgbsv(b.length, l, u, nrhs, ab, a1.length, ipiv, b1, b.length, info);

        const res = extractMNArray(b1, b.length, b[0].length) as number[][];
        expectToBeCloseToArray(res, [
            [0.629, 1.149, 4.713],
            [-0.025, 1.870, -0.726],
            [0.523, 0.615, 2.946],
            [-0.286, -1.434, -2.774],
            [-0.104, -0.524, -1.067],
            [-0.118, -0.591, -0.970],
            [0.220, -0.896, 2.027],
            [-0.079, 1.604, -0.340],
            [0.496, 0.480, 3.599]
        ]);
    });
});

it("slamch", () =>
{
    const e = slamch("E");
    const s = slamch("S");
    const b = slamch("B");
    const p = slamch("P");
    const n = slamch("N");
    const r = slamch("R");
    const m = slamch("M");
    const u = slamch("U");
    const l = slamch("L");
    const o = slamch("O");

    expect(e).toBe(5.960464477539063e-8);
    expect(s).toBe(1.1754943508222875e-38);
    expect(b).toBe(2);
    expect(p).toBe(1.1920928955078125e-7);
    expect(n).toBe(24);
    expect(r).toBe(1);
    expect(m).toBe(-125);
    expect(u).toBe(1.1754943508222875e-38);
    expect(l).toBe(128);
    expect(o).toBe(3.4028234663852886e+38);
});

it("slange", () => {
   const a = [
    [ 1.0,  1.0,  1.0,  1.0,  0.0,  0.0,   0.0,   0.0,   0.0 ],
    [ 1.0,  1.0,  1.0,  1.0,  1.0,  0.0,   0.0,   0.0,   0.0 ],
    [ 4.0,  1.0,  1.0,  1.0,  1.0,  1.0,   0.0,   0.0,   0.0 ],
    [ 0.0,  5.0,  1.0,  1.0,  1.0,  1.0,   1.0,   0.0,   0.0 ],
    [ 0.0,  0.0,  6.0,  1.0,  1.0,  1.0,   1.0,   1.0,   0.0 ],
    [ 0.0,  0.0,  0.0,  7.0,  1.0,  1.0,   1.0,   1.0,   1.0 ],
    [ 0.0,  0.0,  0.0,  0.0,  8.0,  1.0,   1.0,   1.0,   1.0 ],
    [ 0.0,  0.0,  0.0,  0.0,  0.0,  9.0,   1.0,   1.0,   1.0 ],
    [ 0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  10.0,  11.0,  12.0 ]];

   const work = new Float32Array(a.length);
   const a1 = fillArray(a);
   const r1 = slange("1", a.length, a[0].length, a1, a[0].length, work);
   const r2 = slange("I", a.length, a[0].length, a1, a[0].length, work);

   expect(r1).toBe(15);
   expect(r2).toBe(33);
});

it("sgetrf", () =>
{
    const a = [
        [ 1.0,  1.2,  1.4,  1.6,  1.8,  2.0,  2.2,  2.4,  2.6 ],
        [ 1.2,  1.0,  1.2,  1.4,  1.6,  1.8,  2.0,  2.2,  2.4 ],
        [ 1.4,  1.2,  1.0,  1.2,  1.4,  1.6,  1.8,  2.0,  2.2 ],
        [ 1.6,  1.4,  1.2,  1.0,  1.2,  1.4,  1.6,  1.8,  2.0 ],
        [ 1.8,  1.6,  1.4,  1.2,  1.0,  1.2,  1.4,  1.6,  1.8 ],
        [ 2.0,  1.8,  1.6,  1.4,  1.2,  1.0,  1.2,  1.4,  1.6 ],
        [ 2.2,  2.0,  1.8,  1.6,  1.4,  1.2,  1.0,  1.2,  1.4 ],
        [ 2.4,  2.2,  2.0,  1.8,  1.6,  1.4,  1.2,  1.0,  1.2 ],
        [ 2.6,  2.4,  2.2,  2.0,  1.8,  1.6,  1.4,  1.2,  1.0 ]
    ];
    const a1 = fillArray(a);
    const info = new LapackInfo();
    const ipiv = new Int32Array(Math.min(a.length, a[0].length));
    sgetrf(a.length, a[0].length, a1, a.length, ipiv, info);

    const res = extractMNArray(a1, a.length, a[0].length) as number[][];
    for (let i = 0; i < res.length; ++i)
    {
        for (let j = 0; j < res[i].length; ++j )
        {
            res[i][j] = roundTo(res[i][j], 3);
        }
    }

    // console.log(res);
    expectToBeCloseToArray(res, [
        [ 2.6,  2.4,  2.2,  2.,  1.8,   1.6,  1.4       ,  1.2       ,  1.],
        [ 0.3846154 ,  0.27692306,  0.5538461 ,  0.83076924,  1.1076922 ,
            1.3846154 ,  1.6615385 ,  1.9384615 ,  2.2153845 ],
        [ 0.4615385 , -0.38888955,  0.40000033,  0.8000004 ,  1.2000008 ,
            1.6000009 ,  2.000001  ,  2.4000013 ,  2.8000016 ],
        [ 0.53846157, -0.33333376,  0.0000003 ,  0.4000001 ,  0.8000001 ,
            1.2       ,  1.6000001 ,  2.        ,  2.4       ],
        [ 0.61538464, -0.27777842,  0.00000089, -0.00000056,  0.40000004,
            0.8       ,  1.2000002 ,  1.6       ,  2.0000002 ],
        [ 0.6923077 , -0.22222266,  0.0000003 ,  0.00000034, -0.0000003 ,
            0.39999998,  0.7999999 ,  1.1999998 ,  1.5999998 ],
        [ 0.7692308 , -0.16666731,  0.00000089, -0.00000056,  0.0000003 ,
            -0.0000003 ,  0.40000004,  0.8       ,  1.2000002 ],
        [ 0.84615386, -0.11111154,  0.0000003 ,  0.00000032, -0.0000003 ,
            -0.        ,  0.        ,  0.39999992,  0.7999999 ],
        [ 0.923077  , -0.05555663,  0.00000119, -0.0000006 ,  0.0000006 ,
            -0.00000089,  0.00000089, -0.00000089,  0.40000033]]);
});

it("sgetrs", () => {
    const a = [ [ 2.6, 2.4, 2.2, 2, 1.8, 1.6, 1.4, 1.2, 1 ],
        [ 0.385, 0.277, 0.554, 0.831, 1.108, 1.385, 1.662, 1.938, 2.215 ],
        [ 0.462, -0.389, 0.4, 0.8, 1.2, 1.6, 2, 2.4, 2.8 ],
        [ 0.538, -0.333, 0, 0.4, 0.8, 1.2, 1.6, 2, 2.4 ],
        [ 0.615, -0.278, 0, -0, 0.4, 0.8, 1.2, 1.6, 2 ],
        [ 0.692, -0.222, 0, 0, -0, 0.4, 0.8, 1.2, 1.6 ],
        [ 0.769, -0.167, 0, -0, 0, -0, 0.4, 0.8, 1.2 ],
        [ 0.846, -0.111, 0, 0, -0, 0, -0, 0.4, 0.8 ],
        [ 0.923, -0.056, 0, 0, 0, -0, 0, -0, 0.4 ] ];
    const b = [
        [ 93.0,  186.0,  279.0,  372.0,  465.0 ],
        [ 84.4,  168.8,  253.2,  337.6,  422.0 ],
        [ 76.6,  153.2,  229.8,  306.4,  383.0 ],
        [ 70.0,  140.0,  210.0,  280.0,  350.0 ],
        [ 65.0,  130.0,  195.0,  260.0,  325.0 ],
        [ 62.0,  124.0,  186.0,  248.0,  310.0 ],
        [ 61.4,  122.8,  184.2,  245.6,  307.0 ],
        [ 63.6,  127.2,  190.8,  254.4,  318.0 ],
        [ 69.0,  138.0,  207.0,  276.0,  345.0 ]
    ];
    const a1 = fillArray(a);

    const b1 = fillArray(b);

    const ipiv = new Int32Array([9, 9, 9, 9, 9, 9, 9, 9, 9]);
    const info = new LapackInfo();
    sgetrs("N", a[0].length, b[0].length, a1, a.length, ipiv, b1, b.length, info);

    const res = extractMNArray(b1, b.length, b[0].length) as number[][];
    for (let i = 0; i < res.length; ++i)
    {
        for (let j = 0; j < res[i].length; ++j )
        {
            res[i][j] = roundTo(res[i][j], 3);
        }
    }
    expectToBeCloseToArray(res,
        [ [ 1.007, 2.013, 3.02, 4.026, 5.033 ],
            [ 2.075, 4.149, 6.224, 8.299, 10.373 ],
            [ 2.994, 5.987, 8.981, 11.974, 14.968 ],
            [ 3.834, 7.668, 11.502, 15.336, 19.17 ],
            [ 5.166, 10.332, 15.498, 20.664, 25.83 ],
            [ 5.834, 11.668, 17.502, 23.336, 29.17 ],
            [ 7.166, 14.332, 21.498, 28.664, 35.83 ],
            [ 7.834, 15.668, 23.502, 31.336, 39.17 ],
            [ 9.083, 18.167, 27.25, 36.334, 45.417 ] ]
    );
    // console.log(res);
});

it("saxpy", () =>
{
   const x = fillArray([1, 2, 3, 4]);
   const y = fillArray([5, 6, 7, 8]);
   saxpy(x.length, 2, x, 1, y, 1);
   expect(y.join(",")).toBe([ 7, 10, 13, 16 ].join(","));

   const x1 = fillArray([1, 2, 3, 4]);
   const y1 = fillArray([5, 6, 7, 8]);
   saxpy(2, 2, x1, 2, y1, 2);
   expect(y1.join(",")).toBe([ 7, 6, 13, 8 ].join(","));
});

it("strtri", () =>
{
    const a = fillArray([
        [1.0,  3.0,  4.0,  5.0,  6.0],
        [0.0,  2.0,  8.0,  9.0,  1.0],
        [0.0,  0.0,  4.0,  8.0,  4.0],
        [0.0,  0.0,  0.0, -2.0,  6.0],
        [0.0,  0.0,  0.0,  0.0, -1.0]
    ]);
    const info = new LapackInfo();
    strtri("U", "N", 5, a,5, info);
    const res = extractMNArray(a, 5, 5) as number[][];
    expectToBeCloseToArray(res, [
        [ 1, -1.5, 2, 3.75, 35 ],
        [ 0, 0.5, -1, -1.75, -14 ],
        [ 0, 0, 0.25, 1, 7 ],
        [ 0, 0, 0, -0.5, -3 ],
        [ 0, 0, 0, 0, -1 ]
    ])
});

it("sgetri", () =>
{
    const a = fillArray([
        [ 4.0000,  1.0000,  1.0000,  1.0000,   1.0000,   1.0000,   0.0000,   0.0000,   0.0000 ],
        [ 0.0000,  5.0000,  1.0000,  1.0000,   1.0000,   1.0000,   1.0000,   0.0000,   0.0000 ],
        [ 0.0000,  0.0000,  6.0000,  1.0000,   1.0000,   1.0000,   1.0000,   1.0000,   0.0000 ],
        [ 0.0000,  0.0000,  0.0000,  7.0000,   1.0000,   1.0000,   1.0000,   1.0000,   1.0000 ],
        [ 0.0000,  0.0000,  0.0000,  0.0000,   8.0000,   1.0000,   1.0000,   1.0000,   1.0000 ],
        [ 0.0000,  0.0000,  0.0000,  0.0000,   0.0000,   9.0000,   1.0000,   1.0000,   1.0000 ],
        [ 0.0000,  0.0000,  0.0000,  0.0000,   0.0000,   0.0000,  10.0000,  11.0000,  12.0000 ],
        [ 0.2500,  0.1500,  0.1000,  0.0714,   0.0536,  -0.0694,  -0.0306,   0.1806,   0.3111 ],
        [ 0.2500,  0.1500,  0.1000,  0.0714,  -0.0714,  -0.0556,  -0.0194,   0.9385,  -0.0031 ]
    ]);
    const ipiv = new Int32Array([3, 4, 5, 6, 7, 8, 9, 8, 9]);
    const info = new LapackInfo();
    let work = new Float32Array(a.length);
    sgetri(9, a, 9, ipiv, work, work.length, info);

    const res = extractMNArray(a, 9, 9) as number[][];

    for (let i = 0; i < res.length; ++i)
    {
        for (let j = 0; j < res[i].length; ++j )
        {
            res[i][j] = roundTo(res[i][j], 3);
        }
    }
    expectToBeCloseToArray(res, [
        [ 0.341, -0.674, -1.159, -0.495, -0.497, -0.355, 6.781, -0.47, -0.501 ],
        [ 56.377, -51.741, -1.159, -0.495, -0.497, -0.355, 6.781, -0.47, -0.501 ],
        [ -54.758, 51.452, 0.827, 0.496, 0.497, 0.212, -6.685, 0.512, 0.501 ],
        [ -0.994, 0.995, -0, -0, -0, -0, 0.001, -0, 0 ],
        [ -0.994, 0.995, -0, -0, -0, -0, 0.001, -0, 0 ],
        [ -0.994, 0.995, -0, -0, -0, -0, -0.124, 0.125, 0 ],
        [ -224.145, 204.27, -9.928, -5.957, -3.971, -2.835, 67.331, -4.912, -5.008 ],
        [ 555.675, -515.964, -9.928, -5.957, -3.971, -2.835, 67.331, -4.912, -5.008 ],
        [ -322.581, 302.742, 4.96, 2.976, 1.984, 1.416, -39.259, 3.075, 3.006 ]
    ])
});