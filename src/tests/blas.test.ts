import {strmm} from "../blas/strmm";
import {expectToBeCloseToArray, extractMNArray, fillArray, flatten} from "./utils";
import {strmv} from "../blas/strmv";

it("strmm", () =>
{
    const a = fillArray([
        [3.0,  -1.0,   2.0,   2.0,   1.0],
        [0.0,  -2.0,   4.0,  -1.0,   3.0],
        [0.0,   0.0,  -3.0,   0.0,   2.0],
        [0.0,   0.0,   0.0,   4.0,  -2.0],
        [0.0,   0.0,   0.0,   0.0,   1.0],
        [0.0,   0.0,   0.0,   0.0,   0.0],
        [0.0,   0.0,   0.0,   0.0,   0.0]
    ]);
    const b = fillArray([
        [2.0,  3.0,   1.0],
        [5.0,  5.0,   4.0],
        [0.0,  1.0,   2.0],
        [3.0,  1.0,  -3.0],
        [-1.0, 2.0,   1.0],
        [0.0,  0.0,   0.0]
    ]);
    strmm("L", "U", "N", "N", 5, 3, 1.0, a, 7, b, 6);
    const res = extractMNArray(b, 6, 3) as number[][];
    expectToBeCloseToArray(res, [
        [ 6, 10, -2 ],
        [ -16, -1, 6 ],
        [ -2, 1, -4 ],
        [ 14, 0, -14 ],
        [ -1, 2, 1 ],
        [ 0, 0, 0 ]
    ])
});

it("strmv 1", () =>
{
    const a = fillArray([
        [0.0,   0.0,   0.0,   0.0],
        [1.0,   0.0,   0.0,   0.0],
        [2.0,   3.0,   0.0,   0.0],
        [3.0,   4.0,   3.0,   0.0]
    ]);
    const x = fillArray([1.0, 2.0, 3.0, 4.0]);

    strmv("L" , "N" , "U" , 4 , a , 4 , x , 1);
    expect(x.join(",")).toBe("1,3,11,24");
});

it("strmv 2", () =>
{
    const a = fillArray([
        [0.0,   2.0,   3.0,   2.0],
        [0.0,   0.0,   2.0,   5.0],
        [0.0,   0.0,   0.0,   3.0],
        [0.0,   0.0,   0.0,   0.0]
    ]);
    const x = fillArray([5.0, 4.0, 3.0, 2.0]);

    strmv("U" , "T" , "U" , 4 , a , 4 , x , 1);
    expect(x.join(",")).toBe("5,14,26,41");
});