export const fillArray = (a: number[] | number[][], colMajorLayout = true): Float32Array =>
{
    if (a.length === 0) return new Float32Array();
    if (Array.isArray(a[0]))
    {
        const aa = a as number[][];
        const res = new Float32Array(aa.length * aa[0].length);
        let n = 0;
        if ( colMajorLayout )
        {
            for (let j = 0; j < a[0].length; ++j)
            {
                for ( let i = 0; i < a.length; ++i)
                {
                    res[n++] = aa[i][j];
                }
            }
        }
        else
        {
            for ( let i = 0; i < a.length; ++i)
            {
                for (let j = 0; j < a[0].length; ++j)
                {
                    res[n++] = aa[i][j];
                }
            }
        }
        return res;
    }
    else
    {
        return new Float32Array(a as number[]);
    }
};

export const extractMNArray = (a: Float32Array, m: number, n: number = 1, colMajorLayout = true): number[] | number[][] =>
{
    if (n === 1)
    {
        return Array.from(a);
    }
    const res = new Array<Array<number>>(m);
    for (let i = 0; i < m; ++i)
    {
        res[i] = new Array<number>(n);
    }
    let k = 0;
    if (colMajorLayout)
    {
        for ( let j = 0; j < n; ++j )
        {
            for (let i = 0; i < m; ++i)
            {
                res[i][j] = a[k++];
            }
        }
    }
    else
    {
        for (let i = 0; i < m; ++i)
        {
            for ( let j = 0; j < n; ++j )
            {
                res[i][j] = a[k++];
            }
        }
    }
    return res;
};

export const expectToBeCloseToArray = (actual: number[][], expected: number[][]) =>
{
    expect(actual.length).toBe(expected.length);
    if (actual.length > 0)
    {
        for (let i = 0; i < actual.length; ++i)
        {
            for (let j = 0; j < actual[0].length; ++j)
            {
                expect(actual[i][j]).toBeCloseTo(expected[i][j]);
            }
        }
    }
};