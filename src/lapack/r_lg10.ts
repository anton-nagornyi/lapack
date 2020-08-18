const log10e = 0.43429448190325182765;

export const r_lg10 = (x: number): number =>
{
    return ( log10e * Math.log(x) );
};