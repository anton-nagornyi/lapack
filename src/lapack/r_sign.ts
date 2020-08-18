export const r_sign = (a: number, b: number): number =>
{
    const x = (a >= 0 ? a : -a);
    return b >= 0 ? x : -x;
};