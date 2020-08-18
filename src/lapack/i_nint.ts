export const i_nint = (x: number): number =>
{
    return Math.trunc(x >= 0 ? Math.floor(x + .5) : -Math.floor(.5 - x));
};