import {Pointer} from "../pointer";

export const pow_ri = (ap: Pointer<number>, bp: Pointer<number>): number =>
{
    const x = ap.value;
    const n = bp.value;
    return Math.fround(x ** n);
};