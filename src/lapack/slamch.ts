import {slamc2} from "./slamc2";
import {pow_ri} from "./pow_ri";
import {Pointer} from "../pointer";

const t = new Pointer<number>(0);
const rnd = new Pointer<number>(0);
const eps = new Pointer<number>(0);
const base = new Pointer<number>(0);
const emin = new Pointer<number>(0);
const prec = new Pointer<number>(0);
const emax = new Pointer<number>(0);

const rmin = new Pointer<number>(0);
const rmax = new Pointer<number>(0);
const sfmin = new Pointer<number>(0);
let first = true;

export const slamch = (cmach: string): number =>
{
    /* Initialized data */



    /* System generated locals */
    let i__1 = new Pointer<number>(0);

    /* Local variables */

    const it = new Pointer<number>(0);
    const beta = new Pointer<number>(0);
    const imin = new Pointer<number>(0);
    const imax = new Pointer<number>(0);
    const lrnd = new Pointer<boolean>(false);
    let small: number;


    /*  -- LAPACK auxiliary routine (version 3.2) -- */
    /*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
    /*     November 2006 */

    /*     .. Scalar Arguments .. */
    /*     .. */

    /*  Purpose */
    /*  ======= */

    /*  SLAMCH determines single precision machine parameters. */

    /*  Arguments */
    /*  ========= */

    /*  CMACH   (input) CHARACTER*1 */
    /*          Specifies the value to be returned by SLAMCH: */
    /*          = 'E' or 'e',   SLAMCH := eps */
    /*          = 'S' or 's ,   SLAMCH := sfmin */
    /*          = 'B' or 'b',   SLAMCH := base */
    /*          = 'P' or 'p',   SLAMCH := eps*base */
    /*          = 'N' or 'n',   SLAMCH := t */
    /*          = 'R' or 'r',   SLAMCH := rnd */
    /*          = 'M' or 'm',   SLAMCH := emin */
    /*          = 'U' or 'u',   SLAMCH := rmin */
    /*          = 'L' or 'l',   SLAMCH := emax */
    /*          = 'O' or 'o',   SLAMCH := rmax */

    /*          where */

    /*          eps   = relative machine precision */
    /*          sfmin = safe minimum, such that 1/sfmin does not overflow */
    /*          base  = base of the machine */
    /*          prec  = eps*base */
    /*          t     = number of (base) digits in the mantissa */
    /*          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise */
    /*          emin  = minimum exponent before (gradual) underflow */
    /*          rmin  = underflow threshold - base**(emin-1) */
    /*          emax  = largest exponent before overflow */
    /*          rmax  = overflow threshold  - (base**emax)*(1-eps) */

    /* ===================================================================== */

    /*     .. Parameters .. */
    /*     .. */
    /*     .. Local Scalars .. */
    /*     .. */
    /*     .. External Functions .. */
    /*     .. */
    /*     .. External Subroutines .. */
    /*     .. */
    /*     .. Save statement .. */
    /*     .. */
    /*     .. Data statements .. */
    /*     .. */
    /*     .. Executable Statements .. */

    if (first) {
        slamc2(beta, it, lrnd, eps, imin, rmin, imax, rmax);
        base.assign(beta);
        t.assign(it);
        if (lrnd.value)
        {
            rnd.value = 1;
            i__1.value = 1 - it.value;
            eps.value = pow_ri(base, i__1) / 2;
        } else {
            rnd.value = 0;
            i__1.value = 1 - it.value;
            eps.value = pow_ri(base, i__1);
        }

        prec.value = Math.fround(eps.value * base.value);
        emin.assign(imin);
        emax.assign(imax);
        sfmin.assign(rmin);
        small = 1. / rmax.value;
        if (small >= sfmin.value)
        {

            /*           Use SMALL plus a bit, to avoid the possibility of rounding */
            /*           causing overflow when computing  1/sfmin. */

            sfmin.value = small * (eps.value + 1.);
        }
    }

    first = false;
    switch (cmach[0])
    {
        case "E": return eps.value;
        case "S": return sfmin.value;
        case "B": return base.value;
        case "P": return prec.value;
        case "N": return t.value;
        case "R": return rnd.value;
        case "M": return emin.value;
        case "U": return rmin.value;
        case "L": return emax.value;
        case "O": return rmax.value;
        default: throw new Error(`Unknown argument ${cmach}`);
    }

    /*     End of SLAMCH */

}; /* slamch_ */