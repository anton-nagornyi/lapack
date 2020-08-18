import {slassq} from "./slassq";
import {offset} from "../blas/offset";

const c__1 = 1;
const {min, max, abs} = Math;

export const slange = (norm: string, m: number, n: number, a: Float32Array, lda: number,
                       work: Float32Array): number =>
{
    let value = 0.;
    let r__1: number;


    /*  -- LAPACK auxiliary routine (version 3.2) -- */
    /*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
    /*     November 2006 */

    /*     .. Scalar Arguments .. */
    /*     .. */
    /*     .. Array Arguments .. */
    /*     .. */

    /*  Purpose */
    /*  ======= */

    /*  SLANGE  returns the value of the one norm,  or the Frobenius norm, or */
    /*  the  infinity norm,  or the  element of  largest absolute value  of a */
    /*  real matrix A. */

    /*  Description */
    /*  =========== */

    /*  SLANGE returns the value */

    /*     SLANGE = ( max(abs(A(i,j))), NORM = 'M' or 'm' */
    /*              ( */
    /*              ( norm1(A),         NORM = '1', 'O' or 'o' */
    /*              ( */
    /*              ( normI(A),         NORM = 'I' or 'i' */
    /*              ( */
    /*              ( normF(A),         NORM = 'F', 'f', 'E' or 'e' */

    /*  where  norm1  denotes the  one norm of a matrix (maximum column sum), */
    /*  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and */
    /*  normF  denotes the  Frobenius norm of a matrix (square root of sum of */
    /*  squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm. */

    /*  Arguments */
    /*  ========= */

    /*  NORM    (input) CHARACTER*1 */
    /*          Specifies the value to be returned in SLANGE as described */
    /*          above. */

    /*  M       (input) INTEGER */
    /*          The number of rows of the matrix A.  M >= 0.  When M = 0, */
    /*          SLANGE is set to zero. */

    /*  N       (input) INTEGER */
    /*          The number of columns of the matrix A.  N >= 0.  When N = 0, */
    /*          SLANGE is set to zero. */

    /*  A       (input) REAL array, dimension (LDA,N) */
    /*          The m by n matrix A. */

    /*  LDA     (input) INTEGER */
    /*          The leading dimension of the array A.  LDA >= max(M,1). */

    /*  WORK    (workspace) REAL array, dimension (MAX(1,LWORK)), */
    /*          where LWORK >= M when NORM = 'I'; otherwise, WORK is not */
    /*          referenced. */

    /* ===================================================================== */

    /*     .. Parameters .. */
    /*     .. */
    /*     .. Local Scalars .. */
    /*     .. */
    /*     .. External Subroutines .. */
    /*     .. */
    /*     .. External Functions .. */
    /*     .. */
    /*     .. Intrinsic Functions .. */
    /*     .. */
    /*     .. Executable Statements .. */

    /* Parameter adjustments */
    const a_dim1 = lda;
    const a_offset = -(1 + a_dim1);
    const work_offset = -1;

    /* Function Body */
    if (min(m, n) == 0)
    {
        value = 0.;
    } else
    {
        switch (norm)
        {
            case "M":
                /*        Find max(abs(A(i,j))). */

                value = 0.;
                for (let j = 1; j <= n; ++j)
                {
                    for (let i = 1; i <= m; ++i)
                    {
                        /* Computing MAX */
                        let r__1;
                        const r__2 = value;
                        const r__3 = (r__1 = a[a_offset + i + j * a_dim1], abs(r__1));
                        value = max(r__2, r__3);
                        /* L10: */
                    }
                    /* L20: */
                }
                break;
            case "1":
            case "O":
                /*        Find norm1(A). */

                value = 0.;

                for (let j = 1; j <= n; ++j)
                {
                    let sum = 0.;
                    for (let i = 1; i <= m; ++i)
                    {
                        sum += (r__1 = a[a_offset + i + j * a_dim1], abs(r__1));
                        /* L30: */
                    }
                    value = max(value, sum);
                    /* L40: */
                }
                break;
            case "I":
                /*        Find normI(A). */

                for (let i = 1; i <= m; ++i)
                {
                    work[work_offset + i] = 0.;
                    /* L50: */
                }

                for (let j = 1; j <= n; ++j)
                {
                    for (let i = 1; i <= m; ++i)
                    {
                        work[work_offset + i] += (r__1 = a[a_offset + i + j * a_dim1], abs(r__1));
                        /* L60: */
                    }
                    /* L70: */
                }
                value = 0.;
                for (let i = 1; i <= m; ++i)
                {
                    /* Computing MAX */
                    let r__2: number;
                    r__1 = value, r__2 = work[work_offset + i];
                    value = max(r__1, r__2);
                    /* L80: */
                }
                break;
            case "F":
            case "E":
                /*        Find normF(A). */

                const scale = 0.;
                const sum = 1.;

                for (let j = 1; j <= n; ++j)
                {
                    slassq(m, offset(a, a_offset + j * a_dim1 + 1), c__1, scale, sum);
                    /* L90: */
                }
                value = scale * Math.sqrt(sum);
                break;
        }
    }
    return value;

    /*     End of SLANGE */

}; /* slange_ */