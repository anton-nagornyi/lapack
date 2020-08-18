/* Table of constant values */

import {iparmq} from "./iparmq";

const c__1 = 1;
const c_b163 = 0.;
const c_b164 = 1.;
const c__0 = 0;

export const ilaenv = (ispec: number, name__: string, opts: string, n1: number,
                       n2: number, n3: number, n4: number): number =>
{
    /*  -- LAPACK auxiliary routine (version 3.2) -- */
    /*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd.. */
    /*     January 2007 */

    /*     .. Scalar Arguments .. */
    /*     .. */

    /*  Purpose */
    /*  ======= */

    /*  ILAENV is called from the LAPACK routines to choose problem-dependent */
    /*  parameters for the local environment.  See ISPEC for a description of */
    /*  the parameters. */

    /*  ILAENV returns an INTEGER */
    /*  if ILAENV >= 0: ILAENV returns the value of the parameter specified by ISPEC */
    /*  if ILAENV < 0:  if ILAENV = -k, the k-th argument had an illegal value. */

    /*  This version provides a set of parameters which should give good, */
    /*  but not optimal, performance on many of the currently available */
    /*  computers.  Users are encouraged to modify this subroutine to set */
    /*  the tuning parameters for their particular machine using the option */
    /*  and problem size information in the arguments. */

    /*  This routine will not function correctly if it is converted to all */
    /*  lower case.  Converting it to all upper case is allowed. */

    /*  Arguments */
    /*  ========= */

    /*  ISPEC   (input) INTEGER */
    /*          Specifies the parameter to be returned as the value of */
    /*          ILAENV. */
    /*          = 1: the optimal blocksize; if this value is 1, an unblocked */
    /*               algorithm will give the best performance. */
    /*          = 2: the minimum block size for which the block routine */
    /*               should be used; if the usable block size is less than */
    /*               this value, an unblocked routine should be used. */
    /*          = 3: the crossover point (in a block routine, for N less */
    /*               than this value, an unblocked routine should be used) */
    /*          = 4: the number of shifts, used in the nonsymmetric */
    /*               eigenvalue routines (DEPRECATED) */
    /*          = 5: the minimum column dimension for blocking to be used; */
    /*               rectangular blocks must have dimension at least k by m, */
    /*               where k is given by ILAENV(2,...) and m by ILAENV(5,...) */
    /*          = 6: the crossover point for the SVD (when reducing an m by n */
    /*               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds */
    /*               this value, a QR factorization is used first to reduce */
    /*               the matrix to a triangular form.) */
    /*          = 7: the number of processors */
    /*          = 8: the crossover point for the multishift QR method */
    /*               for nonsymmetric eigenvalue problems (DEPRECATED) */
    /*          = 9: maximum size of the subproblems at the bottom of the */
    /*               computation tree in the divide-and-conquer algorithm */
    /*               (used by xGELSD and xGESDD) */
    /*          =10: ieee NaN arithmetic can be trusted not to trap */
    /*          =11: infinity arithmetic can be trusted not to trap */
    /*          12 <= ISPEC <= 16: */
    /*               xHSEQR or one of its subroutines, */
    /*               see IPARMQ for detailed explanation */

    /*  NAME    (input) CHARACTER*(*) */
    /*          The name of the calling subroutine, in either upper case or */
    /*          lower case. */

    /*  OPTS    (input) CHARACTER*(*) */
    /*          The character options to the subroutine NAME, concatenated */
    /*          into a single character string.  For example, UPLO = 'U', */
    /*          TRANS = 'T', and DIAG = 'N' for a triangular routine would */
    /*          be specified as OPTS = 'UTN'. */

    /*  N1      (input) INTEGER */
    /*  N2      (input) INTEGER */
    /*  N3      (input) INTEGER */
    /*  N4      (input) INTEGER */
    /*          Problem dimensions for the subroutine NAME; these may not all */
    /*          be required. */

    /*  Further Details */
    /*  =============== */

    /*  The following conventions have been used when calling ILAENV from the */
    /*  LAPACK routines: */
    /*  1)  OPTS is a concatenation of all of the character options to */
    /*      subroutine NAME, in the same order that they appear in the */
    /*      argument list for NAME, even if they are not used in determining */
    /*      the value of the parameter specified by ISPEC. */
    /*  2)  The problem dimensions N1, N2, N3, N4 are specified in the order */
    /*      that they appear in the argument list for NAME.  N1 is used */
    /*      first, N2 second, and so on, and unused problem dimensions are */
    /*      passed a value of -1. */
    /*  3)  The parameter value returned by ILAENV is checked for validity in */
    /*      the calling subroutine.  For example, ILAENV is used to retrieve */
    /*      the optimal blocksize for STRTRI as follows: */

    /*      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 ) */
    /*      IF( NB.LE.1 ) NB = MAX( 1, N ) */

    /*  ===================================================================== */

    /*     .. Local Scalars .. */
    /*     .. */
    /*     .. Intrinsic Functions .. */
    /*     .. */
    /*     .. External Functions .. */
    /*     .. */
    /*     .. Executable Statements .. */

    switch (ispec)
    {
        case 1:
            return L10(name__, ispec, n2, n4);
        case 2:
            return L10(name__, ispec, n2, n4);
        case 3:
            return L10(name__, ispec, n2, n4);
        case 4:
            return L80();
        case 5:
            return L90();
        case 6:
            return L100(n1, n2);
        case 7:
            return L110();
        case 8:
            return L120();
        case 9:
            return L130();
        case 10:
            return L140();
        case 11:
            return L150();
        case 12:
            return L160(ispec, name__, opts, n1, n2, n3, n4);
        case 13:
            return L160(ispec, name__, opts, n1, n2, n3, n4);
        case 14:
            return L160(ispec, name__, opts, n1, n2, n3, n4);
        case 15:
            return L160(ispec, name__, opts, n1, n2, n3, n4);
        case 16:
            return L160(ispec, name__, opts, n1, n2, n3, n4);
        /*     Invalid value for ISPEC */
        default:
            return -1;
    }

    /*     End of ILAENV */

}; /* ilaenv_ */

const L10 = (name__: string, ispec: number, n2: number, n4: number): number =>
{
    /*     Convert NAME to upper case if the first character is lower case. */

    let ret_val = 1;
    const subnam = name__.toUpperCase();

    let c1 = subnam[0];
    let sname = c1[0] == 'S' || c1[0] == 'D';
    let cname = c1[0] == 'C' || c1[0] == 'Z';
    if (!(cname || sname))
    {
        return ret_val;
    }
    let c2 = subnam.substr(1);
    let c3 = subnam.substr(3);
    let c4 = subnam.substr(4);


    switch (ispec)
    {
        case 1:
            return L50(sname, cname, c2, c3, c4, n2, n4);
        case 2:
            return L60(sname, cname, c2, c3, c4);
        case 3:
            return L70(sname, cname, c2, c3, c4);
        /*     Invalid value for ISPEC */
        default:
            return -1;
    }
};
const L50 = (sname: boolean, cname: boolean, c2: string, c3: string, c4: string, n2: number, n4: number): number =>
{

    /*     ISPEC = 1:  block size */

    /*     In these examples, separate code is provided for setting NB for */
    /*     real and complex.  We assume that NB will take the same value in */
    /*     single or double precision. */

    let nb = 1;

    if (c2 === "GE")
    {
        if (c3 === "TRF")
        {
            if (sname)
            {
                nb = 64;
            } else
            {
                nb = 64;
            }
        } else if (c3 === "QRF" || c3 === "RQF" || c3 === "LQF" || c3 === "QLF")
        {
            if (sname)
            {
                nb = 32;
            } else
            {
                nb = 32;
            }
        } else if (c3 === "HRD")
        {
            if (sname)
            {
                nb = 32;
            } else
            {
                nb = 32;
            }
        } else if (c3 === "BRD")
        {
            if (sname)
            {
                nb = 32;
            } else
            {
                nb = 32;
            }
        } else if (c3 === "TRI")
        {
            if (sname)
            {
                nb = 64;
            } else
            {
                nb = 64;
            }
        }
    } else if (c2 === "PO")
    {
        if (c3 === "TRF")
        {
            if (sname)
            {
                nb = 64;
            } else
            {
                nb = 64;
            }
        }
    } else if (c2 === "SY")
    {
        if (c3 === "TRF")
        {
            if (sname)
            {
                nb = 64;
            } else
            {
                nb = 64;
            }
        } else if (sname && c3 === "TRD")
        {
            nb = 32;
        } else if (sname && c3 === "GST")
        {
            nb = 64;
        }
    } else if (cname && c2 === "HE")
    {
        if (c3 === "TRF")
        {
            nb = 64;
        } else if (c3 === "TRD")
        {
            nb = 32;
        } else if (c3 === "GST")
        {
            nb = 64;
        }
    } else if (sname && c2 === "OR")
    {
        if (c3.startsWith("G"))
        {
            if (c4 === "QR" || c4 === "RQ" || c4 === "LQ" || c4 === "QL" || c4 === "HR" || c4 === "TR" || c4 === "BR")
            {
                nb = 32;
            }
        } else if (c3.startsWith("M"))
        {
            if (c4 === "QR" || c4 === "RQ" || c4 === "LQ" || c4 === "QL" || c4 === "HR" || c4 === "TR" || c4 === "BR")
            {
                nb = 32;
            }
        }
    } else if (cname && c2 === "UN")
    {
        if (c3.startsWith("G"))
        {
            if (c4 === "QR" || c4 === "RQ" || c4 === "LQ" || c4 === "QL" || c4 === "HR" || c4 === "TR" || c4 === "BR")
            {
                nb = 32;
            }
        } else if (c3.startsWith("M"))
        {
            if (c4 === "QR" || c4 === "RQ" || c4 === "LQ" || c4 === "QL" || c4 === "HR" || c4 === "TR" || c4 === "BR")
            {
                nb = 32;
            }
        }
    } else if (c2 === "GB")
    {
        if (c3 === "TRF")
        {
            if (sname)
            {
                if (n4 <= 64)
                {
                    nb = 1;
                } else
                {
                    nb = 32;
                }
            } else
            {
                if (n4 <= 64)
                {
                    nb = 1;
                } else
                {
                    nb = 32;
                }
            }
        }
    } else if (c2 === "PB")
    {
        if (c3 === "TRF")
        {
            if (sname)
            {
                if (n2 <= 64)
                {
                    nb = 1;
                } else
                {
                    nb = 32;
                }
            } else
            {
                if (n2 <= 64)
                {
                    nb = 1;
                } else
                {
                    nb = 32;
                }
            }
        }
    } else if (c2 === "TR")
    {
        if (c3 === "TRI")
        {
            if (sname)
            {
                nb = 64;
            } else
            {
                nb = 64;
            }
        }
    } else if (c2 === "LA")
    {
        if (c3 === "UUM")
        {
            if (sname)
            {
                nb = 64;
            } else
            {
                nb = 64;
            }
        }
    } else if (sname && c2 === "ST")
    {
        if (c3 === "EBZ")
        {
            nb = 1;
        }
    }
    return nb;
};

const L60 = (sname: boolean, cname: boolean, c2: string, c3: string, c4: string) =>
{
    /*     ISPEC = 2:  minimum block size */

    let nbmin = 2;
    if (c2 === "GE")
    {
        if (c3 === "QRF" || c3 === "RQF" || c3 === "LQF" || c3 === "QLF")
        {
            if (sname)
            {
                nbmin = 2;
            } else
            {
                nbmin = 2;
            }
        } else if (c3 === "HRD")
        {
            if (sname)
            {
                nbmin = 2;
            } else
            {
                nbmin = 2;
            }
        } else if (c3 === "BRD")
        {
            if (sname)
            {
                nbmin = 2;
            } else
            {
                nbmin = 2;
            }
        } else if (c3 === "TRI")
        {
            if (sname)
            {
                nbmin = 2;
            } else
            {
                nbmin = 2;
            }
        }
    } else if (c2 === "SY")
    {
        if (c3 === "TRF")
        {
            if (sname)
            {
                nbmin = 8;
            } else
            {
                nbmin = 8;
            }
        } else if (sname && c3 === "TRD")
        {
            nbmin = 2;
        }
    } else if (cname && c2 === "HE")
    {
        if (c3 === "TRD")
        {
            nbmin = 2;
        }
    } else if (sname && c2 === "OR")
    {
        if (c3.startsWith("G"))
        {
            if (c4 === "QR" || c4 === "RQ" || c4 === "LQ" || c4 === "QL" || c4 === "HR" || c4 === "TR" || c4 === "BR")
            {
                nbmin = 2;
            }
        } else if (c3.startsWith("M"))
        {
            if (c4 === "QR" || c4 === "RQ" || c4 === "LQ" || c4 === "QL" || c4 === "HR" || c4 === "TR" || c4 === "BR")
            {
                nbmin = 2;
            }
        }
    } else if (cname && c2 === "UN")
    {
        if (c3.startsWith("G"))
        {
            if (c4 === "QR" || c4 === "RQ" || c4 === "LQ" || c4 === "QL" || c4 === "HR" || c4 === "TR" || c4 === "BR")
            {
                nbmin = 2;
            }
        } else if (c3.startsWith("M"))
        {
            if (c4 === "QR" || c4 === "RQ" || c4 === "LQ" || c4 === "QL" || c4 === "HR" || c4 === "TR" || c4 === "BR")
            {
                nbmin = 2;
            }
        }
    }
    return nbmin;
};

const L70 = (sname: boolean, cname: boolean, c2: string, c3: string, c4: string) =>
{
    /*     ISPEC = 3:  crossover point */

    let nx = 0;
    if (c2 === "GE")
    {
        if (c3 === "QRF" || c3 === "RQF" || c3 === "LQF" || c3 === "QLF")
        {
            if (sname)
            {
                nx = 128;
            } else
            {
                nx = 128;
            }
        } else if (c3 === "HRD")
        {
            if (sname)
            {
                nx = 128;
            } else
            {
                nx = 128;
            }
        } else if (c3 === "BRD")
        {
            if (sname)
            {
                nx = 128;
            } else
            {
                nx = 128;
            }
        }
    } else if (c2 === "SY")
    {
        if (sname && c3 === "TRD")
        {
            nx = 32;
        }
    } else if (cname && c2 === "HE")
    {
        if (c3 === "TRD")
        {
            nx = 32;
        }
    } else if (sname && c2 === "OR")
    {
        if (c3.startsWith("G"))
        {
            if (c4 === "QR" || c4 === "RQ" || c4 === "LQ" || c4 === "QL" || c4 === "HR" || c4 === "TR" || c4 === "BR")
            {
                nx = 128;
            }
        }
    } else if (cname && c2 === "UN")
    {
        if (c3.startsWith("G"))
        {
            if (c4 === "QR" || c4 === "RQ" || c4 === "LQ" || c4 === "QL" || c4 === "HR" || c4 === "TR" || c4 === "BR")
            {
                nx = 128;
            }
        }
    }
    return nx;
};
const L80 = (): number =>
{
    /*     ISPEC = 4:  number of shifts (used by xHSEQR) */
    return 6;
};

const L90 = (): number =>
{
    /*     ISPEC = 5:  minimum column dimension (not used) */
    return 2;
};

const L100 = (n1: number, n2: number): number =>
{
    /*     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD) */
    return Math.min(n1, n2) * 1.6;
};

const L110 = (): number =>
{
    /*     ISPEC = 7:  number of processors (not used) */
    return 1;
};

const L120 = (): number =>
{
    /*     ISPEC = 8:  crossover point for multishift (used by xHSEQR) */

    return 50;
};

const L130 = (): number =>
{
    /*     ISPEC = 9:  maximum size of the subproblems at the bottom of the */
    /*                 computation tree in the divide-and-conquer algorithm */
    /*                 (used by xGELSD and xGESDD) */

    return 25;
};

const L140 = (): number =>
{
    /*     ISPEC = 10: ieee NaN arithmetic can be trusted not to trap */
    return 1;
};

const L150 = (): number =>
{
    /*     ISPEC = 11: infinity arithmetic can be trusted not to trap */
    return 1;
};

const L160 = (ispec: number, name__: string, opts: string, n1: number, n2: number, n3: number, n4: number): number =>
{
    /*     12 <= ISPEC <= 16: xHSEQR or one of its subroutines. */

    return iparmq(ispec, name__, opts, n1, n2, n3, n4);
};