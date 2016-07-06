      SUBROUTINE DSLUCS(N, B, X, NELT, IA, JA, A, ISYM, ITOL, TOL,
     $     ITMAX, ITER, ERR, IERR, IUNIT, RWORK, LENW, IWORK, LENIW )
C *Arguments:
C N      :IN       Integer.
C         Order of the Matrix.
C B      :IN       Double Precision B(N).
C         Right-hand side vector.
C X      :INOUT    Double Precision X(N).
C         On input X is your initial guess for solution vector.
C         On output X is the final approximate solution.
C NELT   :IN       Integer.
C         Number of Non-Zeros stored in A.
C IA     :INOUT    Integer IA(NELT).
C JA     :INOUT    Integer JA(NELT).
C A      :INOUT    Double Precision A(NELT).
C         These arrays should hold the matrix A in either the SLAP
C         Triad format or the SLAP Column format.  See "Description", 
C         below.  If the SLAP Triad format is chosen it is changed 
C         internally to the SLAP Column format.
C ISYM   :IN       Integer.
C         Flag to indicate symmetric storage format.
C         If ISYM=0, all nonzero entries of the matrix are stored.
C         If ISYM=1, the matrix is symmetric, and only the upper
C         or lower triangle of the matrix is stored.
C ITOL   :IN       Integer.
C         Flag to indicate type of convergence criterion.
C         If ITOL=1, iteration stops when the 2-norm of the residual 
C         divided by the 2-norm of the right-hand side is less than TOL.
C         This routine must calculate the residual from R = A*X - B.
C         This is un-natural and hence expensive for this type of iter-
C         ative method.  ITOL=2 is *STRONGLY* recommended.
C         If ITOL=2, iteration stops when the 2-norm of M-inv times the 
C         residual divided by the 2-norm of M-inv times the right hand 
C         side is less than tol, where M-inv time a vector is the pre-
C         conditioning step.  This is the *NATURAL* stopping for this 
C         iterative method and is *STRONGLY* recommended.
C TOL    :IN       Double Precision.
C         Convergence criterion, as described above.
C ITMAX  :IN       Integer.
C         Maximum number of iterations.
C ITER   :OUT      Integer.
C         Number of iterations required to reach convergence, or 
C         ITMAX+1 if convergence criterion could not be achieved in 
C         ITMAX iterations.
C ERR    :OUT      Double Precision.
C         Error estimate of error in final approximate solution, as 
C         defined by ITOL.
C IERR   :OUT      Integer.
C         Return error flag.
C           IERR = 0 => All went well.
C           IERR = 1 => Insufficient storage allocated 
C                       for WORK or IWORK.
C           IERR = 2 => Method failed to converge in 
C                       ITMAX steps.
C           IERR = 3 => Error in user input.  Check input
C                       value of N, ITOL.
C           IERR = 4 => User error tolerance set too tight.
C                       Reset to 500.0*D1MACH(3).  Iteration proceeded.
C           IERR = 5 => Breakdown of the method detected.
C                       $(r0,r) approximately 0.0$.
C           IERR = 6 => Stagnation of the method detected.
C                        $(r0,v) approximately 0.0$.
C           IERR = 7 => Incomplete factorization broke down
C                       and was fudged.  Resulting preconditioning may
C                       be less than the best.
C IUNIT  :IN       Integer.
C         Unit number on which to write the error at each iteration, 
C         if this is desired for monitoring convergence.  If unit 
C         number is 0, no writing will occur.
C RWORK  :WORK     Double Precision RWORK(LENW).
C         Double Precision array used for workspace.  NEL is the 
C         number of non-
C         zeros in the lower triangle of the matrix (including the
C         diagonal).  NU is the number of nonzeros in the upper
C         triangle of the matrix (including the diagonal).
C LENW   :IN       Integer.
C         Length of the double precision workspace, RWORK.  
C         LENW >= NEL+NU+8*N.
C IWORK  :WORK     Integer IWORK(LENIW).
C         Integer array used for workspace.  NEL is the number of non-
C         zeros in the lower triangle of the matrix (including the
C         diagonal).  NU is the number of nonzeros in the upper
C         triangle of the matrix (including the diagonal).
C         Upon return the following locations of IWORK hold information
C         which may be of use to the user:
C         IWORK(9)  Amount of Integer workspace actually used.
C         IWORK(10) Amount of Double Precision workspace actually used.
C LENIW  :IN       Integer.
C         Length of the integer workspace, IWORK.  
C         LENIW >= NEL+NU+4*N+12.
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, ITOL, ITMAX, ITER
      INTEGER IERR, IUNIT, LENW, IWORK(LENIW), LENIW
      DOUBLE PRECISION B(N), X(N), A(NELT), TOL, ERR, RWORK(LENW)
      EXTERNAL DSMV1, DSLUI
      PARAMETER (LOCRB=1, LOCIB=11)
C
C         Change the SLAP input matrix IA, JA, A to SLAP-Column format.
C***FIRST EXECUTABLE STATEMENT  DSLUCS
      IERR = 0
      IF( N.LT.1 .OR. NELT.LT.1 ) THEN
         IERR = 3
         RETURN
      ENDIF
      CALL DS2Y( N, NELT, IA, JA, A, ISYM )
C
C         Count number of Non-Zero elements preconditioner ILU matrix.
C         Then set up the work arrays.
      NL = 0
      NU = 0
      DO 20 ICOL = 1, N
C         Don't count diagonal.
         JBGN = JA(ICOL)+1
         JEND = JA(ICOL+1)-1
         IF( JBGN.LE.JEND ) THEN
CVD$ NOVECTOR
            DO 10 J = JBGN, JEND
               IF( IA(J).GT.ICOL ) THEN
                  NL = NL + 1
                  IF( ISYM.NE.0 ) NU = NU + 1
               ELSE
                  NU = NU + 1
               ENDIF
 10         CONTINUE
         ENDIF
 20   CONTINUE
C         
      LOCIL = LOCIB
      LOCJL = LOCIL + N+1
      LOCIU = LOCJL + NL
      LOCJU = LOCIU + NU
      LOCNR = LOCJU + N+1
      LOCNC = LOCNR + N
      LOCIW = LOCNC + N
C
      LOCL   = LOCRB
      LOCDIN = LOCL + NL
      LOCUU  = LOCDIN + N
      LOCR   = LOCUU + NU
      LOCR0  = LOCR + N
      LOCP   = LOCR0 + N
      LOCQ   = LOCP + N
      LOCU   = LOCQ + N
      LOCV1  = LOCU + N
      LOCV2  = LOCV1 + N
      LOCW   = LOCV2 + N
C
C         Check the workspace allocations.
      CALL DCHKW( 'DSLUCS', LOCIW, LENIW, LOCW, LENW, IERR, ITER, ERR )
      IF( IERR.NE.0 ) RETURN
C
      IWORK(1) = LOCIL
      IWORK(2) = LOCJL
      IWORK(3) = LOCIU
      IWORK(4) = LOCJU
      IWORK(5) = LOCL
      IWORK(6) = LOCDIN
      IWORK(7) = LOCUU
      IWORK(9) = LOCIW
      IWORK(10) = LOCW
C
C         Compute the Incomplete LU decomposition.
      CALL DSILUS( N, NELT, IA, JA, A, ISYM, NL, IWORK(LOCIL),
     $     IWORK(LOCJL), RWORK(LOCL), RWORK(LOCDIN), NU, IWORK(LOCIU),
     $     IWORK(LOCJU), RWORK(LOCUU), IWORK(LOCNR), IWORK(LOCNC) )
C         
C         Perform the incomplete LU preconditioned 
C         BiConjugate Gradient Squared algorithm.
      CALL DCGS1(N, B, X, NELT, IA, JA, A, ISYM, DSMV1,
     $     DSLUI, ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT,
     $     RWORK(LOCR), RWORK(LOCR0), RWORK(LOCP),
     $     RWORK(LOCQ), RWORK(LOCU), RWORK(LOCV1),
     $     RWORK(LOCV2), RWORK, IWORK )
      RETURN
C------------- LAST LINE OF DSLUCS FOLLOWS ----------------------------
      END

      SUBROUTINE DCGS1(N, B, X, NELT, IA, JA, A, ISYM, MATVEC,
     $     MSOLVE, ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, 
     $     R, R0, P, Q, U, V1, V2, RWORK, IWORK)
C***BEGIN PROLOGUE  DCGS
C***DATE WRITTEN   890404   (YYMMDD)
C***REVISION DATE  890404   (YYMMDD)
C***CATEGORY NO.  D2A4
C***KEYWORDS  LIBRARY=SLATEC(SLAP),
C             TYPE=DOUBLE PRECISION(DCGS-D),
C             Non-Symmetric Linear system, Sparse, 
C             Iterative Precondition,  BiConjugate Gradient
C***AUTHOR  Greenbaum, Anne, Courant Institute
C           Seager, Mark K., (LLNL)
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-300
C             Livermore, CA 94550 (415) 423-3141
C             seager@lll-crg.llnl.gov
C***PURPOSE  Preconditioned BiConjugate Gradient Sparse Ax=b solver.
C            Routine to solve a Non-Symmetric linear system Ax = b 
C            using the Preconditioned BiConjugate Gradient method.
C***DESCRIPTION
C *Usage:
C      INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, ITOL, ITMAX
C      INTEGER ITER, IERR, IUNIT, IWORK(USER DEFINABLE)
C      DOUBLE PRECISION B(N), X(N), A(NELT), TOL, ERR, R(N), R0(N), P(N)
C      DOUBLE PRECISION Q(N), U(N), V1(N), V2(N), RWORK(USER DEFINABLE)
C      EXTERNAL MATVEC, MSOLVE
C
C      CALL DCGS(N, B, X, NELT, IA, JA, A, ISYM, MATVEC, 
C     $     MSOLVE, ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, 
C     $     R, R0, P, Q, U, V1, V2, RWORK, IWORK)
C
C *Arguments:
C N      :IN       Integer
C         Order of the Matrix.
C B      :IN       Double Precision B(N).
C         Right-hand side vector.
C X      :INOUT    Double Precision X(N).
C         On input X is your initial guess for solution vector.
C         On output X is the final approximate solution.
C NELT   :IN       Integer.
C         Number of Non-Zeros stored in A.
C IA     :IN       Integer IA(NELT).
C JA     :IN       Integer JA(NELT).
C A      :IN       Double Precision A(NELT).
C         These arrays contain the matrix data structure for A.
C         It could take any form.  See "Description", below
C         for more late breaking details...
C ISYM   :IN       Integer.
C         Flag to indicate symmetric storage format.
C         If ISYM=0, all nonzero entries of the matrix are stored.
C         If ISYM=1, the matrix is symmetric, and only the upper
C         or lower triangle of the matrix is stored.
C MATVEC :EXT      External.
C         Name of a routine which  performs the matrix vector multiply
C         operation  Y = A*X  given A and X.  The  name of  the MATVEC
C         routine must  be declared external  in the  calling program.
C         The calling sequence of MATVEC is:
C             CALL MATVEC( N, X, Y, NELT, IA, JA, A, ISYM )
C         Where N is the number of unknowns, Y is the product A*X upon
C         return,  X is an input  vector.  NELT, IA,  JA,  A and  ISYM
C         define the SLAP matrix data structure: see Description,below.
C MSOLVE :EXT      External.
C         Name of a routine which solves a linear system MZ = R  for Z
C         given R with the preconditioning matrix M (M is supplied via
C         RWORK  and IWORK arrays).   The name  of  the MSOLVE routine
C         must be declared  external  in the  calling   program.   The
C         calling sequence of MSLOVE is:
C             CALL MSOLVE(N, R, Z, NELT, IA, JA, A, ISYM, RWORK, IWORK)
C         Where N is the number of unknowns, R is  the right-hand side
C         vector, and Z is the solution upon return.  NELT,  IA, JA, A
C         and  ISYM define the SLAP  matrix  data structure: see  
C         Description, below.  RWORK is a  double precision array that 
C         can be used
C         to  pass   necessary  preconditioning     information and/or
C         workspace to MSOLVE.  IWORK is an integer work array for the
C         same purpose as RWORK.
C ITOL   :IN       Integer.
C         Flag to indicate type of convergence criterion.
C         If ITOL=1, iteration stops when the 2-norm of the residual 
C         divided by the 2-norm of the right-hand side is less than TOL.
C         This routine must calculate the residual from R = A*X - B.
C         This is un-natural and hence expensive for this type of iter-
C         ative method.  ITOL=2 is *STRONGLY* recommended.
C         If ITOL=2, iteration stops when the 2-norm of M-inv times the 
C         residual divided by the 2-norm of M-inv times the right hand 
C         side is less than tol, where M-inv time a vector is the pre-
C         conditioning step.  This is the *NATURAL* stopping for this 
C         iterative method and is *STRONGLY* recommended.
C         ITOL=11 is often useful for checking and comparing different 
C         routines.  For this case, the user must supply the "exact" 
C         solution or a very accurate approximation (one with an error 
C         much less than tol) through a common block,
C         COMMON /SOLBLK/ SOLN( )
C         if ITOL=11, iteration stops when the 2-norm of the difference 
C         between the iterative approximation and the user-supplied
C         solution divided by the 2-norm of the user-supplied solution 
C         is less than tol.
C TOL    :IN       Double Precision.
C         Convergence criterion, as described above.
C ITMAX  :IN       Integer.
C         Maximum number of iterations.
C ITER   :OUT      Integer.
C         Number of iterations required to reach convergence, or 
C         ITMAX+1 if convergence criterion could not be achieved in 
C         ITMAX iterations.
C ERR    :OUT      Double Precision.
C         Error estimate of error in final approximate solution, as 
C         defined by ITOL.
C IERR   :OUT      Integer.
C         Return error flag.
C           IERR = 0 => All went well.
C           IERR = 1 => Insufficient storage allocated 
C                       for WORK or IWORK.
C           IERR = 2 => Method failed to converge in 
C                       ITMAX steps.
C           IERR = 3 => Error in user input.  Check input
C                       value of N, ITOL.
C           IERR = 4 => User error tolerance set too tight.
C                       Reset to 500.0*D1MACH(3).  Iteration proceeded.
C           IERR = 5 => Breakdown of the method detected.
C                       $(r0,r) approximately 0.0$.
C           IERR = 6 => Stagnation of the method detected.
C                        $(r0,v) approximately 0.0$.
C IUNIT  :IN       Integer.
C         Unit number on which to write the error at each iteration, 
C         if this is desired for monitoring convergence.  If unit 
C         number is 0, no writing will occur.
C R      :WORK     Double Precision R(N).
C R0     :WORK     Double Precision R0(N).
C P      :WORK     Double Precision P(N).
C Q      :WORK     Double Precision Q(N).
C U      :WORK     Double Precision U(N).
C V1     :WORK     Double Precision V1(N).
C V2     :WORK     Double Precision V2(N).
C RWORK  :WORK     Double Precision RWORK(USER DEFINED).
C         Double Precision array that can be used for workspace in 
C         MSOLVE.
C IWORK  :WORK     Integer IWORK(USER DEFINED).
C         Integer array that can be used for workspace in MSOLVE.
C
C *Description
C       This routine does  not care  what matrix data   structure is
C       used for  A and M.  It simply   calls  the MATVEC and MSOLVE
C       routines, with  the arguments as  described above.  The user
C       could write any type of structure and the appropriate MATVEC
C       and MSOLVE routines.  It is assumed  that A is stored in the
C       IA, JA, A  arrays in some fashion and  that M (or INV(M)) is
C       stored  in  IWORK  and  RWORK   in  some fashion.   The SLAP
C       routines SDBCG and DSLUCS are examples of this procedure.
C       
C       Two  examples  of  matrix  data structures  are the: 1) SLAP
C       Triad  format and 2) SLAP Column format.
C       
C       =================== S L A P Triad format ===================
C
C       In  this   format only the  non-zeros are  stored.  They may
C       appear  in *ANY* order.   The user  supplies three arrays of
C       length NELT, where  NELT  is the number  of non-zeros in the
C       matrix:  (IA(NELT), JA(NELT),  A(NELT)).  For each  non-zero
C       the  user puts   the row  and  column index   of that matrix
C       element in the IA and JA arrays.  The  value of the non-zero
C       matrix  element is  placed in  the corresponding location of
C       the A  array.  This is  an extremely easy data  structure to
C       generate.  On  the other hand it  is  not too  efficient  on
C       vector  computers   for the  iterative  solution  of  linear
C       systems.  Hence, SLAP  changes this input  data structure to
C       the SLAP   Column  format for the  iteration (but   does not
C       change it back).
C       
C       Here is an example of the  SLAP Triad   storage format for a
C       5x5 Matrix.  Recall that the entries may appear in any order.
C
C           5x5 Matrix       SLAP Triad format for 5x5 matrix on left.
C                              1  2  3  4  5  6  7  8  9 10 11
C       |11 12  0  0 15|   A: 51 12 11 33 15 53 55 22 35 44 21
C       |21 22  0  0  0|  IA:  5  1  1  3  1  5  5  2  3  4  2
C       | 0  0 33  0 35|  JA:  1  2  1  3  5  3  5  2  5  4  1
C       | 0  0  0 44  0|
C       |51  0 53  0 55|
C
C       =================== S L A P Column format ==================
C       This routine requires  that the  matrix  A be  stored in the
C       SLAP Column format.  In this format the non-zeros are stored
C       counting down columns (except for  the diagonal entry, which
C       must appear first in each  "column")  and are stored  in the
C       double precision array A.   In other words,  for each column
C       in the matrix put the diagonal entry in  A.  Then put in the
C       other non-zero  elements going down  the column (except  the
C       diagonal) in order.   The  IA array holds the  row index for
C       each non-zero.  The JA array holds the offsets  into the IA,
C       A arrays  for  the  beginning  of each   column.   That  is,
C       IA(JA(ICOL)),  A(JA(ICOL)) points   to the beginning  of the
C       ICOL-th   column    in    IA and   A.      IA(JA(ICOL+1)-1),
C       A(JA(ICOL+1)-1) points to  the  end of the   ICOL-th column.
C       Note that we always have  JA(N+1) = NELT+1,  where N is  the
C       number of columns in  the matrix and NELT  is the number  of
C       non-zeros in the matrix.
C       
C       Here is an example of the  SLAP Column  storage format for a
C       5x5 Matrix (in the A and IA arrays '|'  denotes the end of a 
C       column):
C       
C       5x5 Matrix      SLAP Column format for 5x5 matrix on left.
C       1  2  3    4  5    6  7    8    9 10 11
C       |11 12  0  0 15|   A: 11 21 51 | 22 12 | 33 53 | 44 | 55 15 35
C       |21 22  0  0  0|  IA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3
C       | 0  0 33  0 35|  JA:  1  4  6    8  9   12
C       | 0  0  0 44  0|
C       |51  0 53  0 55|
C       
C *Precision:           Double Precision
C *See Also:
C       DSDCGS, DSLUCS
C***REFERENCES  1. P. Sonneveld, ``CGS, a fast Lanczos-type solver
C                 for nonsymmetric linear systems'', Delft University
C                 of Technology Report 84-16, Department of Math-
C                 ematics and Informatics, Julianalaan 132, 2628 BL
C                 Delft, Phone 015-784568.
C
C               2. E.F. Kaasschieter, ``The solution of non-symmetric
C                 linear systems by bi-conjugate gradients or conjugate
C                 gradients squared,''  Delft University of Tech-
C                 nology Report 86-21, Department of Mathematics and 
C                 Informatics, Julianalaan 132, 2628 BL Delft, 
C                 Phone 015-784568.
C***ROUTINES CALLED  MATVEC, MSOLVE, ISDCGS, DDOT, D1MACH
C***END PROLOGUE  DCGS
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, ITOL, ITMAX
      INTEGER ITER, IERR, IUNIT, IWORK(*)
      DOUBLE PRECISION B(N), X(N), A(NELT), TOL, ERR, R(N), R0(N), P(N)
      DOUBLE PRECISION Q(N), U(N), V1(N), V2(N), RWORK(*)
      EXTERNAL MATVEC, MSOLVE
C
C         Check some of the input data.
C***FIRST EXECUTABLE STATEMENT  DCGS
      ITER = 0
      IERR = 0
      IF( N.LT.1 ) THEN
         IERR = 3
         RETURN
      ENDIF
      TOLMIN = 500.0*D1MACH(3)
      IF( TOL.LT.TOLMIN ) THEN
         TOL = TOLMIN
         IERR = 4
      ENDIF
C         
C         Calculate initial residual and pseudo-residual, and check
C         stopping criterion.
      CALL MATVEC(N, X, R, NELT, IA, JA, A, ISYM)
      DO 10 I = 1, N
         V1(I)  = R(I) - B(I)
 10   CONTINUE
      CALL MSOLVE(N, V1, R, NELT, IA, JA, A, ISYM, RWORK, IWORK)
C         
      IF( ISDCGS(N, B, X, NELT, IA, JA, A, ISYM, MATVEC, MSOLVE, 
     $     ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, R, R0, P, Q, 
     $     U, V1, V2, RWORK, IWORK, AK, BK, BNRM, SOLNRM) .NE. 0 )
     $     GO TO 200
      IF( IERR.NE.0 ) RETURN
C
C         Set initial values.
C
      FUZZ = D1MACH(3)**2
      DO 20 I = 1, N
         R0(I) = R(I)
 20   CONTINUE
      RHONM1 = 1.0
C         
C         ***** ITERATION LOOP *****
C         
      DO 100 K=1,ITMAX
         ITER = K
C
C         Calculate coefficient BK and direction vectors U, V and P.
         RHON = DDOT(N, R0, 1, R, 1)
         IF( ABS(RHONM1).LT.FUZZ ) GOTO 998
         BK = RHON/RHONM1
         IF( ITER.EQ.1 ) THEN
            DO 30 I = 1, N
               U(I) = R(I)
               P(I) = R(I)
 30         CONTINUE
         ELSE
            DO 40 I = 1, N
               U(I) = R(I) + BK*Q(I)
               V1(I) = Q(I) + BK*P(I)
 40         CONTINUE
            DO 50 I = 1, N
               P(I) = U(I) + BK*V1(I)
 50         CONTINUE
         ENDIF
C         
C         Calculate coefficient AK, new iterate X, Q
         CALL MATVEC(N, P, V2, NELT, IA, JA, A, ISYM)
         CALL MSOLVE(N, V2, V1, NELT, IA, JA, A, ISYM, RWORK, IWORK)
         SIGMA = DDOT(N, R0, 1, V1, 1)
         IF( ABS(SIGMA).LT.FUZZ ) GOTO 999
         AK = RHON/SIGMA
         AKM = -AK
         DO 60 I = 1, N
            Q(I) = U(I) + AKM*V1(I)
 60      CONTINUE
         DO 70 I = 1, N
            V1(I) = U(I) + Q(I)
 70      CONTINUE
C         X = X - ak*V1.
         CALL DAXPY( N, AKM, V1, 1, X, 1 )
C                     -1
C         R = R - ak*M  *A*V1
         CALL MATVEC(N, V1, V2, NELT, IA, JA, A, ISYM)
         CALL MSOLVE(N, V2, V1, NELT, IA, JA, A, ISYM, RWORK, IWORK)
         CALL DAXPY( N, AKM, V1, 1, R, 1 )
C         
C         check stopping criterion.
         IF( ISDCGS(N, B, X, NELT, IA, JA, A, ISYM, MATVEC, MSOLVE, 
     $        ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, R, R0, P, Q, 
     $        U, V1, V2, RWORK, IWORK, AK, BK, BNRM, SOLNRM) .NE. 0 )
     $        GO TO 200
C
C         Update RHO.
         RHONM1 = RHON
 100  CONTINUE
C         
C         *****   end of loop  *****
C         Stopping criterion not satisfied.
      ITER = ITMAX + 1
      IERR = 2
 200  RETURN
C
C         Breakdown of method detected.
 998  IERR = 5
      RETURN
C
C         Stagnation of method detected.
 999  IERR = 6
      RETURN
C------------- LAST LINE OF DCGS FOLLOWS ----------------------------
      END

*DECK DSILUS
      SUBROUTINE DSILUS(N, NELT, IA, JA, A, ISYM, NL, IL, JL,
     $     L, DINV, NU, IU, JU, U, NROW, NCOL)
C***BEGIN PROLOGUE  DSILUS
C***DATE WRITTEN   890404   (YYMMDD)
C***REVISION DATE  890404   (YYMMDD)
C***CATEGORY NO.  D2A4, D2B4
C***KEYWORDS  LIBRARY=SLATEC(SLAP),
C             TYPE=DOUBLE PRECISION(DSILUS-D),
C             Non-Symmetric Linear system, Sparse, 
C             Iterative Precondition, Incomplete LU Factorization
C***AUTHOR  Greenbaum, Anne, Courant Institute
C           Seager, Mark K., (LLNL)
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-300
C             Livermore, CA 94550 (415) 423-3141
C             seager@lll-crg.llnl.gov
C***PURPOSE  Incomplete LU Decomposition Preconditioner SLAP Set Up.
C            Routine to generate the incomplete LDU decomposition of a 
C            matrix.  The  unit lower triangular factor L is stored by 
C            rows and the  unit upper triangular factor U is stored by 
C            columns.  The inverse of the diagonal matrix D is stored.
C            No fill in is allowed.
C***DESCRIPTION
C *Usage:
C     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM
C     INTEGER NL, IL(N+1), JL(NL), NU, IU(N+1), JU(NU)
C     INTEGER NROW(N), NCOL(N)
C     DOUBLE PRECISION A(NELT), L(NL), U(NU), DINV(N)
C
C     CALL DSILUS( N, NELT, IA, JA, A, ISYM, NL, IL, JL, L, 
C    $    DINV, NU, IU, JU, U, NROW, NCOL )
C
C *Arguments:
C N      :IN       Integer
C         Order of the Matrix.
C NELT   :IN       Integer.
C         Number of elements in arrays IA, JA, and A.
C IA     :IN       Integer IA(NELT).
C JA     :IN       Integer JA(NELT).
C A      :IN       Double Precision A(NELT).
C         These arrays should hold the matrix A in the SLAP Column
C         format.  See "Description", below. 
C ISYM   :IN       Integer.
C         Flag to indicate symmetric storage format.
C         If ISYM=0, all nonzero entries of the matrix are stored.
C         If ISYM=1, the matrix is symmetric, and only the lower 
C         triangle of the matrix is stored.
C NL     :OUT      Integer.
C         Number of non-zeros in the EL array.
C IL     :OUT      Integer IL(N+1).
C JL     :OUT      Integer JL(NL).
C L      :OUT      Double Precision L(NL).
C         IL, JL, L  contain the unit ower  triangular factor of  the
C         incomplete decomposition  of some  matrix stored  in   SLAP
C         Row format.     The   Diagonal  of ones  *IS*  stored.  See
C         "DESCRIPTION", below for more details about the SLAP format.
C NU     :OUT      Integer.
C         Number of non-zeros in the U array.     
C IU     :OUT      Integer IU(N+1).
C JU     :OUT      Integer JU(NU).
C U      :OUT      Double Precision     U(NU).
C         IU, JU, U contain   the unit upper triangular factor of the
C         incomplete  decomposition    of some matrix  stored in SLAP
C         Column  format.   The Diagonal of ones   *IS*  stored.  See 
C         "Description", below  for  more  details  about  the   SLAP 
C         format.
C NROW   :WORK     Integer NROW(N).
C         NROW(I) is the number of non-zero elements in the I-th row
C         of L.
C NCOL   :WORK     Integer NCOL(N).
C         NCOL(I) is the number of non-zero elements in the I-th 
C         column of U.
C
C *Description
C       IL, JL, L should contain the unit  lower triangular factor of
C       the incomplete decomposition of the A matrix  stored in SLAP
C       Row format.  IU, JU, U should contain  the unit upper factor
C       of the  incomplete decomposition of  the A matrix  stored in
C       SLAP Column format This ILU factorization can be computed by
C       the DSILUS routine.  The diagonals (which is all one's) are
C       stored.
C
C       =================== S L A P Column format ==================
C       This routine requires  that the  matrix  A be  stored in the
C       SLAP Column format.  In this format the non-zeros are stored
C       counting down columns (except for  the diagonal entry, which
C       must appear first in each  "column")  and are stored  in the
C       double precision array A.   In other words,  for each column
C       in the matrix put the diagonal entry in  A.  Then put in the
C       other non-zero  elements going down  the column (except  the
C       diagonal) in order.   The  IA array holds the  row index for
C       each non-zero.  The JA array holds the offsets  into the IA,
C       A arrays  for  the  beginning  of each   column.   That  is,
C       IA(JA(ICOL)),  A(JA(ICOL)) points   to the beginning  of the
C       ICOL-th   column    in    IA and   A.      IA(JA(ICOL+1)-1),
C       A(JA(ICOL+1)-1) points to  the  end of the   ICOL-th column.
C       Note that we always have  JA(N+1) = NELT+1,  where N is  the
C       number of columns in  the matrix and NELT  is the number  of
C       non-zeros in the matrix.
C       
C       Here is an example of the  SLAP Column  storage format for a
C       5x5 Matrix (in the A and IA arrays '|'  denotes the end of a 
C       column):
C       
C       5x5 Matrix      SLAP Column format for 5x5 matrix on left.
C       1  2  3    4  5    6  7    8    9 10 11
C       |11 12  0  0 15|   A: 11 21 51 | 22 12 | 33 53 | 44 | 55 15 35
C       |21 22  0  0  0|  IA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3
C       | 0  0 33  0 35|  JA:  1  4  6    8  9   12
C       | 0  0  0 44  0|
C       |51  0 53  0 55|
C       
C       ==================== S L A P Row format ====================
C       This routine requires  that the matrix A  be  stored  in the
C       SLAP  Row format.   In this format  the non-zeros are stored
C       counting across  rows (except for the diagonal  entry, which
C       must appear first in each "row") and  are stored in the 
C       double precision
C       array A.  In other words, for each row in the matrix put the
C       diagonal entry in  A.   Then   put  in the   other  non-zero
C       elements   going  across the  row (except   the diagonal) in
C       order.   The  JA array  holds   the column   index for  each
C       non-zero.   The IA  array holds the  offsets into  the JA, A
C       arrays  for   the   beginning  of   each  row.   That    is,
C       JA(IA(IROW)),  A(IA(IROW)) points  to  the beginning  of the
C       IROW-th row in JA and A.   JA(IA(IROW+1)-1), A(IA(IROW+1)-1)
C       points to the  end of the  IROW-th row.  Note that we always
C       have IA(N+1) =  NELT+1, where  N  is  the number of rows  in
C       the matrix  and NELT  is the  number   of  non-zeros in  the
C       matrix.
C       
C       Here is an example of the SLAP Row storage format for a  5x5
C       Matrix (in the A and JA arrays '|' denotes the end of a row):
C
C           5x5 Matrix         SLAP Row format for 5x5 matrix on left.
C                              1  2  3    4  5    6  7    8    9 10 11
C       |11 12  0  0 15|   A: 11 12 15 | 22 21 | 33 35 | 44 | 55 51 53
C       |21 22  0  0  0|  JA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3
C       | 0  0 33  0 35|  IA:  1  4  6    8  9   12
C       | 0  0  0 44  0|  
C       |51  0 53  0 55|  
C
C *Precision:           Double Precision
C *See Also:
C       SILUR
C***REFERENCES  1. Gene Golub & Charles Van Loan, "Matrix Computations",
C                 John Hopkins University Press; 3 (1983) IBSN 
C                 0-8018-3010-9.
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  DSILUS
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, NL, IL(N+1), JL(NL)
      INTEGER NU, IU(NU), JU(N+1), NROW(N), NCOL(N)
      DOUBLE PRECISION A(NELT), L(NL), DINV(N), U(NU)
C         
C         Count number of elements in each row of the lower triangle.
C***FIRST EXECUTABLE STATEMENT  DSILUS
      DO 10 I=1,N
         NROW(I) = 0
         NCOL(I) = 0
 10   CONTINUE
CVD$R NOCONCUR
CVD$R NOVECTOR
      DO 30 ICOL = 1, N
         JBGN = JA(ICOL)+1
         JEND = JA(ICOL+1)-1
         IF( JBGN.LE.JEND ) THEN
            DO 20 J = JBGN, JEND
               IF( IA(J).LT.ICOL ) THEN
                  NCOL(ICOL) = NCOL(ICOL) + 1
               ELSE
                  NROW(IA(J)) = NROW(IA(J)) + 1
                  IF( ISYM.NE.0 ) NCOL(IA(J)) = NCOL(IA(J)) + 1
               ENDIF
 20         CONTINUE
         ENDIF
 30   CONTINUE
      JU(1) = 1
      IL(1) = 1
      DO 40 ICOL = 1, N
         IL(ICOL+1) = IL(ICOL) + NROW(ICOL)
         JU(ICOL+1) = JU(ICOL) + NCOL(ICOL)
         NROW(ICOL) = IL(ICOL)
         NCOL(ICOL) = JU(ICOL)
 40   CONTINUE
C         
C         Copy the matrix A into the L and U structures.
      DO 60 ICOL = 1, N
         DINV(ICOL) = A(JA(ICOL))
         JBGN = JA(ICOL)+1
         JEND = JA(ICOL+1)-1
         IF( JBGN.LE.JEND ) THEN
            DO 50 J = JBGN, JEND
               IROW = IA(J)
               IF( IROW.LT.ICOL ) THEN
C         Part of the upper triangle.
                  IU(NCOL(ICOL)) = IROW
                  U(NCOL(ICOL)) = A(J)
                  NCOL(ICOL) = NCOL(ICOL) + 1
               ELSE
C         Part of the lower triangle (stored by row).
                  JL(NROW(IROW)) = ICOL
                  L(NROW(IROW)) = A(J)
                  NROW(IROW) = NROW(IROW) + 1
                  IF( ISYM.NE.0 ) THEN
C         Symmetric...Copy lower triangle into upper triangle as well.
                     IU(NCOL(IROW)) = ICOL
                     U(NCOL(IROW)) = A(J)
                     NCOL(IROW) = NCOL(IROW) + 1
                  ENDIF
               ENDIF
 50         CONTINUE
         ENDIF
 60   CONTINUE
C
C         Sort the rows of L and the columns of U.
      DO 110 K = 2, N
         JBGN = JU(K)
         JEND = JU(K+1)-1
         IF( JBGN.LT.JEND ) THEN
            DO 80 J = JBGN, JEND-1
               DO 70 I = J+1, JEND
                  IF( IU(J).GT.IU(I) ) THEN
                     ITEMP = IU(J)
                     IU(J) = IU(I)
                     IU(I) = ITEMP
                     TEMP = U(J)
                     U(J) = U(I)
                     U(I) = TEMP
                  ENDIF
 70            CONTINUE
 80         CONTINUE
         ENDIF
         IBGN = IL(K)
         IEND = IL(K+1)-1
         IF( IBGN.LT.IEND ) THEN
            DO 100 I = IBGN, IEND-1
               DO 90 J = I+1, IEND
                  IF( JL(I).GT.JL(J) ) THEN
                     JTEMP = JU(I)
                     JU(I) = JU(J)
                     JU(J) = JTEMP
                     TEMP = L(I)
                     L(I) = L(J)
                     L(J) = TEMP
                  ENDIF
 90            CONTINUE
 100        CONTINUE
         ENDIF
 110  CONTINUE
C
C         Perform the incomplete LDU decomposition.
      DO 300 I=2,N
C         
C           I-th row of L
         INDX1 = IL(I)
         INDX2 = IL(I+1) - 1
         IF(INDX1 .GT. INDX2) GO TO 200
         DO 190 INDX=INDX1,INDX2
            IF(INDX .EQ. INDX1) GO TO 180
            INDXR1 = INDX1
            INDXR2 = INDX - 1
            INDXC1 = JU(JL(INDX))
            INDXC2 = JU(JL(INDX)+1) - 1
            IF(INDXC1 .GT. INDXC2) GO TO 180
 160        KR = JL(INDXR1)
 170        KC = IU(INDXC1)
            IF(KR .GT. KC) THEN
               INDXC1 = INDXC1 + 1
               IF(INDXC1 .LE. INDXC2) GO TO 170
            ELSEIF(KR .LT. KC) THEN
               INDXR1 = INDXR1 + 1
               IF(INDXR1 .LE. INDXR2) GO TO 160
            ELSEIF(KR .EQ. KC) THEN
               L(INDX) = L(INDX) - L(INDXR1)*DINV(KC)*U(INDXC1)
               INDXR1 = INDXR1 + 1
               INDXC1 = INDXC1 + 1
               IF(INDXR1 .LE. INDXR2 .AND. INDXC1 .LE. INDXC2) GO TO 160
            ENDIF
 180        L(INDX) = L(INDX)/DINV(JL(INDX))
 190     CONTINUE
C         
C         ith column of u
 200     INDX1 = JU(I)
         INDX2 = JU(I+1) - 1
         IF(INDX1 .GT. INDX2) GO TO 260
         DO 250 INDX=INDX1,INDX2
            IF(INDX .EQ. INDX1) GO TO 240
            INDXC1 = INDX1
            INDXC2 = INDX - 1
            INDXR1 = IL(IU(INDX))
            INDXR2 = IL(IU(INDX)+1) - 1
            IF(INDXR1 .GT. INDXR2) GO TO 240
 210        KR = JL(INDXR1)
 220        KC = IU(INDXC1)
            IF(KR .GT. KC) THEN
               INDXC1 = INDXC1 + 1
               IF(INDXC1 .LE. INDXC2) GO TO 220
            ELSEIF(KR .LT. KC) THEN
               INDXR1 = INDXR1 + 1
               IF(INDXR1 .LE. INDXR2) GO TO 210
            ELSEIF(KR .EQ. KC) THEN
               U(INDX) = U(INDX) - L(INDXR1)*DINV(KC)*U(INDXC1)
               INDXR1 = INDXR1 + 1
               INDXC1 = INDXC1 + 1
               IF(INDXR1 .LE. INDXR2 .AND. INDXC1 .LE. INDXC2) GO TO 210
            ENDIF
 240        U(INDX) = U(INDX)/DINV(IU(INDX))
 250     CONTINUE
C         
C         ith diagonal element
 260     INDXR1 = IL(I)
         INDXR2 = IL(I+1) - 1
         IF(INDXR1 .GT. INDXR2) GO TO 300
         INDXC1 = JU(I)
         INDXC2 = JU(I+1) - 1
         IF(INDXC1 .GT. INDXC2) GO TO 300
 270     KR = JL(INDXR1)
 280     KC = IU(INDXC1)
         IF(KR .GT. KC) THEN
            INDXC1 = INDXC1 + 1
            IF(INDXC1 .LE. INDXC2) GO TO 280
         ELSEIF(KR .LT. KC) THEN
            INDXR1 = INDXR1 + 1
            IF(INDXR1 .LE. INDXR2) GO TO 270
         ELSEIF(KR .EQ. KC) THEN
            DINV(I) = DINV(I) - L(INDXR1)*DINV(KC)*U(INDXC1)
            INDXR1 = INDXR1 + 1
            INDXC1 = INDXC1 + 1
            IF(INDXR1 .LE. INDXR2 .AND. INDXC1 .LE. INDXC2) GO TO 270
         ENDIF
C         
 300  CONTINUE
C         
C         replace diagonal lts by their inverses.
CVD$ VECTOR
      DO 430 I=1,N
         DINV(I) = 1./DINV(I)
 430  CONTINUE
C         
      RETURN
C------------- LAST LINE OF DSILUS FOLLOWS ----------------------------
      END

*DECK DCHKW
      SUBROUTINE DCHKW( NAME, LOCIW, LENIW, LOCW, LENW,
     $     IERR, ITER, ERR )
C***BEGIN PROLOGUE  DCHKW
C***DATE WRITTEN   880225   (YYMMDD)
C***REVISION DATE  881213   (YYMMDD)
C***CATEGORY NO.  R2
C***KEYWORDS  LIBRARY=SLATEC(SLAP),
C             TYPE=DOUBLE PRECISION(DCHKW-D),
C             SLAP, Error Checking, Workspace Checking
C***AUTHOR  Seager, Mark K., (LLNL)
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-300
C             Livermore, CA 94550 (415) 423-3141
C             seager@lll-crg.llnl.gov
C***PURPOSE  SLAP WORK/IWORK Array Bounds Checker.
C            This routine checks the work array lengths  and  inter-
C            faces to the SLATEC  error  handler  if  a  problem  is 
C            found.
C***DESCRIPTION
C *Usage:
C     CHARACTER*(*) NAME
C     INTEGER LOCIW, LENIW, LOCW, LENW, IERR, ITER
C     DOUBLE PRECISION ERR
C
C     CALL DCHKW( NAME, LOCIW, LENIW, LOCW, LENW, IERR, ITER, ERR )
C
C *Arguments:
C NAME   :IN       Character*(*).
C         Name of the calling routine.  This is used in the output
C         message, if an error is detected.
C LOCIW  :IN       Integer.
C         Location of the first free element in the integer workspace
C         array.
C LENIW  :IN       Integer.
C         Length of the integer workspace array.
C LOCW   :IN       Integer.
C         Location of the first free element in the double precision 
C         workspace array.
C LENRW  :IN       Integer.
C         Length of the double precision workspace array.
C IERR   :OUT      Integer.
C         Return error flag.
C               IERR = 0 => All went well.
C               IERR = 1 => Insufficient storage allocated for 
C                           WORK or IWORK.
C ITER   :OUT      Integer.
C         Set to 0 if an error is detected.
C ERR    :OUT      Double Precision.
C         Set to a very large number if an error is detected.
C
C *Precision:           Double Precision
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH, XERRWV
C***END PROLOGUE  DCHKW
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      CHARACTER*(*) NAME
      CHARACTER*72 MESG
      INTEGER LOCIW, LENIW, LOCW, LENW, IERR, ITER
      DOUBLE PRECISION ERR, D1MACH
      EXTERNAL D1MACH, XERRWV
C
C         Check the Integer workspace situation.
C***FIRST EXECUTABLE STATEMENT  DCHKW
      IERR = 0
      IF( LOCIW.GT.LENIW ) THEN
         IERR = 1
         ITER = 0
         ERR = D1MACH(2)
         MESG = NAME // ': INTEGER work array too short. '//
     $        ' IWORK needs i1: have allocated i2.'
         CALL XERRWV( MESG, LEN(MESG), 1, 1, 2, LOCIW, LENIW,
     $        0, 0.0, 0.0 )
      ENDIF
C
C         Check the Double Precision workspace situation.
      IF( LOCW.GT.LENW ) THEN
         IERR = 1
         ITER = 0
         ERR = D1MACH(2)
         MESG = NAME // ': DOUBLE PRECISION work array too short. '//
     $        ' RWORK needs i1: have allocated i2.'
         CALL XERRWV( MESG, LEN(MESG), 1, 1, 2, LOCW, LENW,
     $        0, 0.0, 0.0 )
      ENDIF
      RETURN
C------------- LAST LINE OF DCHKW FOLLOWS ----------------------------
      END
*DECK DS2Y
      SUBROUTINE DS2Y(N, NELT, IA, JA, A, ISYM )
C***BEGIN PROLOGUE  DS2Y
C***DATE WRITTEN   871119   (YYMMDD)
C***REVISION DATE  881213   (YYMMDD)
C***CATEGORY NO.  D2A4, D2B4
C***KEYWORDS  LIBRARY=SLATEC(SLAP),
C             TYPE=DOUBLE PRECISION(DS2Y-D),
C             Linear system, SLAP Sparse
C***AUTHOR  Seager, Mark K., (LLNL)
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-300
C             Livermore, CA 94550 (415) 423-3141
C             seager@lll-crg.llnl.gov
C***PURPOSE  SLAP Triad to SLAP Column Format Converter.
C            Routine to convert from the SLAP Triad to SLAP Column
C            format.
C***DESCRIPTION
C *Usage:
C     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM
C     DOUBLE PRECISION A(NELT)
C
C     CALL DS2Y( N, NELT, IA, JA, A, ISYM )
C
C *Arguments:
C N      :IN       Integer
C         Order of the Matrix.
C NELT   :IN       Integer.
C         Number of non-zeros stored in A.
C IA     :INOUT    Integer IA(NELT).
C JA     :INOUT    Integer JA(NELT).
C A      :INOUT    Double Precision A(NELT).
C         These arrays should hold the matrix A in either the SLAP
C         Triad format or the SLAP Column format.  See "LONG 
C         DESCRIPTION", below.  If the SLAP Triad format is used
C         this format is translated to the SLAP Column format by
C         this routine.
C ISYM   :IN       Integer.
C         Flag to indicate symmetric storage format.
C         If ISYM=0, all nonzero entries of the matrix are stored.
C         If ISYM=1, the matrix is symmetric, and only the lower
C         triangle of the matrix is stored.
C
C *Precision:           Double Precision
C
C***LONG DESCRIPTION
C       The Sparse Linear Algebra Package (SLAP) utilizes two matrix
C       data structures: 1) the  SLAP Triad  format or  2)  the SLAP
C       Column format.  The user can hand this routine either of the
C       of these data structures.  If the SLAP Triad format is give
C       as input then this routine transforms it into SLAP Column
C       format.  The way this routine tells which format is given as
C       input is to look at JA(N+1).  If JA(N+1) = NELT+1 then we
C       have the SLAP Column format.  If that equality does not hold
C       then it is assumed that the IA, JA, A arrays contain the 
C       SLAP Triad format.
C       
C       =================== S L A P Triad format ===================
C       This routine requires that the  matrix A be   stored in  the
C       SLAP  Triad format.  In  this format only the non-zeros  are
C       stored.  They may appear in  *ANY* order.  The user supplies
C       three arrays of  length NELT, where  NELT is  the number  of
C       non-zeros in the matrix: (IA(NELT), JA(NELT), A(NELT)).  For
C       each non-zero the user puts the row and column index of that
C       matrix element  in the IA and  JA arrays.  The  value of the
C       non-zero   matrix  element is  placed  in  the corresponding
C       location of the A array.   This is  an  extremely  easy data
C       structure to generate.  On  the  other hand it   is  not too
C       efficient on vector computers for  the iterative solution of
C       linear systems.  Hence,   SLAP changes   this  input    data
C       structure to the SLAP Column format  for  the iteration (but
C       does not change it back).
C       
C       Here is an example of the  SLAP Triad   storage format for a
C       5x5 Matrix.  Recall that the entries may appear in any order.
C
C           5x5 Matrix       SLAP Triad format for 5x5 matrix on left.
C                              1  2  3  4  5  6  7  8  9 10 11
C       |11 12  0  0 15|   A: 51 12 11 33 15 53 55 22 35 44 21
C       |21 22  0  0  0|  IA:  5  1  1  3  1  5  5  2  3  4  2
C       | 0  0 33  0 35|  JA:  1  2  1  3  5  3  5  2  5  4  1
C       | 0  0  0 44  0|
C       |51  0 53  0 55|
C       
C       =================== S L A P Column format ==================
C       This routine requires  that the  matrix  A be  stored in the
C       SLAP Column format.  In this format the non-zeros are stored
C       counting down columns (except for  the diagonal entry, which
C       must appear first in each  "column")  and are stored  in the
C       double precision array A.   In other words,  for each column
C       in the matrix put the diagonal entry in  A.  Then put in the
C       other non-zero  elements going down  the column (except  the
C       diagonal) in order.   The  IA array holds the  row index for
C       each non-zero.  The JA array holds the offsets  into the IA,
C       A arrays  for  the  beginning  of each   column.   That  is,
C       IA(JA(ICOL)),  A(JA(ICOL)) points   to the beginning  of the
C       ICOL-th   column    in    IA and   A.      IA(JA(ICOL+1)-1),
C       A(JA(ICOL+1)-1) points to  the  end of the   ICOL-th column.
C       Note that we always have  JA(N+1) = NELT+1,  where N is  the
C       number of columns in  the matrix and NELT  is the number  of
C       non-zeros in the matrix.
C       
C       Here is an example of the  SLAP Column  storage format for a
C       5x5 Matrix (in the A and IA arrays '|'  denotes the end of a 
C       column):
C       
C       5x5 Matrix      SLAP Column format for 5x5 matrix on left.
C       1  2  3    4  5    6  7    8    9 10 11
C       |11 12  0  0 15|   A: 11 21 51 | 22 12 | 33 53 | 44 | 55 15 35
C       |21 22  0  0  0|  IA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3
C       | 0  0 33  0 35|  JA:  1  4  6    8  9   12
C       | 0  0  0 44  0|
C       |51  0 53  0 55|
C       
C***REFERENCES  (NONE)
C***ROUTINES CALLED  QS2I1D
C***END PROLOGUE  DS2Y
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER N, NELT, IA(NELT), JA(NELT), ISYM
      DOUBLE PRECISION A(NELT)
C
C         Check to see if the (IA,JA,A) arrays are in SLAP Column 
C         format.  If it's not then transform from SLAP Triad.
C***FIRST EXECUTABLE STATEMENT  DS2LT
      IF( JA(N+1).EQ.NELT+1 ) RETURN
C
C         Sort into ascending order by COLUMN (on the ja array).
C         This will line up the columns.
C
      CALL QS2I1D( JA, IA, A, NELT, 1 )
C         
C         Loop over each column to see where the column indicies change 
C         in the column index array ja.  This marks the beginning of the
C         next column.
C         
CVD$R NOVECTOR
      JA(1) = 1
      DO 20 ICOL = 1, N-1
         DO 10 J = JA(ICOL)+1, NELT
            IF( JA(J).NE.ICOL ) THEN
               JA(ICOL+1) = J
               GOTO 20
            ENDIF
 10      CONTINUE
 20   CONTINUE
      JA(N+1) = NELT+1
C         
C         Mark the n+2 element so that future calls to a SLAP routine 
C         utilizing the YSMP-Column storage format will be able to tell.
C         
      JA(N+2) = 0
C
C         Now loop thru the ia(i) array making sure that the Diagonal
C         matrix element appears first in the column.  Then sort the
C         rest of the column in ascending order.
C
      DO 70 ICOL = 1, N
         IBGN = JA(ICOL)
         IEND = JA(ICOL+1)-1
         DO 30 I = IBGN, IEND
            IF( IA(I).EQ.ICOL ) THEN
C         Swap the diag element with the first element in the column.
               ITEMP = IA(I)
               IA(I) = IA(IBGN)
               IA(IBGN) = ITEMP
               TEMP = A(I)
               A(I) = A(IBGN)
               A(IBGN) = TEMP
               GOTO 40
            ENDIF
 30      CONTINUE
 40      IBGN = IBGN + 1
         IF( IBGN.LT.IEND ) THEN
            DO 60 I = IBGN, IEND
               DO 50 J = I+1, IEND
                  IF( IA(I).GT.IA(J) ) THEN
                     ITEMP = IA(I)
                     IA(I) = IA(J)
                     IA(J) = ITEMP
                     TEMP = A(I)
                     A(I) = A(J)
                     A(J) = TEMP
                  ENDIF
 50            CONTINUE
 60         CONTINUE
         ENDIF
 70   CONTINUE
      RETURN
C------------- LAST LINE OF DS2Y FOLLOWS ----------------------------
      END
*DECK QS2I1D
      SUBROUTINE QS2I1D( IA, JA, A, N, KFLAG )
C***BEGIN PROLOGUE  QS2I1D
C***DATE WRITTEN   761118   (YYMMDD)
C***REVISION DATE  890125   (YYMMDD)
C***CATEGORY NO.  N6A2A
C***KEYWORDS  LIBRARY=SLATEC(SLAP),
C             TYPE=INTEGER(QS2I1D-I),
C             QUICKSORT,DOUBLETON QUICKSORT,SORT,SORTING
C***AUTHOR  Jones, R. E., (SNLA)
C           Kahaner, D. K., (NBS)
C           Seager, M. K., (LLNL) seager@lll-crg.llnl.gov
C           Wisniewski, J. A., (SNLA)
C***PURPOSE  Sort an integer array also moving an integer and DP array
C            This routine sorts the integer  array  IA and makes the
C            same interchanges   in the integer   array  JA  and the
C            double precision array A.  The  array IA may be  sorted
C            in increasing order or decreas- ing  order.  A slightly
C            modified QUICKSORT algorithm is used.
C
C***DESCRIPTION
C     Written by Rondall E Jones
C     Modified by John A. Wisniewski to use the Singleton QUICKSORT
C     algorithm. date 18 November 1976.
C
C     Further modified by David K. Kahaner
C     National Bureau of Standards
C     August, 1981
C
C     Even further modification made to bring the code up to the 
C     Fortran 77 level and make it more readable and to carry
C     along one integer array and one double precision array during 
C     the sort by
C     Mark K. Seager
C     Lawrence Livermore National Laboratory
C     November, 1987
C     This routine was adapted from the ISORT routine.
C
C     ABSTRACT
C         This routine sorts an integer array IA and makes the same
C         interchanges in the integer array JA and the double precision
C          array A.  
C         The array a may be sorted in increasing order or decreasing 
C         order.  A slightly modified quicksort algorithm is used.
C
C     DESCRIPTION OF PARAMETERS
C        IA - Integer array of values to be sorted.
C        JA - Integer array to be carried along.
C         A - Double Precision array to be carried along.
C         N - Number of values in integer array IA to be sorted.
C     KFLAG - Control parameter
C           = 1 means sort IA in INCREASING order.
C           =-1 means sort IA in DECREASING order.
C
C***REFERENCES
C     Singleton, R. C., Algorithm 347, "An Efficient Algorithm for 
C     Sorting with Minimal Storage", cacm, Vol. 12, No. 3, 1969, 
C     Pp. 185-187.
C***ROUTINES CALLED  XERROR
C***END PROLOGUE  QS2I1D
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
CVD$R NOVECTOR
CVD$R NOCONCUR
      DIMENSION IL(21),IU(21)
      INTEGER   IA(N),JA(N),IT,IIT,JT,JJT
      DOUBLE PRECISION A(N), TA, TTA
C
C***FIRST EXECUTABLE STATEMENT  QS2I1D
      NN = N
      IF (NN.LT.1) THEN
         CALL XERROR ( 'QS2I1D- the number of values to be sorted was no
     $T POSITIVE.',59,1,1)
         RETURN
      ENDIF
      IF( N.EQ.1 ) RETURN
      KK = IABS(KFLAG)
      IF ( KK.NE.1 ) THEN
         CALL XERROR ( 'QS2I1D- the sort control parameter, k, was not 1
     $ OR -1.',55,2,1)
         RETURN
      ENDIF
C
C     Alter array IA to get decreasing order if needed.
C
      IF( KFLAG.LT.1 ) THEN
         DO 20 I=1,NN
            IA(I) = -IA(I)
 20      CONTINUE
      ENDIF
C
C     Sort IA and carry JA and A along.
C     And now...Just a little black magic...
      M = 1
      I = 1
      J = NN
      R = .375
 210  IF( R.LE.0.5898437 ) THEN
         R = R + 3.90625E-2
      ELSE
         R = R-.21875
      ENDIF
 225  K = I
C
C     Select a central element of the array and save it in location 
C     it, jt, at.
C
      IJ = I + IDINT( DBLE(J-I)*R )
      IT = IA(IJ)
      JT = JA(IJ)
      TA = A(IJ)
C
C     If first element of array is greater than it, interchange with it.
C
      IF( IA(I).GT.IT ) THEN
         IA(IJ) = IA(I)
         IA(I)  = IT
         IT     = IA(IJ)
         JA(IJ) = JA(I)
         JA(I)  = JT
         JT     = JA(IJ)
         A(IJ)  = A(I)
         A(I)   = TA
         TA     = A(IJ)
      ENDIF
      L=J
C                           
C     If last element of array is less than it, swap with it.
C
      IF( IA(J).LT.IT ) THEN
         IA(IJ) = IA(J)
         IA(J)  = IT
         IT     = IA(IJ)
         JA(IJ) = JA(J)
         JA(J)  = JT
         JT     = JA(IJ)
         A(IJ)  = A(J)
         A(J)   = TA
         TA     = A(IJ)
C
C     If first element of array is greater than it, swap with it.
C
         IF ( IA(I).GT.IT ) THEN
            IA(IJ) = IA(I)
            IA(I)  = IT
            IT     = IA(IJ)
            JA(IJ) = JA(I)
            JA(I)  = JT
            JT     = JA(IJ)
            A(IJ)  = A(I)
            A(I)   = TA
            TA     = A(IJ)
         ENDIF
      ENDIF
C
C     Find an element in the second half of the array which is 
C     smaller than it.
C
  240 L=L-1
      IF( IA(L).GT.IT ) GO TO 240
C
C     Find an element in the first half of the array which is 
C     greater than it.
C
  245 K=K+1
      IF( IA(K).LT.IT ) GO TO 245
C
C     Interchange these elements.
C
      IF( K.LE.L ) THEN
         IIT   = IA(L)
         IA(L) = IA(K)
         IA(K) = IIT
         JJT   = JA(L)
         JA(L) = JA(K)
         JA(K) = JJT
         TTA   = A(L)
         A(L)  = A(K)
         A(K)  = TTA
         GOTO 240
      ENDIF
C
C     Save upper and lower subscripts of the array yet to be sorted.
C
      IF( L-I.GT.J-K ) THEN
         IL(M) = I
         IU(M) = L
         I = K
         M = M+1
      ELSE
         IL(M) = K
         IU(M) = J
         J = L
         M = M+1
      ENDIF
      GO TO 260
C
C     Begin again on another portion of the unsorted array.
C                                  
  255 M = M-1
      IF( M.EQ.0 ) GO TO 300
      I = IL(M)
      J = IU(M)
  260 IF( J-I.GE.1 ) GO TO 225
      IF( I.EQ.J ) GO TO 255
      IF( I.EQ.1 ) GO TO 210
      I = I-1
  265 I = I+1
      IF( I.EQ.J ) GO TO 255
      IT = IA(I+1)
      JT = JA(I+1)
      TA =  A(I+1)
      IF( IA(I).LE.IT ) GO TO 265
      K=I
  270 IA(K+1) = IA(K)
      JA(K+1) = JA(K)
      A(K+1)  =  A(K)
      K = K-1
      IF( IT.LT.IA(K) ) GO TO 270
      IA(K+1) = IT
      JA(K+1) = JT
      A(K+1)  = TA
      GO TO 265
C
C     Clean up, if necessary.
C
  300 IF( KFLAG.LT.1 ) THEN
         DO 310 I=1,NN
            IA(I) = -IA(I)
 310     CONTINUE
      ENDIF
      RETURN
C------------- LAST LINE OF QS2I1D FOLLOWS ----------------------------
      END
C======================================================================C
C
      subroutine xerror(messg,nmessg,nerr,level)
c***description
c
c     abstract
c        xerror processes a diagnostic message, in a manner
c        determined by the value of level and the current value
c        of the library error control flag, kontrl.
c        (see subroutine xsetf for details.)
c
c     description of parameters
c      --input--
c        messg - the hollerith message to be processed, containing
c                no more than 72 characters.
c        nmessg- the actual number of characters in messg.
c        nerr  - the error number associated with this message.
c                nerr must not be zero.
c        level - error category.
c                =2 means this is an unconditionally fatal error.
c                =1 means this is a recoverable error.  (i.e., it is
c                   non-fatal if xsetf has been appropriately called.)
c                =0 means this is a warning message only.
c                =-1 means this is a warning message which is to be
c                   printed at most once, regardless of how many
c                   times this call is executed.
c
      character*(*) messg
c***first executable statement  xerror
      call xerrwv(messg,nmessg,nerr,level,0,0,0,0,0.,0.)
      return
      end

C======================================================================C
C
      subroutine xerrwv(messg,nmessg,nerr,level,ni,i1,i2,nr,r1,r2)
c***purpose  process an error message allowing 2 integer and 2 real
c            values to be included in the message.
c***description
c
c     abstract
c        xerrwv processes a diagnostic message, in a manner
c        determined by the value of level and the current value
c        of the library error control flag, kontrl.
c        (see subroutine xsetf for details.)
c        in addition, up to two integer values and two real
c        values may be printed along with the message.
c
c     description of parameters
c      --input--
c        messg - the hollerith message to be processed.
c        nmessg- the actual number of characters in messg.
c        nerr  - the error number associated with this message.
c                nerr must not be zero.
c        level - error category.
c                =2 means this is an unconditionally fatal error.
c                =1 means this is a recoverable error.  (i.e., it is
c                   non-fatal if xsetf has been appropriately called.)
c                =0 means this is a warning message only.
c                =-1 means this is a warning message which is to be
c                   printed at most once, regardless of how many
c                   times this call is executed.
c        ni    - number of integer values to be printed. (0 to 2)
c        i1    - first integer value.
c        i2    - second integer value.
c        nr    - number of real values to be printed. (0 to 2)
c        r1    - first real value.
c        r2    - second real value.
c
c     examples
c        call xerrwv('smooth -- num (=i1) was zero.',29,1,2,
c    1   1,num,0,0,0.,0.)
c        call xerrwv('quadxy -- requested error (r1) less than minimum (
c    1r2).,54,77,1,0,0,0,2,errreq,errmin)
c
c     latest revision ---  1 august 1985
c     written by ron jones, with slatec common math library subcommittee
c***references  jones r.e., kahaner d.k., 'xerror, the slatec error-
c                 handling package', sand82-0800, sandia laboratories,
c                 1982.
c***routines called  fdump,i1mach,j4save,xerabt,xerctl,xerprt,xersav,
c                    xgetua
c***end prologue  xerrwv
      character*(*) messg
      character*20 lfirst
      character*37 form
      dimension lun(5)
c     get flags
c***first executable statement  xerrwv
      lkntrl = j4save(2,0,.false.)
      maxmes = j4save(4,0,.false.)
c     check for valid input
      if ((nmessg.gt.0).and.(nerr.ne.0).and.
     1    (level.ge.(-1)).and.(level.le.2)) go to 10
         if (lkntrl.gt.0) call xerprt('fatal error in...',17)
         call xerprt('xerror -- invalid input',23)
c        if (lkntrl.gt.0) call fdump
         if (lkntrl.gt.0) call xerprt('job abort due to fatal error.',
     1  29)
         if (lkntrl.gt.0) call xersav(' ',0,0,0,kdummy)
CC         call xerabt('xerror -- invalid input',23)
         return
   10 continue
c     record message
      junk = j4save(1,nerr,.true.)
      call xersav(messg,nmessg,nerr,level,kount)
c     let user override
      lfirst = messg
      lmessg = nmessg
      lerr = nerr
      llevel = level
      call xerctl(lfirst,lmessg,lerr,llevel,lkntrl)
c     reset to original values
      lmessg = nmessg
      lerr = nerr
      llevel = level
      lkntrl = max0(-2,min0(2,lkntrl))
      mkntrl = iabs(lkntrl)
c     decide whether to print message
      if ((llevel.lt.2).and.(lkntrl.eq.0)) go to 100
      if (((llevel.eq.(-1)).and.(kount.gt.min0(1,maxmes)))
     1.or.((llevel.eq.0)   .and.(kount.gt.maxmes))
     2.or.((llevel.eq.1)   .and.(kount.gt.maxmes).and.(mkntrl.eq.1))
     3.or.((llevel.eq.2)   .and.(kount.gt.max0(1,maxmes)))) go to 100
         if (lkntrl.le.0) go to 20
            call xerprt(' ',1)
c           introduction
            if (llevel.eq.(-1)) call xerprt
     1('warning message...this message will only be printed once.',57)
            if (llevel.eq.0) call xerprt('warning in...',13)
            if (llevel.eq.1) call xerprt
     1      ('recoverable error in...',23)
            if (llevel.eq.2) call xerprt('fatal error in...',17)
   20    continue
c        message
         call xerprt(messg,lmessg)
         call xgetua(lun,nunit)
         isizei = log10(float(i1mach(9))) + 1.0
         isizef = log10(float(i1mach(10))**i1mach(11)) + 1.0
         do 50 kunit=1,nunit
            iunit = lun(kunit)
            if (iunit.eq.0) iunit = i1mach(4)
            do 22 i=1,min(ni,2)
               write (form,21) i,isizei
   21          format ('(11x,21hin above message, i',i1,'=,i',i2,')   ')
               if (i.eq.1) write (iunit,form) i1
               if (i.eq.2) write (iunit,form) i2
   22       continue
            do 24 i=1,min(nr,2)
               write (form,23) i,isizef+10,isizef
   23          format ('(11x,21hin above message, r',i1,'=,e',
     1         i2,'.',i2,')')
               if (i.eq.1) write (iunit,form) r1
               if (i.eq.2) write (iunit,form) r2
   24       continue
            if (lkntrl.le.0) go to 40
c              error number
               write (iunit,30) lerr
   30          format (15h error number =,i10)
   40       continue
   50    continue
c        trace-back
c        if (lkntrl.gt.0) call fdump
  100 continue
      ifatal = 0
      if ((llevel.eq.2).or.((llevel.eq.1).and.(mkntrl.eq.2)))
     1ifatal = 1
c     quit here if message is not fatal
      if (ifatal.le.0) return
      if ((lkntrl.le.0).or.(kount.gt.max0(1,maxmes))) go to 120
c        print reason for abort
         if (llevel.eq.1) call xerprt
     1   ('job abort due to unrecovered error.',35)
         if (llevel.eq.2) call xerprt
     1   ('job abort due to fatal error.',29)
c        print error summary
         call xersav(' ',-1,0,0,kdummy)
  120 continue
c     abort
      if ((llevel.eq.2).and.(kount.gt.max0(1,maxmes))) lmessg = 0
CC      call xerabt(messg,lmessg)
      return
      end

C======================================================================C
C
      subroutine xersav(messg,nmessg,nerr,level,icount)
c***purpose  record that an error has occurred.
c***description
c
c     abstract
c        record that this error occurred.
c
c     description of parameters
c     --input--
c       messg, nmessg, nerr, level are as in xerror,
c       except that when nmessg=0 the tables will be
c       dumped and cleared, and when nmessg is less than zero the
c       tables will be dumped and not cleared.
c     --output--
c       icount will be the number of times this message has
c       been seen, or zero if the table has overflowed and
c       does not contain this message specifically.
c       when nmessg=0, icount will not be altered.
c
c     written by ron jones, with slatec common math library subcommittee
c     latest revision ---  1 august 1985
      integer lun(5)
      character*(*) messg
      character*20 mestab(10),mes
      dimension nertab(10),levtab(10),kount(10)
      save mestab,nertab,levtab,kount,kountx
c     next two data statements are necessary to provide a blank
c     error table initially
      data kount(1),kount(2),kount(3),kount(4),kount(5),
     1     kount(6),kount(7),kount(8),kount(9),kount(10)
     2     /0,0,0,0,0,0,0,0,0,0/
      data kountx/0/
c***first executable statement  xersav
      if (nmessg.gt.0) go to 80
c     dump the table
         if (kount(1).eq.0) return
c        print to each unit
         call xgetua(lun,nunit)
         do 60 kunit=1,nunit
            iunit = lun(kunit)
            if (iunit.eq.0) iunit = i1mach(4)
c           print table header
            write (iunit,10)
   10       format (32h0          error message summary/
     1      51h message start             nerr     level     count)
c           print body of table
            do 20 i=1,10
               if (kount(i).eq.0) go to 30
               write (iunit,15) mestab(i),nertab(i),levtab(i),kount(i)
   15          format (1x,a20,3i10)
   20       continue
   30       continue
c           print number of other errors
            if (kountx.ne.0) write (iunit,40) kountx
   40       format (41h0other errors not individually tabulated=,i10)
            write (iunit,50)
   50       format (1x)
   60    continue
         if (nmessg.lt.0) return
c        clear the error tables
         do 70 i=1,10
   70       kount(i) = 0
         kountx = 0
         return
   80 continue
c     process a message...
c     search for this messg, or else an empty slot for this messg,
c     or else determine that the error table is full.
      mes = messg
      do 90 i=1,10
         ii = i
         if (kount(i).eq.0) go to 110
         if (mes.ne.mestab(i)) go to 90
         if (nerr.ne.nertab(i)) go to 90
         if (level.ne.levtab(i)) go to 90
         go to 100
   90 continue
c     three possible cases...
c     table is full
         kountx = kountx+1
         icount = 1
         return
c     message found in table
  100    kount(ii) = kount(ii) + 1
         icount = kount(ii)
         return
c     empty slot found for new message
  110    mestab(ii) = mes
         nertab(ii) = nerr
         levtab(ii) = level
         kount(ii)  = 1
         icount = 1
         return
      end
      subroutine xgetf(kontrl)
c***purpose  return the current value of the error control flag.
c***description
c
c   abstract
c        xgetf returns the current value of the error control flag
c        in kontrl.  see subroutine xsetf for flag value meanings.
c        (kontrl is an output parameter only.)
c
c     written by ron jones, with slatec common math library subcommittee
c     latest revision ---  7 june 1978
c***first executable statement  xgetf
      kontrl = j4save(2,0,.false.)
      return
      end
C
      subroutine xgetua(iunita,n)
c***purpose  return unit number(s) to which error messages are being
c            sent.
c***description
c
c     abstract
c        xgetua may be called to determine the unit number or numbers
c        to which error messages are being sent.
c        these unit numbers may have been set by a call to xsetun,
c        or a call to xsetua, or may be a default value.
c
c     description of parameters
c      --output--
c        iunit - an array of one to five unit numbers, depending
c                on the value of n.  a value of zero refers to the
c                default unit, as defined by the i1mach machine
c                constant routine.  only iunit(1),...,iunit(n) are
c                defined by xgetua.  the values of iunit(n+1),...,
c                iunit(5) are not defined (for n .lt. 5) or altered
c                in any way by xgetua.
c        n     - the number of units to which copies of the
c                error messages are being sent.  n will be in the
c                range from 1 to 5.
c
c     latest revision ---  19 mar 1980
c     written by ron jones, with slatec common math library subcommittee
      dimension iunita(5)
c***first executable statement  xgetua
      n = j4save(5,0,.false.)
      do 30 i=1,n
         index = i+4
         if (i.eq.1) index = 3
         iunita(i) = j4save(index,0,.false.)
   30 continue
      return
      end
C
      function j4save(iwhich,ivalue,iset)
c***begin prologue  j4save
c***refer to  xerror
c***routines called  (none)
c***description
c
c     abstract
c        j4save saves and recalls several global variables needed
c        by the library error handling routines.
c
c     description of parameters
c      --input--
c        iwhich - index of item desired.
c                = 1 refers to current error number.
c                = 2 refers to current error control flag.
c                 = 3 refers to current unit number to which error
c                    messages are to be sent.  (0 means use standard.)
c                 = 4 refers to the maximum number of times any
c                     message is to be printed (as set by xermax).
c                 = 5 refers to the total number of units to which
c                     each error message is to be written.
c                 = 6 refers to the 2nd unit for error messages
c                 = 7 refers to the 3rd unit for error messages
c                 = 8 refers to the 4th unit for error messages
c                 = 9 refers to the 5th unit for error messages
c        ivalue - the value to be set for the iwhich-th parameter,
c                 if iset is .true. .
c        iset   - if iset=.true., the iwhich-th parameter will be
c                 given the value, ivalue.  if iset=.false., the
c                 iwhich-th parameter will be unchanged, and ivalue
c                 is a dummy parameter.
c      --output--
c        the (old) value of the iwhich-th parameter will be returned
c        in the function value, j4save.
c
c     written by ron jones, with slatec common math library subcommittee
c    adapted from bell laboratories port library error handler
c     latest revision ---  1 august 1985
c***references  jones r.e., kahaner d.k., 'xerror, the slatec error-
c                 handling package', sand82-0800, sandia laboratories,
c                 1982.
c***end prologue  j4save
      logical iset
      integer iparam(9)
      save iparam
      data iparam(1),iparam(2),iparam(3),iparam(4)/0,2,0,10/
      data iparam(5)/1/
      data iparam(6),iparam(7),iparam(8),iparam(9)/0,0,0,0/
c***first executable statement  j4save
      j4save = iparam(iwhich)
      if (iset) iparam(iwhich) = ivalue
      return
      end
C
      subroutine xerclr
c***purpose  reset current error number to zero.
c***description
c
c     abstract
c        this routine simply resets the current error number to zero.
c        this may be necessary to do in order to determine that
c        a certain error has occurred again since the last time
c        numxer was referenced.
c
c     written by ron jones, with slatec common math library subcommittee
c***references  jones r.e., kahaner d.k., 'xerror, the slatec error-
c                 handling package', sand82-0800, sandia laboratories,
c                 1982.
c***routines called  j4save
c***end prologue  xerclr
c***first executable statement  xerclr
      junk = j4save(1,0,.true.)
      return
      end
C
      subroutine xerdmp
c***purpose  print the error tables and then clear them.
c***description
c
c     abstract
c        xerdmp prints the error tables, then clears them.
c
c     written by ron jones, with slatec common math library subcommittee
c     latest revision ---  7 june 1978
c***references  jones r.e., kahaner d.k., 'xerror, the slatec error-
c                 handling package', sand82-0800, sandia laboratories,
c                 1982.
c***routines called  xersav
c***end prologue  xerdmp
c***first executable statement  xerdmp
      call xersav(' ',0,0,0,kount)
      return
      end
C
      subroutine xermax(max)
c***purpose  set maximum number of times any error message is to be
c            printed.
c***description
c
c     abstract
c        xermax sets the maximum number of times any message
c        is to be printed.  that is, non-fatal messages are
c        not to be printed after they have occured max times.
c        such non-fatal messages may be printed less than
c        max times even if they occur max times, if error
c        suppression mode (kontrl=0) is ever in effect.
c
c     description of parameter
c      --input--
c        max - the maximum number of times any one message
c              is to be printed.
c
c     written by ron jones, with slatec common math library subcommittee
c     latest revision ---  7 june 1978
c***references  jones r.e., kahaner d.k., 'xerror, the slatec error-
c                 handling package', sand82-0800, sandia laboratories,
c                 1982.
c***routines called  j4save
c***end prologue  xermax
c***first executable statement  xermax
      junk = j4save(4,max,.true.)
      return
      end
C
      subroutine xgetun(iunit)
c***purpose  return the (first) output file to which error messages
c            are being sent.
c***description
c
c     abstract
c        xgetun gets the (first) output file to which error messages
c        are being sent.  to find out if more than one file is being
c        used, one must use the xgetua routine.
c
c     description of parameter
c      --output--
c        iunit - the logical unit number of the  (first) unit to
c                which error messages are being sent.
c                a value of zero means that the default file, as
c                defined by the i1mach routine, is being used.
c
c     written by ron jones, with slatec common math library subcommittee
c     latest revision --- 23 may 1979
c***references  jones r.e., kahaner d.k., 'xerror, the slatec error-
c                 handling package', sand82-0800, sandia laboratories,
c                 1982.
c***routines called  j4save
c***end prologue  xgetun
c***first executable statement  xgetun
      iunit = j4save(3,0,.false.)
      return
      end
C
      subroutine xsetf(kontrl)
c***purpose  set the error control flag.
c***description
c
c     abstract
c        xsetf sets the error control flag value to kontrl.
c        (kontrl is an input parameter only.)
c        the following table shows how each message is treated,
c        depending on the values of kontrl and level.  (see xerror
c        for description of level.)
c
c        if kontrl is zero or negative, no information other than the
c        message itself (including numeric values, if any) will be
c        printed.  if kontrl is positive, introductory messages,
c        trace-backs, etc., will be printed in addition to the message.
c
c              iabs(kontrl)
c        level        0              1              2
c        value
c          2        fatal          fatal          fatal
c
c          1     not printed      printed         fatal
c
c          0     not printed      printed        printed
c
c         -1     not printed      printed        printed
c                                  only           only
c                                  once           once
c
c     written by ron jones, with slatec common math library subcommittee
c     latest revision ---  19 mar 1980
c***first executable statement  xsetf
      if ((kontrl.ge.(-2)).and.(kontrl.le.2)) go to 10
         call xerrwv('xsetf  -- invalid value of kontrl (i1).',33,1,2,
     1  1,kontrl,0,0,0.,0.)
         return
   10 junk = j4save(2,kontrl,.true.)
      return
      end
C
      subroutine xsetua(iunita,n)
c***purpose  set logical unit numbers (up to 5) to which error
c            messages are to be sent.
c***description
c
c     abstract
c        xsetua may be called to declare a list of up to five
c        logical units, each of which is to receive a copy of
c        each error message processed by this package.
c        the purpose of xsetua is to allow simultaneous printing
c        of each error message on, say, a main output file,
c        an interactive terminal, and other files such as graphics
c        communication files.
c
c     description of parameters
c      --input--
c        iunit - an array of up to five unit numbers.
c                normally these numbers should all be different
c                (but duplicates are not prohibited.)
c        n     - the number of unit numbers provided in iunit
c                must have 1 .le. n .le. 5.
c
c     written by ron jones, with slatec common math library subcommittee
c     latest revision ---  19 mar 1980
      dimension iunita(5)
c***first executable statement  xsetua
      if ((n.ge.1).and.(n.le.5)) go to 10
         call xerrwv('xsetua -- invalid value of n (i1).',34,1,2,
     1  1,n,0,0,0.,0.)
         return
   10 continue
      do 20 i=1,n
         index = i+4
         if (i.eq.1) index = 3
         junk = j4save(index,iunita(i),.true.)
   20 continue
      junk = j4save(5,n,.true.)
      return
      end
C
      subroutine xsetun(iunit)
c***purpose  set output file to which error messages are to be sent.
c***description
c
c     abstract
c        xsetun sets the output file to which error messages are to
c        be sent.  only one file will be used.  see xsetua for
c        how to declare more than one file.
c
c     description of parameter
c      --input--
c        iunit - an input parameter giving the logical unit number
c                to which error messages are to be sent.
c
c     written by ron jones, with slatec common math library subcommittee
c     latest revision ---  7 june 1978
c***first executable statement  xsetun
      junk = j4save(3,iunit,.true.)
      junk = j4save(5,1,.true.)
      return
      end
C
C======================================================================C
C
      INTEGER FUNCTION I1MACH(I)
C
C  I/O UNIT NUMBERS.
C
C    I1MACH( 1) = THE STANDARD INPUT UNIT.
C
C    I1MACH( 2) = THE STANDARD OUTPUT UNIT.
C
C    I1MACH( 3) = THE STANDARD PUNCH UNIT.
C
C    I1MACH( 4) = THE STANDARD ERROR MESSAGE UNIT.
C
C  WORDS.
C
C    I1MACH( 5) = THE NUMBER OF BITS PER INTEGER STORAGE UNIT.
C
C    I1MACH( 6) = THE NUMBER OF CHARACTERS PER INTEGER STORAGE UNIT.
C
C  INTEGERS.
C
C    ASSUME INTEGERS ARE REPRESENTED IN THE S-DIGIT, BASE-A FORM
C
C               SIGN ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) )
C
C               WHERE 0 .LE. X(I) .LT. A FOR I=0,...,S-1.
C
C    I1MACH( 7) = A, THE BASE.
C
C    I1MACH( 8) = S, THE NUMBER OF BASE-A DIGITS.
C
C    I1MACH( 9) = A**S - 1, THE LARGEST MAGNITUDE.
C
C  FLOATING-POINT NUMBERS.
C
C    ASSUME FLOATING-POINT NUMBERS ARE REPRESENTED IN THE T-DIGIT,
C    BASE-B FORM
C
C               SIGN (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
C
C               WHERE 0 .LE. X(I) .LT. B FOR I=1,...,T,
C               0 .LT. X(1), AND EMIN .LE. E .LE. EMAX.
C
C    I1MACH(10) = B, THE BASE.
C
C  SINGLE-PRECISION
C
C    I1MACH(11) = T, THE NUMBER OF BASE-B DIGITS.
C
C    I1MACH(12) = EMIN, THE SMALLEST EXPONENT E.
C
C    I1MACH(13) = EMAX, THE LARGEST EXPONENT E.
C
C  DOUBLE-PRECISION
C
C    I1MACH(14) = T, THE NUMBER OF BASE-B DIGITS.
C
C    I1MACH(15) = EMIN, THE SMALLEST EXPONENT E.
C
C    I1MACH(16) = EMAX, THE LARGEST EXPONENT E.
C
C  TO ALTER THIS FUNCTION FOR A PARTICULAR ENVIRONMENT,
C  THE DESIRED SET OF DATA STATEMENTS SHOULD BE ACTIVATED BY
C  REMOVING THE C FROM COLUMN 1.  ALSO, THE VALUES OF
C  I1MACH(1) - I1MACH(4) SHOULD BE CHECKED FOR CONSISTENCY
C  WITH THE LOCAL OPERATING SYSTEM.
C  ON RARE MACHINES A STATIC STATEMENT MAY NEED TO BE ADDED.
C  (BUT PROBABLY MORE SYSTEMS PROHIBIT IT THAN REQUIRE IT.)
C
      INTEGER IMACH(16),OUTPUT
C
      EQUIVALENCE (IMACH(4),OUTPUT)
C
C     MACHINE CONSTANTS FOR IEEE ARITHMETIC MACHINES (E.G., AT&T 3B
C     SERIES COMPUTERS AND 8087-BASED MACHINES LIKE THE IBM PC).
C
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /    7 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   32 /
C      DATA IMACH( 6) /    4 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   31 /
C      DATA IMACH( 9) / 2147483647 /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   24 /
C      DATA IMACH(12) / -125 /
C      DATA IMACH(13) /  128 /
C      DATA IMACH(14) /   53 /
C      DATA IMACH(15) / -1021 /
C      DATA IMACH(16) /  1024 /
C
C     MACHINE CONSTANTS FOR THE INTERDATA 8/32
C     WITH THE UNIX SYSTEM FORTRAN 77 COMPILER.
C
C     FOR THE INTERDATA FORTRAN VII COMPILER REPLACE
C     THE Z'S SPECIFYING HEX CONSTANTS WITH Y'S.
C
C      DATA IMACH( 1) /   5 /
C      DATA IMACH( 2) /   6 /
C      DATA IMACH( 3) /   6 /
C      DATA IMACH( 4) /   6 /
C      DATA IMACH( 5) /  32 /
C      DATA IMACH( 6) /   4 /
C      DATA IMACH( 7) /   2 /
C      DATA IMACH( 8) /  31 /
C      DATA IMACH( 9) / Z'7FFFFFFF' /
C      DATA IMACH(10) /  16 /
C      DATA IMACH(11) /   6 /
C      DATA IMACH(12) / -64 /
C      DATA IMACH(13) /  62 /
C      DATA IMACH(14) /  14 /
C      DATA IMACH(15) / -64 /
C      DATA IMACH(16) /  62 /
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRANS SUPPORTING
C     32-BIT INTEGER ARITHMETIC.
C
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /    7 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   32 /
C      DATA IMACH( 6) /    4 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   31 /
C      DATA IMACH( 9) / 2147483647 /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   24 /
C      DATA IMACH(12) / -127 /
C      DATA IMACH(13) /  127 /
C      DATA IMACH(14) /   56 /
C      DATA IMACH(15) / -127 /
C      DATA IMACH(16) /  127 /
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRANS SUPPORTING
C     16-BIT INTEGER ARITHMETIC.
C
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /    7 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   16 /
C      DATA IMACH( 6) /    2 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   15 /
C      DATA IMACH( 9) / 32767 /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   24 /
C      DATA IMACH(12) / -127 /
C      DATA IMACH(13) /  127 /
C      DATA IMACH(14) /   56 /
C      DATA IMACH(15) / -127 /
C      DATA IMACH(16) /  127 /
C
C     MACHINE CONSTANTS FOR THE SUN MICROSYSTEMS UNIX F77 COMPILER.
C
      DATA IMACH( 1) /     5 /
      DATA IMACH( 2) /     6 /
      DATA IMACH( 3) /     6 /
      DATA IMACH( 4) /     0 /
      DATA IMACH( 5) /    32 /
      DATA IMACH( 6) /     4 /
      DATA IMACH( 7) /     2 /
      DATA IMACH( 8) /    32 /
      DATA IMACH( 9) /2147483647/
      DATA IMACH(10) /     2 /
      DATA IMACH(11) /    24 /
      DATA IMACH(12) /  -126 /
      DATA IMACH(13) /   128 /
      DATA IMACH(14) /    53 /
      DATA IMACH(15) / -1022 /
      DATA IMACH(16) /  1024 /
C
      IF (I .LT. 1  .OR.  I .GT. 16) GO TO 999
      I1MACH=IMACH(I)
      RETURN
  999 WRITE(OUTPUT,1999) I
 1999 FORMAT(' I1MACH - I OUT OF BOUNDS',I10)
      STOP
      END

C======================================================================C
C
      REAL*8 FUNCTION D1MACH(I)
C
C  DOUBLE-PRECISION MACHINE CONSTANTS
C
C  D1MACH( 1) = B**(EMIN-1), THE SMALLEST POSITIVE MAGNITUDE.
C
C  D1MACH( 2) = B**EMAX*(1 - B**(-T)), THE LARGEST MAGNITUDE.
C
C  D1MACH( 3) = B**(-T), THE SMALLEST RELATIVE SPACING.
C
C  D1MACH( 4) = B**(1-T), THE LARGEST RELATIVE SPACING.
C
C  D1MACH( 5) = LOG10(B)
C
C  TO ALTER THIS FUNCTION FOR A PARTICULAR ENVIRONMENT,
C  THE DESIRED SET OF DATA STATEMENTS SHOULD BE ACTIVATED BY
C  REMOVING THE C FROM COLUMN 1.
C  ON RARE MACHINES A STATIC STATEMENT MAY NEED TO BE ADDED.
C  (BUT PROBABLY MORE SYSTEMS PROHIBIT IT THAN REQUIRE IT.)
C
C  WHERE POSSIBLE, OCTAL OR HEXADECIMAL CONSTANTS HAVE BEEN USED
C  TO SPECIFY THE CONSTANTS EXACTLY WHICH HAS IN SOME CASES
C  REQUIRED THE USE OF EQUIVALENT INTEGER ARRAYS.
C
      INTEGER SMALL(4)
      INTEGER LARGE(4)
      INTEGER RIGHT(4)
      INTEGER DIVER(4)
      INTEGER LOG10(4)
C
      REAL*8 DMACH(5)
C
      EQUIVALENCE (DMACH(1),SMALL(1))
      EQUIVALENCE (DMACH(2),LARGE(1))
      EQUIVALENCE (DMACH(3),RIGHT(1))
      EQUIVALENCE (DMACH(4),DIVER(1))
      EQUIVALENCE (DMACH(5),LOG10(1))
C
C     MACHINE CONSTANTS FOR THE IBM PC AND OTHER 8087-ARITHMETIC MICROS
C
C      DATA SMALL(1),SMALL(2) /          0,    1048576 /
C      DATA LARGE(1),LARGE(2) /         -1, 2146435071 /
C      DATA RIGHT(1),RIGHT(2) /          0, 1017118720 /
C      DATA DIVER(1),DIVER(2) /          0, 1018167296 /
C      DATA LOG10(1),LOG10(2) / 1352628735, 1070810131 /
C
C     MACHINE CONSTANTS FOR THE INTERDATA 8/32
C     WITH THE UNIX SYSTEM FORTRAN 77 COMPILER.
C
C     FOR THE INTERDATA FORTRAN VII COMPILER REPLACE
C     THE Z'S SPECIFYING HEX CONSTANTS WITH Y'S.
C
C      DATA SMALL(1),SMALL(2) / Z'00100000', Z'00000000' /
C      DATA LARGE(1),LARGE(2) / Z'7EFFFFFF', Z'FFFFFFFF' /
C      DATA RIGHT(1),RIGHT(2) / Z'33100000', Z'00000000' /
C      DATA DIVER(1),DIVER(2) / Z'34100000', Z'00000000' /
C      DATA LOG10(1),LOG10(2) / Z'41134413', Z'509F79FF' /
C
C     MACHINE CONSTANTS FOR THE SUN MICROSYSTEMS UNIX F77 COMPILER.
C
      DATA DMACH(1) / 2.22507385850720D-38 /
      DATA DMACH(2) / 1.79769313486231D+38 /
      DATA DMACH(3) / 1.1101827117665D-16 /
      DATA DMACH(4) / 2.2203654423533D-16 /
      DATA DMACH(5) / 3.01029995663981E-1 /
C
      IF (I .LT. 1  .OR.  I .GT. 5) GOTO 999
      D1MACH = DMACH(I)
      RETURN
  999 WRITE(I1MACH(2),1999) I
 1999 FORMAT(' D1MACH - I OUT OF BOUNDS',I10)
      STOP
      END
C======================================================================C
C
      SUBROUTINE DAXPY(N,DA,DX,INCX,DY,INCY)
C
C     OVERWRITE DOUBLE PRECISION DY WITH DOUBLE PRECISION DA*DX + DY.
C     FOR I = 0 TO N-1, REPLACE  DY(LY+I*INCY) WITH DA*DX(LX+I*INCX) +
C       DY(LY+I*INCY), WHERE LX = 1 IF INCX .GE. 0, ELSE LX = (-INCX)*N,
C       AND LY IS DEFINED IN A SIMILAR WAY USING INCY.
C
      REAL*8 DX(1),DY(1),DA
      IF(N.LE.0.OR.DA.EQ.0.D0) RETURN
      IF(INCX.EQ.INCY) IF(INCX-1) 5,20,60
    5 CONTINUE
C
C        CODE FOR NONEQUAL OR NONPOSITIVE INCREMENTS.
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DY(IY) = DY(IY) + DA*DX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 4.
C
   20 M = MOD(N,4)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DY(I) = DY(I) + DA*DX(I)
   30 CONTINUE
      IF( N .LT. 4 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,4
        DY(I) = DY(I) + DA*DX(I)
        DY(I + 1) = DY(I + 1) + DA*DX(I + 1)
        DY(I + 2) = DY(I + 2) + DA*DX(I + 2)
        DY(I + 3) = DY(I + 3) + DA*DX(I + 3)
   50 CONTINUE
      RETURN
C
C        CODE FOR EQUAL, POSITIVE, NONUNIT INCREMENTS.
C
   60 CONTINUE
      NS = N*INCX
          DO 70 I=1,NS,INCX
          DY(I) = DA*DX(I) + DY(I)
   70     CONTINUE
      RETURN
      END

C======================================================================C
C
      REAL*8 FUNCTION DDOT(N,DX,INCX,DY,INCY)
C
C     RETURNS THE DOT PRODUCT OF DOUBLE PRECISION DX AND DY.
C     DDOT = SUM FOR I = 0 TO N-1 OF  DX(LX+I*INCX) * DY(LY+I*INCY)
C     WHERE LX = 1 IF INCX .GE. 0, ELSE LX = (-INCX)*N, AND LY IS
C     DEFINED IN A SIMILAR WAY USING INCY.
C
      REAL*8 DX(1),DY(1)
      DDOT = 0.D0
      IF(N.LE.0)RETURN
      IF(INCX.EQ.INCY) IF(INCX-1) 5,20,60
    5 CONTINUE
C
C         CODE FOR UNEQUAL OR NONPOSITIVE INCREMENTS.
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
         DDOT = DDOT + DX(IX)*DY(IY)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1.
C
C
C        CLEAN-UP LOOP SO REMAINING VECTOR LENGTH IS A MULTIPLE OF 5.
C
   20 M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
         DDOT = DDOT + DX(I)*DY(I)
   30 CONTINUE
      IF( N .LT. 5 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
         DDOT = DDOT + DX(I)*DY(I) + DX(I+1)*DY(I+1) +
     1   DX(I + 2)*DY(I + 2) + DX(I + 3)*DY(I + 3) + DX(I + 4)*DY(I + 4)
   50 CONTINUE
      RETURN
C
C         CODE FOR POSITIVE EQUAL INCREMENTS .NE.1.
C
   60 CONTINUE
      NS = N*INCX
          DO 70 I=1,NS,INCX
          DDOT = DDOT + DX(I)*DY(I)
   70     CONTINUE
      RETURN
      END

C======================================================================C
C======================================================================C
C
      SUBROUTINE DSMV1( N, X, Y, NELT, IA, JA, A, ISYM )
C***PURPOSE  SLAP Column Format Sparse Matrix Vector Product.
C            Routine to calculate the sparse matrix vector product:
C            Y = A*X.
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER N, NELT, IA(NELT), JA(NELT), ISYM
      REAL*8 A(NELT), X(N), Y(N)
C
C         Zero out the result vector.
C***FIRST EXECUTABLE STATEMENT  DSMV
      DO 10 I = 1, N
         Y(I) = 0.0D0
 10   CONTINUE
C
C         Multiply by A.
C
      DO 30 ICOL = 1, N
         IBGN = JA(ICOL)
         IEND = JA(ICOL+1)-1
         DO 20 I = IBGN, IEND
            Y(IA(I)) = Y(IA(I)) + A(I)*X(ICOL)
 20      CONTINUE
 30   CONTINUE
C
      IF( ISYM.EQ.1 ) THEN
C
C         The matrix is non-symmetric.  Need to get the other half in...
C         This loops assumes that the diagonal is the first entry in
C         each column.
C
         DO 50 IROW = 1, N
            JBGN = JA(IROW)+1
            JEND = JA(IROW+1)-1
            IF( JBGN.GT.JEND ) GOTO 50
            DO 40 J = JBGN, JEND
               Y(IROW) = Y(IROW) + A(J)*X(IA(J))
 40         CONTINUE
 50      CONTINUE
      ENDIF
      RETURN
C------------- LAST LINE OF DSMV FOLLOWS ----------------------------
      END

C======================================================================C
      SUBROUTINE DSLUI(N, B, X, NELT, IA, JA, A, ISYM, RWORK, IWORK)
C***PURPOSE  SLAP for LDU Factorization.
C            This routine  acts as an  interface between  the   SLAP
C            generic MSLOVE calling convention and the routine  that 
C            actually computes:     -1
C                              (LDU)  B = X.
      IMPLICIT REAL*8(A-H,O-Z)
!      INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, IWORK(10)
!      REAL*8 B(N), X(N), A(NELT), RWORK(1)
      INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, IWORK(*)
      REAL*8 B(N), X(N), A(NELT), RWORK(*)
C
C         Pull out the locations of the arrays holding the ILU
C         factorization.
C***FIRST EXECUTABLE STATEMENT  DSLUI
      LOCIL  = IWORK(1)
      LOCJL  = IWORK(2)
      LOCIU  = IWORK(3)
      LOCJU  = IWORK(4)
      LOCL   = IWORK(5)
      LOCDIN = IWORK(6)
      LOCU   = IWORK(7)
C
C         Solve the system LUx = b
      CALL DSLUI2(N, B, X, IWORK(LOCIL), IWORK(LOCJL), RWORK(LOCL),
     $     RWORK(LOCDIN), IWORK(LOCIU), IWORK(LOCJU), RWORK(LOCU) )
C         
      RETURN
C------------- LAST LINE OF DSLUI FOLLOWS ----------------------------
      END

C======================================================================C

      SUBROUTINE DSLUI2(N, B, X, IL, JL, L, DINV, IU, JU, U )
C***PURPOSE  SLAP Back solve for LDU Factorization.
C            Routine  to  solve a system of the form  L*D*U X  =  B,
C            where L is a unit  lower  triangular  matrix,  D  is  a 
C            diagonal matrix, and U is a unit upper triangular matrix.
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER N, IL(1), JL(1), IU(1), JU(1)
      REAL*8 B(N), X(N), L(1), DINV(N), U(1)
C         
C         Solve  L*Y = B,  storing result in X, L stored by rows.
C***FIRST EXECUTABLE STATEMENT  DSLUI2
      DO 10 I = 1, N
         X(I) = B(I)
 10   CONTINUE                 
      DO 30 IROW = 2, N
         JBGN = IL(IROW)
         JEND = IL(IROW+1)-1
         IF( JBGN.LE.JEND ) THEN
            DO 20 J = JBGN, JEND
               X(IROW) = X(IROW) - L(J)*X(JL(J))
 20         CONTINUE
         ENDIF
 30   CONTINUE
C         
C         Solve  D*Z = Y,  storing result in X.
      DO 40 I=1,N
         X(I) = X(I)*DINV(I)
 40   CONTINUE
C         
C         Solve  U*X = Z, U stored by columns.
      DO 60 ICOL = N, 2, -1
         JBGN = JU(ICOL)
         JEND = JU(ICOL+1)-1
         IF( JBGN.LE.JEND ) THEN
            DO 50 J = JBGN, JEND
               X(IU(J)) = X(IU(J)) - U(J)*X(ICOL)
 50         CONTINUE
         ENDIF
 60   CONTINUE
C         
      RETURN
C------------- LAST LINE OF DSLUI2 FOLLOWS ----------------------------
      END

C======================================================================C
      subroutine xerctl(messg1,nmessg,nerr,level,kontrl)
c***purpose  allow user control over handling of errors.
c***description
c
c     abstract
c        allows user control over handling of individual errors.
c        just after each message is recorded, but before it is
c        processed any further (i.e., before it is printed or
c        a decision to abort is made), a call is made to xerctl.
c        if the user has provided his own version of xerctl, he
c        can then override the value of kontrol used in processing
c        this message by redefining its value.
c        kontrl may be set to any value from -2 to 2.
c        the meanings for kontrl are the same as in xsetf, except
c        that the value of kontrl changes only for this message.
c        if kontrl is set to a value outside the range from -2 to 2,
c        it will be moved back into that range.
c
c     description of parameters
c
c      --input--
c        messg1 - the first word (only) of the error message.
c        nmessg - same as in the call to xerror or xerrwv.
c        nerr   - same as in the call to xerror or xerrwv.
c        level  - same as in the call to xerror or xerrwv.
c        kontrl - the current value of the control flag as set
c                 by a call to xsetf.
c
c      --output--
c        kontrl - the new value of kontrl.  if kontrl is not
c                 defined, it will remain at its original value.
c                 this changed value of control affects only
c                 the current occurrence of the current message.
c***references  jones r.e., kahaner d.k., 'xerror, the slatec error-
c                 handling package', sand82-0800, sandia laboratories,
c                 1982.
c***routines called  (none)
c***end prologue  xerctl
      character*20 messg1
c***first executable statement  xerctl
      return
      end
C======================================================================C
C
      REAL*8 FUNCTION DNRM21 ( N, DX, INCX)
      INTEGER          NEXT
      REAL*8   DX(1), CUTLO, CUTHI, HITEST, SUM, XMAX,ZERO,ONE
      DATA   ZERO, ONE /0.0D0, 1.0D0/
C
C     EUCLIDEAN NORM OF THE N-VECTOR STORED IN DX() WITH STORAGE
C     INCREMENT INCX .
C     IF    N .LE. 0 RETURN WITH RESULT = 0.
C     IF N .GE. 1 THEN INCX MUST BE .GE. 1
C
C           C.L.LAWSON, 1978 JAN 08
C
C     FOUR PHASE METHOD     USING TWO BUILT-IN CONSTANTS THAT ARE
C     HOPEFULLY APPLICABLE TO ALL MACHINES.
C         CUTLO = MAXIMUM OF  DSQRT(U/EPS)  OVER ALL KNOWN MACHINES.
C         CUTHI = MINIMUM OF  DSQRT(V)      OVER ALL KNOWN MACHINES.
C     WHERE
C         EPS = SMALLEST NO. SUCH THAT EPS + 1. .GT. 1.
C         U   = SMALLEST POSITIVE NO.   (UNDERFLOW LIMIT)
C         V   = LARGEST  NO.            (OVERFLOW  LIMIT)
C
C     BRIEF OUTLINE OF ALGORITHM..
C
C     PHASE 1    SCANS ZERO COMPONENTS.
C     MOVE TO PHASE 2 WHEN A COMPONENT IS NONZERO AND .LE. CUTLO
C     MOVE TO PHASE 3 WHEN A COMPONENT IS .GT. CUTLO
C     MOVE TO PHASE 4 WHEN A COMPONENT IS .GE. CUTHI/M
C     WHERE M = N FOR X() REAL AND M = 2*N FOR COMPLEX.
C
C     VALUES FOR CUTLO AND CUTHI..
C     FROM THE ENVIRONMENTAL PARAMETERS LISTED IN THE IMSL CONVERTER
C     DOCUMENT THE LIMITING VALUES ARE AS FOLLOWS..
C     CUTLO, S.P.   U/EPS = 2**(-102) FOR  HONEYWELL.  CLOSE SECONDS ARE
C                   UNIVAC AND DEC AT 2**(-103)
C                   THUS CUTLO = 2**(-51) = 4.44089E-16
C     CUTHI, S.P.   V = 2**127 FOR UNIVAC, HONEYWELL, AND DEC.
C                   THUS CUTHI = 2**(63.5) = 1.30438E19
C     CUTLO, D.P.   U/EPS = 2**(-67) FOR HONEYWELL AND DEC.
C                   THUS CUTLO = 2**(-33.5) = 8.23181D-11
C     CUTHI, D.P.   SAME AS S.P.  CUTHI = 1.30438D19
C     DATA CUTLO, CUTHI / 8.232D-11,  1.304D19 /
C     DATA CUTLO, CUTHI / 4.441E-16,  1.304E19 /
      DATA CUTLO, CUTHI / 8.232D-11,  1.304D19 /
C
      IF(N .GT. 0) GO TO 10
         DNRM21  = ZERO
         GO TO 300
C
   10 ASSIGN 30 TO NEXT
      SUM = ZERO
      NN = N * INCX
C                                                 BEGIN MAIN LOOP
      I = 1
   20    GO TO NEXT,(30, 50, 70, 110)
   30 IF( DABS(DX(I)) .GT. CUTLO) GO TO 85
      ASSIGN 50 TO NEXT
      XMAX = ZERO
C
C                        PHASE 1.  SUM IS ZERO
C
   50 IF( DX(I) .EQ. ZERO) GO TO 200
      IF( DABS(DX(I)) .GT. CUTLO) GO TO 85
C
C                                PREPARE FOR PHASE 2.
      ASSIGN 70 TO NEXT
      GO TO 105
C
C                                PREPARE FOR PHASE 4.
C
  100 I = J
      ASSIGN 110 TO NEXT
      SUM = (SUM / DX(I)) / DX(I)
  105 XMAX = DABS(DX(I))
      GO TO 115
C
C                   PHASE 2.  SUM IS SMALL.
C                             SCALE TO AVOID DESTRUCTIVE UNDERFLOW.
C
   70 IF( DABS(DX(I)) .GT. CUTLO ) GO TO 75
C
C                     COMMON CODE FOR PHASES 2 AND 4.
C                     IN PHASE 4 SUM IS LARGE.  SCALE TO AVOID OVERFLOW.
C
  110 IF( DABS(DX(I)) .LE. XMAX ) GO TO 115
         SUM = ONE + SUM * (XMAX / DX(I))**2
         XMAX = DABS(DX(I))
         GO TO 200
C
  115 SUM = SUM + (DX(I)/XMAX)**2
      GO TO 200
C
C
C                  PREPARE FOR PHASE 3.
C
   75 SUM = (SUM * XMAX) * XMAX
C
C
C     FOR REAL OR D.P. SET HITEST = CUTHI/N
C     FOR COMPLEX      SET HITEST = CUTHI/(2*N)
C
   85 HITEST = CUTHI/FLOAT( N )
C
C                   PHASE 3.  SUM IS MID-RANGE.  NO SCALING.
C
      DO 95 J =I,NN,INCX
      IF(DABS(DX(J)) .GE. HITEST) GO TO 100
   95    SUM = SUM + DX(J)**2
      DNRM21 = DSQRT( SUM )
      GO TO 300
C
  200 CONTINUE
      I = I + INCX
      IF ( I .LE. NN ) GO TO 20
C
C              END OF MAIN LOOP.
C
C              COMPUTE SQUARE ROOT AND ADJUST FOR SCALING.
C
      DNRM21 = XMAX * DSQRT(SUM)
  300 CONTINUE
      RETURN
      END

C======================================================================C
      subroutine xerprt(messg,nmessg)
c***purpose  print error messages.
c***description
c
c     abstract
c        print the hollerith message in messg, of length nmessg,
c        on each file indicated by xgetua.
c     latest revision ---  1 august 1985
c***references  jones r.e., kahaner d.k., 'xerror, the slatec error-
c                 handling package', sand82-0800, sandia laboratories,
c                 1982.
c***routines called  i1mach,xgetua
c***end prologue  xerprt
      integer lun(5)
      character*(*) messg
c     obtain unit numbers and write line to each unit
c***first executable statement  xerprt
      call xgetua(lun,nunit)
      lenmes = len(messg)
      do 20 kunit=1,nunit
         iunit = lun(kunit)
         if (iunit.eq.0) iunit = i1mach(4)
         do 10 ichar=1,lenmes,72
            last = min0(ichar+71 , lenmes)
            write (iunit,'(1x,a)') messg(ichar:last)
   10    continue
   20 continue
      return
      end
C
C======================================================================C
*DECK ISDCGS
      FUNCTION ISDCGS(N, B, X, NELT, IA, JA, A, ISYM, MATVEC, MSOLVE, 
     $     ITOL, TOL, ITMAX, ITER, ERR, IERR, IUNIT, R, R0, P, Q, U, 
     $     V1, V2, RWORK, IWORK, AK, BK, BNRM, SOLNRM)
C***BEGIN PROLOGUE  ISDCGS
C***REFER TO  DCGS, DSDCGS, DSLUCS
C***DATE WRITTEN   890404   (YYMMDD)
C***REVISION DATE  890404   (YYMMDD)
C***CATEGORY NO.  D2A4
C***KEYWORDS  LIBRARY=SLATEC(SLAP),
C             TYPE=DOUBLE PRECISION(ISDCGS-D),
C             Non-Symmetric Linear system, Sparse,
C             Iterative Precondition, Stop Test
C***AUTHOR  Greenbaum, Anne, Courant Institute
C           Seager, Mark K., (LLNL)
C             Lawrence Livermore National Laboratory
C             PO BOX 808, L-300
C             Livermore, CA 94550 (415) 423-3141
C             seager@lll-crg.llnl.gov
C***PURPOSE  Preconditioned BiConjugate Gradient Stop Test.
C            This routine calculates the stop test for the BiConjugate
C            Gradient iteration scheme.  It returns a nonzero if the
C            error estimate (the type of which is determined by ITOL)
C            is less than the user specified tolerance TOL.
C***DESCRIPTION
C *Usage:
C     INTEGER  N, NELT, IA(NELT), JA(NELT), ISYM, ITOL, ITMAX, ITER
C     INTEGER  IERR, IUNIT, IWORK(USER DEFINED)
C     DOUBLE PRECISION B(N), X(N), A(N), TOL, ERR, R(N), R0(N), P(N)
C     DOUBLE PRECISION Q(N), U(N), V1(N), V2(N)
C     DOUBLE PRECISION RWORK(USER DEFINED), AK, BK, BNRM, SOLNRM
C     EXTERNAL MATVEC, MSOLVE
C
C     IF( ISDCGS(N, B, X, NELT, IA, JA, A, ISYM, MATVEC, MSOLVE, ITOL,
C    $     TOL, ITMAX, ITER, ERR, IERR, IUNIT, R, R0, P, Q, U, V1, 
C    $     V2, RWORK, IWORK, AK, BK, BNRM, SOLNRM) .NE. 0 ) 
C    $     THEN ITERATION DONE
C
C *Arguments:
C N      :IN       Integer
C         Order of the Matrix.
C B      :IN       Double Precision B(N).
C         Right-hand side vector.
C X      :INOUT    Double Precision X(N).
C         On input X is your initial guess for solution vector.
C         On output X is the final approximate solution.
C NELT   :IN       Integer.
C         Number of Non-Zeros stored in A.
C IA     :IN       Integer IA(NELT).
C JA     :IN       Integer JA(NELT).
C A      :IN       Double Precision A(NELT).
C         These arrays contain the matrix data structure for A.
C         It could take any form.  See "LONG DESCRIPTION", in
C         the SLAP routine DCGS for more late breaking details...
C ISYM   :IN       Integer.
C         Flag to indicate symmetric storage format.
C         If ISYM=0, all nonzero entries of the matrix are stored.
C         If ISYM=1, the matrix is symmetric, and only the upper
C         or lower triangle of the matrix is stored.
C MATVEC :EXT      External.
C         Name of a routine which  performs the matrix vector multiply
C         operation  Y = A*X  given A and X.  The  name of  the MATVEC
C         routine must  be declared external  in the  calling program.
C         The calling sequence of MATVEC is:
C             CALL MATVEC( N, X, Y, NELT, IA, JA, A, ISYM )
C         Where N is the number of unknowns, Y is the product A*X upon
C         return,  X is an input  vector.  NELT, IA,  JA,  A and  ISYM
C         define the SLAP matrix data structure: see LONG DESCRIPTION,
C         below.
C MSOLVE :EXT      External.
C         Name of a routine which solves a linear system MZ = R  for Z
C         given R with the preconditioning matrix M (M is supplied via
C         RWORK  and IWORK arrays).   The name  of  the MSOLVE routine
C         must be declared  external  in the  calling   program.   The
C         calling sequence of MSLOVE is:
C             CALL MSOLVE(N, R, Z, NELT, IA, JA, A, ISYM, RWORK, IWORK)
C         Where N is the number of unknowns, R is  the right-hand side
C         vector, and Z is the solution upon return.  NELT,  IA, JA, A
C         and  ISYM define the SLAP  matrix  data structure: see  LONG
C         DESCRIPTION, below.  RWORK is a  double precision array that 
C         can be used
C         to  pass   necessary  preconditioning     information and/or
C         workspace to MSOLVE.  IWORK is an integer work array for the
C         same purpose as RWORK.
C ITOL   :IN       Integer.
C         Flag to indicate type of convergence criterion.
C         If ITOL=1, iteration stops when the 2-norm of the residual 
C         divided by the 2-norm of the right-hand side is less than TOL.
C         This routine must calculate the residual from R = A*X - B.
C         This is un-natural and hence expensive for this type of iter-
C         ative method.  ITOL=2 is *STRONGLY* recommended.
C         If ITOL=2, iteration stops when the 2-norm of M-inv times the 
C         residual divided by the 2-norm of M-inv times the right hand 
C         side is less than tol, where M-inv time a vector is the pre-
C         conditioning step.  This is the *NATURAL* stopping for this 
C         iterative method and is *STRONGLY* recommended.
C         ITOL=11 is often useful for checking and comparing different 
C         routines.  For this case, the user must supply the "exact" 
C         solution or a very accurate approximation (one with an error 
C         much less than TOL) through a common block,
C                     COMMON /SOLBLK/ SOLN(1)
C         if ITOL=11, iteration stops when the 2-norm of the difference 
C         between the iterative approximation and the user-supplied
C         solution divided by the 2-norm of the user-supplied solution 
C         is less than TOL.  Note that this requires the user to set up
C         the "COMMON /SOLBLK/ SOLN(LENGTH)" in the calling routine. 
C         The routine with this declaration should be loaded before the
C         stop test so that the correct length is used by the loader.  
C         This procedure is not standard Fortran and may not work 
C         correctly on your system (although it has worked on every
C         system the authors have tried).  If ITOL is not 11 then this
C         common block is indeed standard Fortran.
C TOL    :IN       Double Precision.
C         Convergence criterion, as described above.
C ITMAX  :IN       Integer.
C         Maximum number of iterations.
C ITER   :OUT      Integer.
C         Number of iterations required to reach convergence, or 
C         ITMAX+1 if convergence criterion could not be achieved in 
C         ITMAX iterations.
C ERR    :OUT      Double Precision.
C         Error estimate of error in final approximate solution, as 
C         defined by ITOL.
C IERR   :OUT      Integer.
C         Error flag.  IERR is set to 3 if ITOL is not on of the 
C         acceptable values, see above. 
C IUNIT  :IN       Integer.
C         Unit number on which to write the error at each iteration, 
C         if this is desired for monitoring convergence.  If unit 
C         number is 0, no writing will occur.
C R      :IN       Double Precision R(N).
C         The residual r = b - Ax.
C R0     :WORK     Double Precision R0(N).
C P      :DUMMY    Double Precision P(N).
C Q      :DUMMY    Double Precision Q(N).
C U      :DUMMY    Double Precision U(N).
C V1     :DUMMY    Double Precision V1(N).
C V2     :WORK     Double Precision V2(N).
C         If ITOL.eq.1 then V2 is used to hold A * X - B on every call.
C         If ITOL.eq.2 then V2 is used to hold M-inv * B on the first
C         call.
C         If ITOL.eq.11 then V2 is used to X - SOLN.
C RWORK  :WORK     Double Precision RWORK(USER DEFINED).
C         Double Precision array that can be used for workspace in 
C         MSOLVE.
C IWORK  :WORK     Integer IWORK(USER DEFINED).
C         Integer array that can be used for workspace in MSOLVE.
C AK     :IN       Double Precision.
C         Current iterate BiConjugate Gradient iteration parameter.
C BK     :IN       Double Precision.
C         Current iterate BiConjugate Gradient iteration parameter.
C BNRM   :INOUT    Double Precision.
C         Norm of the right hand side.  Type of norm depends on ITOL.
C         Calculated only on the first call.
C SOLNRM :INOUT    Double Precision.
C         2-Norm of the true solution, SOLN.  Only computed and used
C         if ITOL = 11.
C
C *Function Return Values:
C       0 : Error estimate (determined by ITOL) is *NOT* less than the 
C           specified tolerance, TOL.  The iteration must continue.
C       1 : Error estimate (determined by ITOL) is less than the 
C           specified tolerance, TOL.  The iteration can be considered
C           complete.
C
C *Precision:           Double Precision
C***REFERENCES  (NONE)
C***ROUTINES CALLED  MATVEC, MSOLVE, DNRM2
C***COMMON BLOCKS    SOLBLK
C***END PROLOGUE  ISDCGS
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, ITOL, ITMAX
      INTEGER ITER, IERR, IUNIT, IWORK(1)
      DOUBLE PRECISION B(N), X(N), A(NELT), TOL, ERR, R(N), R0(N), P(N)
      DOUBLE PRECISION Q(N), U(N), V1(N), V2(N), RWORK(1)
      DOUBLE PRECISION AK, BK, BNRM, SOLNRM
      COMMON /SOLBLK/ SOLN(1)
      EXTERNAL MATVEC, MSOLVE
C         
C***FIRST EXECUTABLE STATEMENT  ISDCGS
      ISDCGS = 0
C         
      IF( ITOL.EQ.1 ) THEN
C         err = ||Residual||/||RightHandSide|| (2-Norms).
         IF(ITER .EQ. 0) BNRM = DNRM21(N, B, 1)
         CALL MATVEC(N, X, V2, NELT, IA, JA, A, ISYM )
         DO 5 I = 1, N
            V2(I) = V2(I) - B(I)
 5       CONTINUE
         ERR = DNRM21(N, V2, 1)/BNRM
      ELSE IF( ITOL.EQ.2 ) THEN
C                  -1              -1
C         err = ||M  Residual||/||M  RightHandSide|| (2-Norms).
         IF(ITER .EQ. 0) THEN
            CALL MSOLVE(N, B, V2, NELT, IA, JA, A, ISYM, RWORK, IWORK)
            BNRM = DNRM21(N, V2, 1)
         ENDIF
         ERR = DNRM21(N, R, 1)/BNRM
      ELSE IF( ITOL.EQ.11 ) THEN
C         err = ||x-TrueSolution||/||TrueSolution|| (2-Norms).
         IF(ITER .EQ. 0) SOLNRM = DNRM21(N, SOLN, 1)
         DO 10 I = 1, N
            V2(I) = X(I) - SOLN(I)
 10      CONTINUE
         ERR = DNRM21(N, V2, 1)/SOLNRM
      ELSE
C
C         If we get here ITOL is not one of the acceptable values.
         ERR = 1.0E10
         IERR = 3
      ENDIF
C         
C         Print the error and Coeficients AK, BK on each step,
C         if desired.
      IF(IUNIT .NE. 0) THEN
         IF( ITER.EQ.0 ) THEN
            WRITE(IUNIT,1000) N, ITOL
         ENDIF
         WRITE(IUNIT,1010) ITER, ERR, AK, BK
      ENDIF
      IF(ERR .LE. TOL) ISDCGS = 1
C         
      RETURN
 1000 FORMAT(' Preconditioned BiConjugate Gradient Squared for ',
     $     'N, ITOL = ',I5, I5,
     $     /' ITER','   Error Estimate','            Alpha',
     $     '             Beta')
 1010 FORMAT(1X,I4,1X,E16.7,1X,E16.7,1X,E16.7)
C------------- LAST LINE OF ISDCGS FOLLOWS ----------------------------
      END
C======================================================================C
