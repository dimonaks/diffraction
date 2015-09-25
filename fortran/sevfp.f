      SUBROUTINE SEVFP(N, U, X, Y, B, C, D, S, SP)
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER N, I, J, K
C      REAL*8 U, X(N), Y(N), B(N), C(N), D(N)
      DIMENSION X(N), Y(N), B(N), C(N), D(N)	
      DATA I /1/
      IF (I .GE. N) I = 1
      IF (U .LT. X(I)) GOTO 10
      IF (U .LE. X(I+1)) GOTO 30
   10 I = 1
      J = N + 1
   20 K = (I + J)  /  2
      IF (U .LT. X(K)) J = K
      IF (U .GE. X(K)) I = K
      IF (J .GT. I + 1) GOTO 20
   30 DX   = U - X(I)
      SP   = B(I) + DX*(2.D0*C(I) + 3.D0*DX*D(I))
      S    = Y(I) + DX*(B(I) + DX*(C(I) + DX*D(I)))
      RETURN
      END
 