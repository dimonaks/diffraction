      SUBROUTINE SPLINE(N, X, Y, B, C, D)
C
C        Подготовительная программа для сплайн-интерполяции.
C
C        N    INT            IN   число точек сетки
C        X    REAL*8    N    IN   точки сетки
C        Y     - " -    N    IN   значения функции в точках
C      B,C,D   - " -    N    OUT  служебные массивы
C
      IMPLICIT none	
      INTEGER*4 N
      REAL*8 X(N), Y(N), B(N), C(N), D(N)
      INTEGER*4 NM1, IB, I
      REAL*8 T
C
      NM1 = N - 1
      IF (N .LT. 2) RETURN
      IF (N .LT. 3) GOTO 50
C
      D(1) = X(2) - X(1)
      C(2) = (Y(2) - Y(1)) / D(1)
        DO 10 I = 2, NM1
        D(I) = X(I+1) - X(I)
        B(I) = 2.d0*(D(I-1) + D(I))
        C(I+1) = (Y(I+1) - Y(I)) / D(I)
        C(I) = C(I+1) - C(I)
   10   CONTINUE
C
      B(1) = - D(1)
      B(N) = - D(N-1)
      C(1) = 0.D0
      C(N) = 0.D0
      IF (N .EQ. 3) GOTO 15
      C(1) = C(3) / (X(4) - X(2))  -  C(2) / (X(3) - X(1))
      C(N) = C(N-1) / (X(N) - X(N-2))  -  C(N-2) / (X(N-1) - X(N-3))
      C(1) = C(1) * D(1) * D(1) / (X(4) - X(1))
      C(N) = - C(N) * D(N-1) * D(N-1) / (X(N) - X(N-3))
C
   15   DO 20 I = 2, N
        T    = D(I-1) / B(I-1)
        B(I) = B(I) - T*D(I-1)
   20   C(I) = C(I) - T*C(I-1)
C
      C(N) = C(N) / B(N)
        DO 30 IB = 1, NM1
        I = N - IB
   30   C(I) = (C(I) - D(I)*C(I+1))  /  B(I)
C
      B(N) = (Y(N) - Y(NM1)) / D(NM1) + D(NM1)*(C(NM1) + 2.d0*C(N))
        DO 40 I = 1, NM1
        B(I) = (Y(I+1) - Y(I)) / D(I)  -  D(I)*(C(I+1) + 2.d0*C(I))
        D(I) = (C(I+1) - C(I))  /  D(I)
   40   C(I) = 3.d0*C(I)
      C(N) = 3.d0*C(N)
      D(N) = D(N-1)
      RETURN
C
   50 B(1) = (Y(2) - Y(1))  /  (X(2) - X(1))
      C(1) = 0.d0
      D(1) = 0.d0
      B(2) = B(1)
      C(2) = 0.d0
      D(2) = 0.d0
      RETURN
      END
  
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
