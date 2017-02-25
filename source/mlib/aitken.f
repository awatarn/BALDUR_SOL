      subroutine Aitken(X,Xbest,K)
C
C  Aitken takes the vector X of length K and repeatedly applies 
C  the Aitken acceleration of convergence algorithm forming a 
C  triangular table of the form A_(Iteration,Level), where
C
C        A_11=X_1   
C        A_21=X_2
C        A_31=X_3    A_12
C        A_41=X_4    A_22
C        A_51=X_5    A_32   A_13  
C        A_61=X_6    A_42   A_23
C        A_71=X_7    A_52   A_33  A_14
C        A_81=X_8    A_62   A_43  A_24
C        A_91=X_9    A_72   A_53  A_34   A_15  etc.
C
C  The extreme right-lowest element is identified as Xbest
C 
C================================================================
      COMPLEX X(90), A(90,90), Xbest(90)
      INTEGER K, ir, ic, LoopEnd
      SAVE A
      IF(K.LE.2)THEN
          write(*,'(A,I3)')'In Aitken, K is too small, K=',K
          RETURN
      ENDIF
      do ir = 1,K
          A(ir,1) = X(ir)
      enddo
      LoopEnd = (K+1)/2 
      do ic = 1,Loopend-1        
         ir = K-2*ic              
         A(ir,ic+1)=A(ir+2,ic) - (A(ir+2,ic)-A(ir+1,ic))**2/
     &              (A(ir+2,ic)-2*A(ir+1,ic)+A(ir,ic))
      enddo
      Xbest(K) = A(K-2*LoopEnd+2,LoopEnd)
      RETURN
      END

