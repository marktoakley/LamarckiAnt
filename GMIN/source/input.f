C GPL License Info {{{
C   GMIN: A program for finding global minima
C   Copyright (C) 1999-2006 David J. Wales
C   This file is part of GMIN.
C
C   GMIN is free software; you can redistribute it and/or modify
C   it under the terms of the GNU General Public License as published by
C   the Free Software Foundation; either version 2 of the License, or
C   (at your option) any later version.
C
C   GMIN is distributed in the hope that it will be useful,
C   but WITHOUT ANY WARRANTY; without even the implied warranty of
C   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C   GNU General Public License for more details.
C
C   You should have received a copy of the GNU General Public License
C   along with this program; if not, write to the Free Software
C   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
C
C  INPUT  E3   ( 11:41:40 Tuesday 7 September 1993 )
C }}}
 
C  Doxygen Comments {{{
C> \name SUBROUTINE INPUT(LOGICAL END)
C>
C>  Read next input record from unit IR into buffer, and analyse for data
C>  items, null items (enclosed by commas), and comment (enclosed in
C>  parentheses, which may be nested, and occurring where spaces may
C>  occur).
C>
C>  E1  Attempts to handle stream-switching commands in the data.
C>      #include file-name
C>         switches input to be from the file specified;
C>      #revert
C>         (or end-of-file) reverts to the previous file.
C>    Also
C>      #concat "string"
C>         sets the concatenation string; e.g.
C>         #concat "\"
C> \param ITEM    is the number of the last item read from the buffer
C> \param NITEM   is the number of items in the buffer
C> \param CHAR(I) is the Ith character in the buffer
C> \param LOC(I)  is the starting position of the Ith item
C> \param LINE    is the number of lines (records) read
C>
C>  If \b SKIPBL is set to \b TRUE, lines containing no items other than
C>  comment are skipped automatically, so that INPUT will always return
C>  a line with at least one item (possibly null).  If SKIPBL is FALSE
C>  (default) no skipping occurs and the next data line is returned
C>  regardless of what it contains.
C>
C>  If \b CLEAR is set to \b TRUE (default) then an attempt to read more than
C>  NITEM items from a line will cause zeros or blanks to be returned. If
C>  CLEAR is FALSE, the variables in question are left unaltered.
C>
C>  NCR is the number of single characters which have been read from the
C>  current item.  It is always zero unless READCH has been used; if non-
C>  zero, ITEM refers to the current item from which characters are being
C>  read.
C>
C>  NERROR specifies the treatment of errors found while reading numbers:
C>    0  Hard error - print message and stop (default).
C>    1  Soft error - print message and return zero.
C>    2  Soft error - no message, return zero, set NERROR to -1.
C>  If NERROR is set to 2 it is important to test and reset it after reading
C>  a number.
C>
C>  IR is the input stream from which the data are to be read (default 5).
C>
C>  If ECHO is TRUE, the input line will be reflected to the system output
C>  stream.
C>  LAST gives the position of the last non-blank character in the line.
C>
C>  CONCAT is the line concatenation string: if it is found as the last
C>  non-blank character string on a line, the following line is read and
C>  appended to the current line, overwriting the concatenation string.
C>  This procedure is repeated if necessary until a line is found that does
C>  not end with the concatenation string.
C> 
C>  The concatenation string may be modified by changing the DATA statement
C>  above. The integer LC must be set to the number of non-blank characters
C>  in the concatenation string. If it is set to zero, no concatenation
C>  occurs. The concatenation string may also be changed by the #concat
C>  directive in the data file.
C> 
C }}}
      SUBROUTINE INPUT(END, IR)
C  Comments, Declarations {{{
C
C  Read next input record from unit IR into buffer, and analyse for data
C  items, null items (enclosed by commas), and comment (enclosed in
C  parentheses, which may be nested, and occurring where spaces may
C  occur).
C
C  E1  Attempts to handle stream-switching commands in the data.
C      #include file-name
C         switches input to be from the file specified;
C      #revert
C         (or end-of-file) reverts to the previous file.
C    Also
C      #concat "string"
C         sets the concatenation string; e.g.
C         #concat "\"
      CHARACTER CHAR
      COMMON /BUFFER/ CHAR(800)
      INTEGER IR
      LOGICAL SKIPBL, CLEAR, ECHO, CAT
      COMMON /BUFINF/ ITEM, NITEM, LOC(80), LINE, SKIPBL, CLEAR, NCR,
     &                NERROR, ECHO, LAST, CAT

      LOGICAL BLANK, TCOMMA, END
      CHARACTER(LEN=8) CONCAT
      DATA CONCAT /'+++'/, LC /3/
      SAVE CONCAT, LC
 
C     ITEM    is the number of the last item read from the buffer
C     NITEM   is the number of items in the buffer
C     CHAR(I) is the Ith character in the buffer
C     LOC(I)  is the starting position of the Ith item
C     LINE    is the number of lines (records) read
C
C  If SKIPBL is set to TRUE, lines containing no items other than
C
C  comment are skipped automatically, so that INPUT will always return
C  a line with at least one item (possibly null).  If SKIPBL is FALSE
C  (default) no skipping occurs and the next data line is returned
C  regardless of what it contains.
C
C  If CLEAR is set to TRUE (default) then an attempt to read more than
C  NITEM items from a line will cause zeros or blanks to be returned. If
C  CLEAR is FALSE, the variables in question are left unaltered.
C
C  NCR is the number of single characters which have been read from the
C  current item.  It is always zero unless READCH has been used; if non-
C  zero, ITEM refers to the current item from which characters are being
C  read.
C
C  NERROR specifies the treatment of errors found while reading numbers:
C    0  Hard error - print message and stop (default).
C    1  Soft error - print message and return zero.
C    2  Soft error - no message, return zero, set NERROR to -1.
C  If NERROR is set to 2 it is important to test and reset it after reading
C  a number.
 
C  IR is the input stream from which the data are to be read (default 5).
 
C  If ECHO is TRUE, the input line will be reflected to the system output
C  stream.
 
C  LAST gives the position of the last non-blank character in the line.
 
C  CONCAT is the line concatenation string: if it is found as the last
C  non-blank character string on a line, the following line is read and
C  appended to the current line, overwriting the concatenation string.
C  This procedure is repeated if necessary until a line is found that does
C  not end with the concatenation string.
 
C  The concatenation string may be modified by changing the DATA statement
C  above. The integer LC must be set to the number of non-blank characters
C  in the concatenation string. If it is set to zero, no concatenation
C  occurs. The concatenation string may also be changed by the #concat
C  directive in the data file.
 
      INTEGER LEVEL
      SAVE LEVEL, IR0
      CHARACTER(LEN=80) W, F

      CHARACTER SPACE, BRA, KET, COMMA, SQUOTE, DQUOTE, TERM
      DATA LEVEL /0/
      DATA SPACE /' '/, BRA /'('/, KET /')'/, COMMA /','/,
     *    SQUOTE /''''/, DQUOTE /'"'/
C }}}

C {{{ 
      ECHO=.FALSE.

      END=.FALSE.
      CAT=.FALSE.
10    M=1
20    LAST=M+79
      READ (IR,1001,END=900) (CHAR(I), I=M,LAST)
      IF (ECHO) PRINT '(1X,80A1)', (CHAR(I), I=M,LAST)
1001  FORMAT (80A1)
      CAT=.FALSE.
      LINE=LINE+1
C  Find last non-blank character
30    IF (CHAR(LAST) .EQ. SPACE) THEN
        LAST=LAST-1
        IF (LAST .GT. 1) GO TO 30
      ENDIF
C  Look for concatenation string
      IF (LC .GT. 0 .AND. LAST .GE. LC) THEN
        L=LC
        M=LAST
40      IF (CHAR(M) .EQ. CONCAT(L:L)) THEN
          IF (L .EQ. 1) THEN
            CAT=.TRUE.
            GO TO 20
          ENDIF
          L=L-1
          M=M-1
          GO TO 40
        ENDIF
      ENDIF
      CAT=.FALSE.

C  Logical line assembled. First look for input directives
      L=0
50    L=L+1
      IF (CHAR(L) .EQ. SPACE .AND. L .LT. LAST) GO TO 50
      IF (CHAR(L) .NE. '#') GO TO 70
      W=' '
      I=0
52    I=I+1
      W(I:I)=CHAR(L)
      L=L+1
      IF (CHAR(L) .NE. SPACE) GO TO 52
      CALL MYUPCASE(W)
60    L=L+1
      F=' '
      IF (L .GT. LAST) GO TO 67
      IF (CHAR(L) .EQ. SPACE) GO TO 60
      IF (CHAR(L) .EQ. SQUOTE .OR. CHAR(L) .EQ. DQUOTE) THEN
        TERM=CHAR(L)
        L=L+1
      ELSE
        TERM=SPACE
      ENDIF
      I=0
65    IF (L .GT. LAST .OR. CHAR(L) .EQ. TERM) GO TO 67
      I=I+1
      F(I:I)=CHAR(L)
      L=L+1
      GO TO 65
67    IF (W .EQ. '#INCLUDE') THEN
        IF (F .EQ. ' ') THEN
           PRINT '(/A)', 'No filename given in #include directive'
           STOP
         ENDIF
        IF (LEVEL .EQ. 0) IR0=IR
        LEVEL=LEVEL+1
        IR=90+LEVEL
        OPEN (IR,FILE=F,STATUS='UNKNOWN')
      ELSE IF (W .EQ. '#CONCAT') THEN
        CONCAT=F
        LC=I
      ELSE IF (W .EQ. '#REVERT') THEN
        CLOSE(IR)
        LEVEL=LEVEL-1
        IF (LEVEL .EQ. 0) THEN
          IR=IR0
        ELSE
          IR=90+LEVEL
        ENDIF
      ELSE IF (W .EQ. '#ECHO') THEN
        ECHO=.TRUE.
         IF (F .EQ. 'OFF') ECHO=.FALSE.
      ELSE
        GO TO 70
      ENDIF
      GO TO 10

C  Analyse input 
70    ITEM=0
      NITEM=0
      NCR=0
      L=0
 
80    TCOMMA=.TRUE.
90    BLANK=.TRUE.
      TERM=SPACE
 
100   L=L+1
      IF (L .GT. LAST) GO TO 800
      IF (BLANK) GO TO 200
 
C  Reading through an item, seeking terminator.
C  (Next line of code suppressed: comment not allowed within an item)
C     IF (CHAR(L) .EQ. BRA) GO TO 400
      TCOMMA=CHAR(L) .EQ. COMMA .AND. TERM .EQ. SPACE
      BLANK=TCOMMA .OR. CHAR(L) .EQ. TERM
      GO TO 100
 
C  Looking for next item
200   TERM=SPACE
      BLANK=CHAR(L) .EQ. SPACE
      IF (BLANK) GO TO 100
      IF (CHAR(L) .EQ. BRA) GO TO 400
 
C  Item found
      NITEM=NITEM+1
      LOC(NITEM)=L
 
C  Quoted string?
      IF (CHAR(L) .EQ. SQUOTE .OR. CHAR(L) .EQ. DQUOTE) TERM=CHAR(L)
 
C  Null item?
      IF (CHAR(L) .NE. COMMA) GO TO 100
      IF (TCOMMA) LOC(NITEM)=0
      GO TO 80
 
C  Comment (enclosed in parentheses, possibly nested)
400   NEST=1
410   L=L+1
      IF (L .GT. LAST) GO TO 800
      IF (CHAR(L) .EQ. BRA) NEST=NEST+1
      IF (CHAR(L) .EQ. KET) NEST=NEST-1
      IF (NEST .GT. 0) GO TO 410
      GO TO 90
 
C  End of card -- was it blank?
800   IF (SKIPBL .AND. NITEM .EQ. 0) GO TO 10
      RETURN
 
C  End of data
900   IF (CAT) THEN
        PRINT '(A)', 'Apparently concatenating at end-of-file'
        CALL REPORT('Unexpected end of data file',.TRUE.)
      ENDIF
      IF (LEVEL .GT. 0) THEN
C  Revert to previous input
        CLOSE(IR)
        LEVEL=LEVEL-1
        IF (LEVEL .EQ. 0) THEN
          IR=IR0
        ELSE
          IR=90+LEVEL
        ENDIF
        GO TO 10
      ENDIF
      DO L=1,LAST
        CHAR(L)=SPACE
      END DO
      ITEM=0
      NITEM=-1
      END=.TRUE.
      RETURN
C }}} 
      END
 
C-----------------------------------------------------------------------
C> \name BLOCK DATA INBLK 
      BLOCK DATA INBLK

      CHARACTER CHAR
      COMMON /BUFFER/ CHAR(800)
      LOGICAL SKIPBL, CLEAR, ECHO, CAT
      COMMON /BUFINF/ ITEM, NITEM, LOC(80), LINE, SKIPBL, CLEAR, NCR,
     &                NERROR, ECHO, LAST, CAT

      DATA ITEM /0/, NITEM /0/, CHAR /800*' '/, LOC /80*0/, LINE /0/
      DATA SKIPBL /.FALSE./, CLEAR /.TRUE./, NCR /0/, NERROR /0/
      DATA ECHO /.FALSE./
      END
 
C-----------------------------------------------------------------------
 
C     SUBROUTINE STREAM(N)
 
C  Set the input stream for subsequent data to be N.
 
C     LOGICAL SKIPBL, CLEAR, ECHO, CAT
C     COMMON /BUFINF/ ITEM, NITEM, LOC(80), LINE, SKIPBL, CLEAR, NCR,
C    &                NERROR, ECHO, LAST, CAT
C
C     IR=N
C     RETURN
C     END
 
C-----------------------------------------------------------------------
 
      SUBROUTINE READF(A)
C> \name READF
C> \brief Read the next item from the buffer as a real number
C> \param A (DP)
C>
C Declarations {{{ 
      INTEGER P, STATE
      LOGICAL NULL, XNULL
      DOUBLE PRECISION A, AD, B, TEN, SIGN
      CHARACTER CHAR
      COMMON /BUFFER/ CHAR(800)
      LOGICAL SKIPBL, CLEAR, ECHO, CAT
      COMMON /BUFINF/ ITEM, NITEM, LOC(80), LINE, SKIPBL, CLEAR, NCR,
     &                NERROR, ECHO, LAST, CAT

      CHARACTER C, PLUS, MINUS, DOT, SPACE, STAR, D, E, COMMA,
     &    S1, S2, DIGIT(10)
      DATA DIGIT /'1', '2', '3', '4', '5', '6', '7', '8', '9', '0'/
      DATA PLUS /'+'/, MINUS /'-'/, DOT /'.'/, SPACE /' '/, STAR /'*'/,
     *    D /'D'/, E /'E'/, COMMA /','/
      DATA TEN /10D0/
C }}}

C Subroutine body {{{ 
      IF (CLEAR) A=0D0
 
C  If the item is null, or if no item remains, A is unchanged
      IF (ITEM .GE. NITEM) GO TO 110
      ITEM=ITEM+1
      NCR=0
      P=LOC(ITEM)
      IF (P .EQ. 0) GO TO 110
 
      B=0D0
      SIGN=1D0
      NULL=.TRUE.
      XNULL=.FALSE.
      KXIMP=0
      KX=0
      KXSIGN=1
      STATE=1
      GO TO 12
 
10    P=P+1
      IF (P .GT. LAST) GO TO 100
 
C  Terminators
12    C=CHAR(P)
      IF (C .EQ. SPACE .OR. C .EQ. COMMA) GO TO 100
 
C  Digits
      DO 20 I=1,10
      IF (C .EQ. DIGIT(I)) GO TO (70,72,71,80,81), STATE
20    CONTINUE
 
C  State numbers:  1  nothing read yet
C                  2  reading digits of integral part
C                  3  reading decimal fraction
C                  4  expecting exponent
C                  5  reading digits of exponent
 
C  Number control characters -- branch on state
 
      IF (C .EQ. DOT) GO TO (30,30,99,99,99), STATE
      IF (C .EQ. PLUS) GO TO (42,52,52,52,99), STATE
      IF (C .EQ. MINUS) GO TO (41,51,51,51,99), STATE
      IF (C .EQ. D .OR. C .EQ. E) GO TO (99,60,60,99,99), STATE
 
C  Error
99    IF (NERROR .LE. 1) THEN
        WRITE (6,'(1X,A,A,A)') 'Invalid character "', C,
     &                             '" while reading number.         '
        I2=MIN(LAST,P+20)
        I1=MAX(1,I2-70)
        S1=' '
        IF (I1 .GT. 1) S1='...'
        S2=' '
        IF (I2 .LT. LAST) S2='...'
        PRINT '(1X,A/1X,A,1X,80A,1X,A)',
     &           'Current input line is:', S1, (CHAR(I), I=I1,I2), S2
        PRINT '(5X,80A1)', (' ', I=I1,P), '*'
      ENDIF
      IF (NERROR .LE. 0) THEN
        STOP
      ELSE IF (NERROR .EQ. 1) THEN
        WRITE (6,'(1X,A)') 'Item replaced by zero.'
        A=0D0
        RETURN
      ELSE
        NERROR=-1
        A=0D0
        RETURN
      ENDIF
 
 
30    STATE=3
      GO TO 10
 
41    SIGN=-1D0
42    STATE=2
      GO TO 10
 
51    KXSIGN=-1
52    STATE=5
      IF (NULL) GO TO 99
      XNULL=.TRUE.
      GO TO 10
 
60    STATE=4
      GO TO 10
 
 
70    STATE=2
      GO TO 72
71    KXIMP=KXIMP-1
72    AD=I
      IF (I .EQ. 10) AD=0D0
      B=B*TEN+AD
      NULL=.FALSE.
      GO TO 10
 
 
80    STATE=5
81    IF (I .EQ. 10) I=0
      KX=KX*10+I
      XNULL=.FALSE.
      GO TO 10
 
 
100   IF (NULL .OR. XNULL) GO TO 99
      KX=KX*KXSIGN+KXIMP
      A=SIGN*B*(TEN**KX)
 
110   RETURN
C }}} 
      END
 
C-----------------------------------------------------------------------
 
      SUBROUTINE INPUTI(I)
      DOUBLE PRECISION A
      A=I
      CALL READF(A)
      I=A
      RETURN
      END
 
C-----------------------------------------------------------------------
 
C Doxygen comments {{{ 
C> \name READK
C> 
C>  The arrays C, M and N comprise a tree, each node of which consists of
C>  a character C(I) and two integers M(I) and N(I). Initially I=1. The
C>  character C(I) is tested against the current character in the input.
C>  If it matches, the next input character is selected and M(I) is
C>  examined; otherwise N(I) is examined. If a positive value is found,
C>  it is a pointer to the next node to be tested; if negative or zero,
C>  it is minus the key value to be returned in K. If a keyword is
C>  present in the input but does not match a keyword in the tree, then
C>  zero is returned; if a null item is present or if there are no more
C>  items on the line, then -1 is returned. Note that this arrangement
C>  allows the code in the tree to be set up in many ways to suit the
C>  requirements of the problem. Abbreviations may be defined
C>  economically, and the keyword can be required to terminate with a
C>  space or other symbol, or it can be self-terminating.
 
C>  Example:  PRINT and PRT are to return key 19, and must be terminated
C>  by space.  P alone (i.e. not present as the first letter of PRINT or
C>  PRT) is to return 7.  Initial spaces are to be ignored.  The code
C>  table might begin
 
C>        1   2   3   4   5   6   7   8
 
C>   C   ' '  P   R   I   N   T  ' '  ...
C>   M    1   3   4   5   6   7  -19  ...
C>   N    2   8  -7   6   0   0   0   ...
 
C>  Such code tables can be created automatically by the procedure
C>  ENCODE, which produces the DATA statements required to initialize
C>  the array.
C }}}
      SUBROUTINE READK(C,M,N,K)
C Original comments {{{ 
C  The arrays C, M and N comprise a tree, each node of which consists of
C  a character C(I) and two integers M(I) and N(I). Initially I=1. The
C  character C(I) is tested against the current character in the input.
C  If it matches, the next input character is selected and M(I) is
C  examined; otherwise N(I) is examined. If a positive value is found,
C  it is a pointer to the next node to be tested; if negative or zero,
C  it is minus the key value to be returned in K. If a keyword is
C  present in the input but does not match a keyword in the tree, then
C  zero is returned; if a null item is present or if there are no more
C  items on the line, then -1 is returned. Note that this arrangement
C  allows the code in the tree to be set up in many ways to suit the
C  requirements of the problem. Abbreviations may be defined
C  economically, and the keyword can be required to terminate with a
C  space or other symbol, or it can be self-terminating.
 
C  Example:  PRINT and PRT are to return key 19, and must be terminated
C  by space.  P alone (i.e. not present as the first letter of PRINT or
C  PRT) is to return 7.  Initial spaces are to be ignored.  The code
C  table might begin
 
C        1   2   3   4   5   6   7   8
 
C   C   ' '  P   R   I   N   T  ' '  ...
C   M    1   3   4   5   6   7  -19  ...
C   N    2   8  -7   6   0   0   0   ...
 
C  Such code tables can be created automatically by the procedure
C  ENCODE, which produces the DATA statements required to initialize
C  the array.
C }}}
 
C Declarations {{{ 
      INTEGER M(*), N(*), TP, P
      CHARACTER C(*), SPACE
      CHARACTER CHAR
      COMMON /BUFFER/ CHAR(800)
      LOGICAL SKIPBL, CLEAR, ECHO, CAT
      COMMON /BUFINF/ ITEM, NITEM, LOC(80), LINE, SKIPBL, CLEAR, NCR,
     &                NERROR, ECHO, LAST, CAT
      DATA SPACE /' '/
C }}} 

C Subroutine body {{{
      K=-1
      IF (ITEM .GE. NITEM) RETURN
      ITEM=ITEM+1
      NCR=0
      P=LOC(ITEM)
      IF (P .EQ. 0) RETURN
      K=0
      TP=1
      GO TO 20
 
C  Advance character pointer
10    P=P+1
20    IF (P .LE. LAST) GO TO 30
C  End of line
      IF (TP .EQ. 1) K=-1
      RETURN
 
30    CALL MYUPCASE(CHAR(P))
40    IF(CHAR(P) .EQ. C(TP)) GO TO 50
C  No match -- advance tree pointer
      TP=N(TP)
C  Zero pointer -- undecodeable
      IF (TP .EQ. 0) RETURN
C  Positive value -- test same input character again
      IF (TP .GT. 0) GO TO 40
C  Negative fail pointer: the input contains a keyword which is an initial
C  substring of some other keyword.
      GO TO 70
 
C  Matched -- advance tree pointer to next character of command
50    TP=M(TP)
      IF (TP .GT. 0) GO TO 10
 
C  Last character of keyword found.
70    K=-TP
C }}}
      RETURN
      END
 
C-----------------------------------------------------------------------
 
      SUBROUTINE READA(M)
      IMPLICIT NONE
 
C  Read characters and pack them into the character variable M.
C  If the first character is a single or double quote, the string is
C  terminated by a matching quote and the quotes are removed; otherwise
C  it is terminated by space or comma.  If necessary, M is filled
C  out with spaces.
C Declarations {{{ 
      CHARACTER SPACE, BLANK, COMMA, SQUOTE, DQUOTE, TERM
      CHARACTER*(*) M
      CHARACTER CHAR
      INTEGER I, K, NCR, NERROR, ITEM, NITEM, LOC, LAST, L, LINE
      COMMON /BUFFER/ CHAR(800)
      LOGICAL SKIPBL, CLEAR, ECHO, CAT
      COMMON /BUFINF/ ITEM, NITEM, LOC(80), LINE, SKIPBL, CLEAR, NCR,
     &                NERROR, ECHO, LAST, CAT
      LOGICAL QUOTE
      DATA SPACE /' '/, COMMA /','/, SQUOTE /''''/, DQUOTE /'"'/
C }}}

C Subroutine body {{{ 

      IF (ITEM .GE. NITEM .AND. .NOT. CLEAR) RETURN
      K=0
      IF (ITEM .GE. NITEM) GOTO 65
      ITEM=ITEM+1
      NCR=0
      L=LOC(ITEM)
      IF (L .EQ. 0) RETURN
 
      QUOTE=(CHAR(L) .EQ. SQUOTE .OR. CHAR(L) .EQ. DQUOTE)
      IF (.NOT. QUOTE) THEN
         K=K+1
         GOTO 10
      ENDIF
      TERM=CHAR(L)
      QUOTE=.TRUE.
      L=L+1
      GO TO 11
 
10    CONTINUE
C
C  The new Sun compiler refuses to execute the following line - hence rewritten
C
C     K=K+1
C
11    IF (L+K .GT. LAST) GO TO 50

C     IF (.NOT. QUOTE .AND. CHAR(L+K) .NE. SPACE
C    *                .AND. CHAR(L+K) .NE. COMMA) GO TO 10
C     IF (QUOTE .AND. CHAR(L+K) .NE. TERM) GO TO 10

      IF ((.NOT. QUOTE .AND. CHAR(L+K) .NE. SPACE
     *                .AND. CHAR(L+K) .NE. COMMA).OR.
     *    (QUOTE .AND. CHAR(L+K) .NE. TERM)) THEN
         K=K+1
         GOTO 10
      ENDIF
 
50    IF (K .EQ. 0) RETURN
      K=MIN0(K,LEN(M))
      DO I=1,K
        M(I:I)=CHAR(L+I-1)
      END DO
65    DO I=K+1,LEN(M)
        M(I:I)=SPACE
      END DO
C }}}
      RETURN
      END
 
C-----------------------------------------------------------------------
 
      SUBROUTINE READU(M)
      CHARACTER*(*) M
 
      CALL READA(M)
      CALL MYUPCASE(M)
      RETURN

      ENTRY READL(M)
      CALL READA(M)
      CALL LOCASE(M)
      RETURN
      END
 
C-----------------------------------------------------------------------
 
      SUBROUTINE READCH(M)
 
C  Read a single character from the next (or the current) data item.
C  No account is taken of special characters.

C Declarations {{{ 
      CHARACTER CHAR
      COMMON /BUFFER/ CHAR(800)
      LOGICAL SKIPBL, CLEAR, ECHO, CAT
      COMMON /BUFINF/ ITEM, NITEM, LOC(80), LINE, SKIPBL, CLEAR, NCR,
     &                NERROR, ECHO, LAST, CAT
      CHARACTER*(*) M
C }}} 
      M=' '
      IF (ITEM .GE. NITEM .AND. NCR .EQ. 0) RETURN
      IF (NCR .EQ. 0) ITEM=ITEM+1
      L=LOC(ITEM)
      M=CHAR(L+NCR)
      NCR=NCR+1
      RETURN
      END
 
C-----------------------------------------------------------------------
 
!      SUBROUTINE GETF(A)
!C  Read the next item as a double-precision number, reading new data
!C  records if necessary
!C Declarations  {{{
!      DOUBLE PRECISION A
!
!      CHARACTER CHAR
!      COMMON /BUFFER/ CHAR(800)
!      LOGICAL SKIPBL, CLEAR, ECHO, CAT
!      COMMON /BUFINF/ ITEM, NITEM, LOC(80), LINE, SKIPBL, CLEAR, NCR,
!     &                NERROR, ECHO, LAST, CAT
!
!      LOGICAL END
!C }}}
! 
!10    IF (ITEM .LT. NITEM) THEN
!        CALL READF(A)
!        RETURN
!      ENDIF
!      CALL INPUT(END)
!      IF (END) THEN
!        WRITE (6,1001)
!1001    FORMAT ('0End of file while attempting to read a number')
!        STOP
!      ENDIF
!      GO TO 10
! 
!      END
 
C-----------------------------------------------------------------------
 
!      SUBROUTINE GETS(S)
!  Get a single-precision number
!      DOUBLE PRECISION A
!      REAL S
!      A=S
!      CALL GETF(A)
!      S=A
!      RETURN
!      END
 
C-----------------------------------------------------------------------
 
!      SUBROUTINE GETI(I)
!C  Get an integer
!      DOUBLE PRECISION A
!      A=I
!      CALL GETF(A)
!      I=A
!      RETURN
!      END
 
C-----------------------------------------------------------------------
!> \name GETA
!> \brief Get a character string 
!      SUBROUTINE GETA(M)
!  Get a character string
!  Declarations {{{
!      CHARACTER*(*) M
!
!      CHARACTER CHAR
!      COMMON /BUFFER/ CHAR(800)
!      LOGICAL SKIPBL, CLEAR, ECHO, CAT
!      COMMON /BUFINF/ ITEM, NITEM, LOC(80), LINE, SKIPBL, CLEAR, NCR,
!     &                NERROR, ECHO, LAST, CAT
!
! }}}
!      LOGICAL END
!
!10    IF (ITEM .LT. NITEM) THEN
!        CALL READA(M)
!        RETURN
!      ENDIF
!      CALL INPUT(END)
!      IF (END) THEN
!        WRITE (6,1001)
!1001    FORMAT
!     &    ('0End of file while attempting to read a character string')
!        STOP
!      ENDIF
!      GO TO 10
! 
!      END
 
C-----------------------------------------------------------------------
 
      SUBROUTINE READSNGL(S)
C>  \name READSNGL
C>  \brief Read a number from the current record into the single-precision
C>  \brief variable S
C  Read a number from the current record into the single-precision
C  variable S
      DOUBLE PRECISION A
      REAL S
      A=S
      CALL READF(A)
      S=A
      RETURN
      END
 
C-----------------------------------------------------------------------
 
      SUBROUTINE READI(I)
C> \name READI
C> \brief Read an integer from the current record
      DOUBLE PRECISION A
      A=I
      CALL READF(A)
      I=A
      RETURN
      END
 
C-----------------------------------------------------------------------
 
      SUBROUTINE REREAD(K)
C> \name REREAD 
C>
C>  K>0  Reread from item K
C>  K<0  Go back |K| items
C>  K=0  Same as K=-1, i.e. reread last item.
 
      CHARACTER CHAR
      COMMON /BUFFER/ CHAR(800)
      LOGICAL SKIPBL, CLEAR, ECHO, CAT
      COMMON /BUFINF/ ITEM, NITEM, LOC(80), LINE, SKIPBL, CLEAR, NCR,
     &                NERROR, ECHO, LAST, CAT
 
      IF (K .LT. 0) ITEM=ITEM+K
      IF (K .EQ. 0) ITEM=ITEM-1
      IF (K .GT. 0) ITEM=K-1
      IF (ITEM .LT. 0) ITEM=0
      NCR=0
      RETURN
      END
 
C-----------------------------------------------------------------------
 
      SUBROUTINE MYUPCASE(WORD)
      CHARACTER WORD*(*)
 
      CHARACTER UC*26, LC*26
      DATA UC /'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
      DATA LC /'abcdefghijklmnopqrstuvwxyz'/
 
      DO 10 I=1,LEN(WORD)
      K=INDEX(LC,WORD(I:I))
      IF (K .NE. 0) WORD(I:I)=UC(K:K)
10    CONTINUE
      RETURN
 
      ENTRY LOCASE(WORD)
 
      DO 20 I=1,LEN(WORD)
      K=INDEX(UC,WORD(I:I))
      IF (K .NE. 0) WORD(I:I)=LC(K:K)
20    CONTINUE
      RETURN
 
      END
C  REPORT  B1  ( 16:10:41 Thursday 30 April 1992 )

      SUBROUTINE REPORT(C,REFLCT)
C Declarations {{{
      CHARACTER C*(*)
      LOGICAL REFLCT

      CHARACTER BUFF
      COMMON /BUFFER/ BUFF(800)
      LOGICAL SKIPBL, CLEAR, ECHO, CAT
      COMMON /BUFINF/ ITEM, NITEMS, LOC(80), LINE, SKIPBL, CLEAR, NCR,
     &                NERROR, ECHO, LAST, CAT

      INTEGER PRINT
      LOGICAL SWITCH, ONLINE
      COMMON /TEST/ SWITCH(8), PRINT
      EQUIVALENCE (ONLINE, SWITCH(4))

      CHARACTER(LEN=3) S1, S2
C  }}}

      PRINT '(1X,A)', C
      IF (REFLCT) THEN
        L=LOC(ITEM)
        I2=MIN(LAST,L+20)
        I1=MAX(1,I2-70)
        S1=' '
        IF (I1 .GT. 1) S1='...'
        S2=' '
        IF (I2 .LT. LAST) S2='...'
C       PRINT '(1X,4(A,I3))',
C    &      'L =', L, '   LAST =', LAST, '   I1 =', I1, '   I2 =', I2
        PRINT '(1X,A/1X,A,1X,80A,1X,A)',
     &           'Current input line is:', S1, (BUFF(I), I=I1,I2), S2
        PRINT '(4X,80A1)', (' ', I=I1,L), '*'
      ENDIF
C     IF (.NOT. ONLINE) STOP
      END
