	  ºU  Ï   k820309    ë          13.0        ÖDW                                                                                                           
       MWwaterCoeff_Define.f90 MWWATERCOEFF_DEFINE              MWWATERCOEFF_TYPE MWWATERCOEFF_ASSOCIATED MWWATERCOEFF_DESTROY MWWATERCOEFF_CREATE MWWATERCOEFF_INSPECT MWWATERCOEFF_VALIDRELEASE MWWATERCOEFF_INFO MWWATERCOEFF_DEFINEVERSION MWWATERCOEFF_READFILE MWWATERCOEFF_WRITEFILE i@                      @                              
       FP LONG DOUBLE                      @                              
       SUCCESS FAILURE INFORMATION DISPLAY_MESSAGE                      @                              
       u@LESSTHAN                      @                              
       FILE_OPEN FILE_EXISTS                      @                              
       OPEN_BINARY_FILE WRITEGATTS_BINARY_FILE READGATTS_BINARY_FILE                      @                              
                            @                              
                                                                  #FITCOEFF_1D_EQUAL    #FITCOEFF_2D_EQUAL    #FITCOEFF_3D_EQUAL    #MWWATERCOEFF_EQUAL    %         @    @                                                   #FITCOEFF_1D_EQUAL%ALL 	   #X 
   #Y                  @                            	     ALL           
   @                              
     X              #FITCOEFF_1D_TYPE              
   @                                   X              #FITCOEFF_1D_TYPE    %         @    @                                                    #FITCOEFF_2D_EQUAL%ALL    #X    #Y                  @                                 ALL           
   @                                   x              #FITCOEFF_2D_TYPE              
   @                                   x              #FITCOEFF_2D_TYPE    %         @    @                                                   #FITCOEFF_3D_EQUAL%ALL    #X    #Y                  @                                 ALL           
   @                                                 #FITCOEFF_3D_TYPE              
   @                                                 #FITCOEFF_3D_TYPE    %         @    @X                                                      #X    #Y              
@ @@                                   (             #MWWATERCOEFF_TYPE              
  @@                                   (             #MWWATERCOEFF_TYPE                                                              o #EQUALTO_REAL_SINGLE    #EQUALTO_REAL_DOUBLE !   #EQUALTO_COMPLEX_SINGLE '   #EQUALTO_COMPLEX_DOUBLE ,   %         H    @                                                    #EQUALTO_REAL_SINGLE%ABS    #EQUALTO_REAL_SINGLE%SPACING    #EQUALTO_REAL_SINGLE%MAX    #X    #Y                   @                                 ABS               @                                 SPACING               @                                 MAX           
   @                                  	                
   @                                   	      %         H    @                          !                          #EQUALTO_REAL_DOUBLE%ABS "   #EQUALTO_REAL_DOUBLE%SPACING #   #EQUALTO_REAL_DOUBLE%MAX $   #X %   #Y &                 @                            "     ABS               @                            #     SPACING               @                            $     MAX           
   @                             %     
                
   @                             &     
      %         H    @                          '                          #EQUALTO_COMPLEX_SINGLE%REAL (   #EQUALTO_COMPLEX_SINGLE%AIMAG )   #X *   #Y +                 @              @             (     REAL               @                            )     AIMAG           
   @                             *                     
   @                             +           %         H    @                          ,                          #EQUALTO_COMPLEX_DOUBLE%REAL -   #EQUALTO_COMPLEX_DOUBLE%AIMAG .   #X /   #Y 0                 @              @             -     REAL               @                            .     AIMAG           
   @                             /                     
   @                             0                                                                 u #FILE_UNIT_EXISTS 1   #FILE_NAME_EXISTS 3   %         @    @                           1                           #FILEID 2             
   @                              2           %         @    @                          3                           #FILENAME 4             
   @                             4                    1                                                       u #FILE_OPEN_BY_UNIT 5   #FILE_OPEN_BY_NAME 7   %         @    @                           5                           #FILEID 6             
   @                              6           %         @    @                          7                           #FILENAME 8             
   @                             8                    1                   @              A                'X                    #IS_ALLOCATED 9   #RELEASE :   #VERSION ;   #DIMENSIONS <   #C =                $                              9                                                                                                                    $                             :                                                                                               1                 $                             ;                                                                                               1                $                             <                              p          p            p                                                  ê                                                          0            $                             =                             
            &                                                             @              A                '                    #IS_ALLOCATED >   #RELEASE ?   #VERSION @   #DIMENSIONS A   #C B                $                              >                                                                                                                    $                             ?                                                                                               1                 $                             @                                                                                               1                $                             A                              p          p            p                                                  ê                                                          0            $                             B                             
            &                   &                   &                                                             @              A           C     '`                   #IS_ALLOCATED D   #RELEASE E   #VERSION F   #N_ANGLES G   #N_FREQUENCIES H   #N_TEMPERATURES I   #N_WIND_SPEEDS J   #ANGLE K   #FREQUENCY L   #TEMPERATURE M   #WIND_SPEED N   #EV O   #EH P                $                              D                                                                                                                    $                             E                                                                                               1                 $                             F                                                                                               1                 $                             G                                                                                               0                 $                             H                                                                                               0                 $                             I                                                                                               0                 $                             J                                                                                               0                $                             K                              
            &                                                       $                             L            h              	   
            &                                                       $                             M            °              
   
            &                                                       $                             N            ø                 
            &                                                       $                             O            @                
            &                   &                   &                   &                                                       $                             P            Ð                
            &                   &                   &                   &                                                             @               A                'x                    #IS_ALLOCATED Q   #RELEASE R   #VERSION S   #DIMENSIONS T   #C U                $                              Q                                                                                                                    $                             R                                                                                               1                 $                             S                                                                                               1                $                             T                              p          p            p                                                  ê                                                          0            $                             U                             
            &                   &                                                            !                            V                                                                          !                            W                                                                          !                            X                                                                          !                            Y                                                       0                 !                            Z                                                      3                 !                            [                                                      1#         @      !                           \                   #DISPLAY_MESSAGE%TRIM ]   #DISPLAY_MESSAGE%PRESENT ^   #ROUTINE_NAME _   #MESSAGE `   #ERROR_STATE a   #MESSAGE_LOG b                 @                            ]     TRIM               @                            ^     PRESENT           
   @                             _                    1           
   @                             `                    1           
   @                              a                     
  @                             b                    1 %         @      !                          c                          #OPEN_BINARY_FILE%PRESENT d   #OPEN_BINARY_FILE%TRIM e   #FILENAME f   #FILEID g   #FOR_OUTPUT h   #NO_CHECK i                 @                            d     PRESENT               @                            e     TRIM           
   @                             f                    1             @                              g                      
  @                              h                     
  @                              i           %         @      !                          j                          #WRITEGATTS_BINARY_FILE%PRESENT k   #FID l   #WRITE_MODULE m   #CREATED_ON n   #TITLE o   #HISTORY p   #COMMENT q                 @                            k     PRESENT           
   @                              l                     
  @                             m                    1           
  @                             n                    1           
  @                             o                    1           
  @                             p                    1           
  @                             q                    1 %         @      !                          r                           #FID s   #WRITE_MODULE t   #CREATED_ON u   #TITLE v   #HISTORY w   #COMMENT x             
   @                              s                       @                             t                     1             @                             u                     1             @                             v                     1             @                             w                     1             @                             x                     1 %         @                               y                           #SELF z             
   @                              z     `             #MWWATERLUT_TYPE C   #         @                                  {                   #MWWATERLUT_INSPECT%PRESENT |   #SELF }   #PAUSE ~                 @                            |     PRESENT           
   @                              }     `             #MWWATERLUT_TYPE C             
  @                              ~           %         @                                                          #MWWATERLUT_READFILE%PRESENT    #MWWATERLUT_READFILE%TRIM    #MWWATERLUT    #FILENAME    #NO_CLOSE    #QUIET    #TITLE    #HISTORY    #COMMENT    #DEBUG                  @                                 PRESENT               @                                 TRIM             @                                   `              #MWWATERLUT_TYPE C             
   @                                                 1           
  @                                                   
  @                                                     @                                                  1             @                                                  1             @                                                  1           
  @                                         %         @                                                          #MWWATERLUT_WRITEFILE%PRESENT    #MWWATERLUT_WRITEFILE%TRIM    #MWWATERLUT    #FILENAME    #NO_CLOSE    #QUIET    #TITLE    #HISTORY    #COMMENT    #DEBUG                  @                                 PRESENT               @                                 TRIM           
   @                                   `             #MWWATERLUT_TYPE C             
   @                                                 1           
  @                                                   
  @                                                   
  @                                                 1           
  @                                                 1           
  @                                                 1           
  @                                                           @                               '(             
      #IS_ALLOCATED    #RELEASE    #VERSION    #FCCOEFF    #FRCOEFF    #RCCOEFF    #AZCOEFF    #SSCCOEFF    #LSCCOEFF    #LUT                 $                                                                                                                                                  $                                                                                                                            1                 $                                                                                                                            1                  $                                   X                     #FITCOEFF_1D_TYPE                  $                                   X       h              #FITCOEFF_1D_TYPE                  $                                          À              #FITCOEFF_3D_TYPE                  $                                          P             #FITCOEFF_3D_TYPE                  $                                   X       à             #FITCOEFF_1D_TYPE                  $                                          8      	       #FITCOEFF_3D_TYPE                  $                                   `      È      
       #MWWATERLUT_TYPE C   %         @                                                          #SELF               
   @                                    (             #MWWATERCOEFF_TYPE    #         @                                 ¡                    #SELF ¢             D  @                              ¢     (              #MWWATERCOEFF_TYPE    #         @                                 £                    #SELF ¤             
D  @                              ¤     (              #MWWATERCOEFF_TYPE    #         @                                   ¥                   #MWWATERCOEFF_INSPECT%PRESENT ¦   #SELF §   #PAUSE ¨                 @                            ¦     PRESENT           
@ @@                              §     (             #MWWATERCOEFF_TYPE              
 @@                              ¨           %         @                                ©                           #SELF ª             
   @                              ª     (             #MWWATERCOEFF_TYPE    #         @                                  «                   #MWWATERCOEFF_INFO%LEN_TRIM ¬   #MWWATERCOEFF_INFO%LEN ­   #MWWATERCOEFF_INFO%MIN ®   #MWWATERCOEFF_INFO%ACHAR ¯   #SELF °   #INFO ±                 @                            ¬     LEN_TRIM               @                            ­     LEN               @                            ®     MIN               @                            ¯     ACHAR           
   @                              °     (             #MWWATERCOEFF_TYPE              D @@                             ±                     1 #         @                                   ²                    #ID ³             D  @                             ³                     1 %         @                                 ´                          #MWWATERCOEFF_READFILE%TRIM µ   #MWWATERCOEFF_READFILE%PRESENT ¶   #MWWATERCOEFF ·   #FILENAME ¸   #NO_CLOSE ¹   #QUIET º   #TITLE »   #HISTORY ¼   #COMMENT ½   #DEBUG ¾                 @                            µ     TRIM               @                            ¶     PRESENT           D @@                              ·     (              #MWWATERCOEFF_TYPE              
@ @@                             ¸                    1           
 @@                              ¹                     
B @@                              º                     F @@                             »                     1           F @@                             ¼                     1           F @@                             ½                     1           
B @@                              ¾           %         @                                 ¿                          #MWWATERCOEFF_WRITEFILE%TRIM À   #MWWATERCOEFF_WRITEFILE%PRESENT Á   #MWWATERCOEFF Â   #FILENAME Ã   #NO_CLOSE Ä   #QUIET Å   #TITLE Æ   #HISTORY Ç   #COMMENT È   #DEBUG É                 @                            À     TRIM               @                            Á     PRESENT           
@ @@                              Â     (             #MWWATERCOEFF_TYPE              
@ @@                             Ã                    1           
 @@                              Ä                     
B @@                              Å                     
 @@                             Æ                    1           
 @@                             Ç                    1           
 @@                             È                    1           
B @@                              É                  4      fn#fn )   Ô   ñ   b   uapp(MWWATERCOEFF_DEFINE    Å  O   J  TYPE_KINDS       l   J  MESSAGE_HANDLER &     K   J   COMPARE_FLOAT_NUMBERS    Ë  V   J  FILE_UTILITY $   !  ~   J  BINARY_FILE_UTILITY       @   J   FITCOEFF_DEFINE "   ß  @   J   MWWATERLUT_DEFINE            i@ 2   ¼  y      FITCOEFF_1D_EQUAL+FITCOEFF_DEFINE :   5  <      FITCOEFF_1D_EQUAL%ALL+FITCOEFF_DEFINE=ALL 4   q  ^   e   FITCOEFF_1D_EQUAL%X+FITCOEFF_DEFINE 4   Ï  ^   e   FITCOEFF_1D_EQUAL%Y+FITCOEFF_DEFINE 2   -  y      FITCOEFF_2D_EQUAL+FITCOEFF_DEFINE :   ¦  <      FITCOEFF_2D_EQUAL%ALL+FITCOEFF_DEFINE=ALL 4   â  ^   e   FITCOEFF_2D_EQUAL%X+FITCOEFF_DEFINE 4   @  ^   e   FITCOEFF_2D_EQUAL%Y+FITCOEFF_DEFINE 2     y      FITCOEFF_3D_EQUAL+FITCOEFF_DEFINE :     <      FITCOEFF_3D_EQUAL%ALL+FITCOEFF_DEFINE=ALL 4   S  ^   e   FITCOEFF_3D_EQUAL%X+FITCOEFF_DEFINE 4   ±  ^   e   FITCOEFF_3D_EQUAL%Y+FITCOEFF_DEFINE #   	  ^      MWWATERCOEFF_EQUAL %   m	  _   a   MWWATERCOEFF_EQUAL%X %   Ì	  _   a   MWWATERCOEFF_EQUAL%Y 0   +
  ª      u@EQUALTO+COMPARE_FLOAT_NUMBERS :   Õ
  ¹      EQUALTO_REAL_SINGLE+COMPARE_FLOAT_NUMBERS B     <      EQUALTO_REAL_SINGLE%ABS+COMPARE_FLOAT_NUMBERS=ABS J   Ê  @      EQUALTO_REAL_SINGLE%SPACING+COMPARE_FLOAT_NUMBERS=SPACING B   
  <      EQUALTO_REAL_SINGLE%MAX+COMPARE_FLOAT_NUMBERS=MAX <   F  @   e   EQUALTO_REAL_SINGLE%X+COMPARE_FLOAT_NUMBERS <     @   e   EQUALTO_REAL_SINGLE%Y+COMPARE_FLOAT_NUMBERS :   Æ  ¹      EQUALTO_REAL_DOUBLE+COMPARE_FLOAT_NUMBERS B     <      EQUALTO_REAL_DOUBLE%ABS+COMPARE_FLOAT_NUMBERS=ABS J   »  @      EQUALTO_REAL_DOUBLE%SPACING+COMPARE_FLOAT_NUMBERS=SPACING B   û  <      EQUALTO_REAL_DOUBLE%MAX+COMPARE_FLOAT_NUMBERS=MAX <   7  @   e   EQUALTO_REAL_DOUBLE%X+COMPARE_FLOAT_NUMBERS <   w  @   e   EQUALTO_REAL_DOUBLE%Y+COMPARE_FLOAT_NUMBERS =   ·  ¡      EQUALTO_COMPLEX_SINGLE+COMPARE_FLOAT_NUMBERS G   X  =      EQUALTO_COMPLEX_SINGLE%REAL+COMPARE_FLOAT_NUMBERS=REAL I     >      EQUALTO_COMPLEX_SINGLE%AIMAG+COMPARE_FLOAT_NUMBERS=AIMAG ?   Ó  @   e   EQUALTO_COMPLEX_SINGLE%X+COMPARE_FLOAT_NUMBERS ?     @   e   EQUALTO_COMPLEX_SINGLE%Y+COMPARE_FLOAT_NUMBERS =   S  ¡      EQUALTO_COMPLEX_DOUBLE+COMPARE_FLOAT_NUMBERS G   ô  =      EQUALTO_COMPLEX_DOUBLE%REAL+COMPARE_FLOAT_NUMBERS=REAL I   1  >      EQUALTO_COMPLEX_DOUBLE%AIMAG+COMPARE_FLOAT_NUMBERS=AIMAG ?   o  @   e   EQUALTO_COMPLEX_DOUBLE%X+COMPARE_FLOAT_NUMBERS ?   ¯  @   e   EQUALTO_COMPLEX_DOUBLE%Y+COMPARE_FLOAT_NUMBERS -   ï  l       gen@FILE_EXISTS+FILE_UTILITY .   [  \      FILE_UNIT_EXISTS+FILE_UTILITY 5   ·  @   e   FILE_UNIT_EXISTS%FILEID+FILE_UTILITY .   ÷  ^      FILE_NAME_EXISTS+FILE_UTILITY 7   U  L   e   FILE_NAME_EXISTS%FILENAME+FILE_UTILITY +   ¡  n       gen@FILE_OPEN+FILE_UTILITY /     \      FILE_OPEN_BY_UNIT+FILE_UTILITY 6   k  @   e   FILE_OPEN_BY_UNIT%FILEID+FILE_UTILITY /   «  ^      FILE_OPEN_BY_NAME+FILE_UTILITY 8   	  L   e   FILE_OPEN_BY_NAME%FILENAME+FILE_UTILITY 1   U         FITCOEFF_1D_TYPE+FITCOEFF_DEFINE >   è  ¤   a   FITCOEFF_1D_TYPE%IS_ALLOCATED+FITCOEFF_DEFINE 9     ¥   a   FITCOEFF_1D_TYPE%RELEASE+FITCOEFF_DEFINE 9   1  ¥   a   FITCOEFF_1D_TYPE%VERSION+FITCOEFF_DEFINE <   Ö  ù   a   FITCOEFF_1D_TYPE%DIMENSIONS+FITCOEFF_DEFINE 3   Ï     a   FITCOEFF_1D_TYPE%C+FITCOEFF_DEFINE 1   c         FITCOEFF_3D_TYPE+FITCOEFF_DEFINE >   ö  ¤   a   FITCOEFF_3D_TYPE%IS_ALLOCATED+FITCOEFF_DEFINE 9     ¥   a   FITCOEFF_3D_TYPE%RELEASE+FITCOEFF_DEFINE 9   ?  ¥   a   FITCOEFF_3D_TYPE%VERSION+FITCOEFF_DEFINE <   ä  ù   a   FITCOEFF_3D_TYPE%DIMENSIONS+FITCOEFF_DEFINE 3   Ý  Ä   a   FITCOEFF_3D_TYPE%C+FITCOEFF_DEFINE 2   ¡        MWWATERLUT_TYPE+MWWATERLUT_DEFINE ?   °  ¤   a   MWWATERLUT_TYPE%IS_ALLOCATED+MWWATERLUT_DEFINE :   T  ¥   a   MWWATERLUT_TYPE%RELEASE+MWWATERLUT_DEFINE :   ù  ¥   a   MWWATERLUT_TYPE%VERSION+MWWATERLUT_DEFINE ;      ¥   a   MWWATERLUT_TYPE%N_ANGLES+MWWATERLUT_DEFINE @   C!  ¥   a   MWWATERLUT_TYPE%N_FREQUENCIES+MWWATERLUT_DEFINE A   è!  ¥   a   MWWATERLUT_TYPE%N_TEMPERATURES+MWWATERLUT_DEFINE @   "  ¥   a   MWWATERLUT_TYPE%N_WIND_SPEEDS+MWWATERLUT_DEFINE 8   2#     a   MWWATERLUT_TYPE%ANGLE+MWWATERLUT_DEFINE <   Æ#     a   MWWATERLUT_TYPE%FREQUENCY+MWWATERLUT_DEFINE >   Z$     a   MWWATERLUT_TYPE%TEMPERATURE+MWWATERLUT_DEFINE =   î$     a   MWWATERLUT_TYPE%WIND_SPEED+MWWATERLUT_DEFINE 5   %  Ü   a   MWWATERLUT_TYPE%EV+MWWATERLUT_DEFINE 5   ^&  Ü   a   MWWATERLUT_TYPE%EH+MWWATERLUT_DEFINE 1   :'         FITCOEFF_2D_TYPE+FITCOEFF_DEFINE >   Í'  ¤   a   FITCOEFF_2D_TYPE%IS_ALLOCATED+FITCOEFF_DEFINE 9   q(  ¥   a   FITCOEFF_2D_TYPE%RELEASE+FITCOEFF_DEFINE 9   )  ¥   a   FITCOEFF_2D_TYPE%VERSION+FITCOEFF_DEFINE <   »)  ù   a   FITCOEFF_2D_TYPE%DIMENSIONS+FITCOEFF_DEFINE 3   ´*  ¬   a   FITCOEFF_2D_TYPE%C+FITCOEFF_DEFINE "   `+  p       DOUBLE+TYPE_KINDS    Ð+  p       FP+TYPE_KINDS     @,  p       LONG+TYPE_KINDS (   °,  q       SUCCESS+MESSAGE_HANDLER (   !-  q       FAILURE+MESSAGE_HANDLER ,   -  q       INFORMATION+MESSAGE_HANDLER 0   .  À       DISPLAY_MESSAGE+MESSAGE_HANDLER :   Ã.  =      DISPLAY_MESSAGE%TRIM+MESSAGE_HANDLER=TRIM @    /  @      DISPLAY_MESSAGE%PRESENT+MESSAGE_HANDLER=PRESENT =   @/  L   e   DISPLAY_MESSAGE%ROUTINE_NAME+MESSAGE_HANDLER 8   /  L   e   DISPLAY_MESSAGE%MESSAGE+MESSAGE_HANDLER <   Ø/  @   e   DISPLAY_MESSAGE%ERROR_STATE+MESSAGE_HANDLER <   0  L   e   DISPLAY_MESSAGE%MESSAGE_LOG+MESSAGE_HANDLER 5   d0  Á       OPEN_BINARY_FILE+BINARY_FILE_UTILITY E   %1  @      OPEN_BINARY_FILE%PRESENT+BINARY_FILE_UTILITY=PRESENT ?   e1  =      OPEN_BINARY_FILE%TRIM+BINARY_FILE_UTILITY=TRIM >   ¢1  L   e   OPEN_BINARY_FILE%FILENAME+BINARY_FILE_UTILITY <   î1  @   e   OPEN_BINARY_FILE%FILEID+BINARY_FILE_UTILITY @   .2  @   e   OPEN_BINARY_FILE%FOR_OUTPUT+BINARY_FILE_UTILITY >   n2  @   e   OPEN_BINARY_FILE%NO_CHECK+BINARY_FILE_UTILITY ;   ®2  Ä       WRITEGATTS_BINARY_FILE+BINARY_FILE_UTILITY K   r3  @      WRITEGATTS_BINARY_FILE%PRESENT+BINARY_FILE_UTILITY=PRESENT ?   ²3  @   e   WRITEGATTS_BINARY_FILE%FID+BINARY_FILE_UTILITY H   ò3  L   e   WRITEGATTS_BINARY_FILE%WRITE_MODULE+BINARY_FILE_UTILITY F   >4  L   e   WRITEGATTS_BINARY_FILE%CREATED_ON+BINARY_FILE_UTILITY A   4  L   e   WRITEGATTS_BINARY_FILE%TITLE+BINARY_FILE_UTILITY C   Ö4  L   e   WRITEGATTS_BINARY_FILE%HISTORY+BINARY_FILE_UTILITY C   "5  L   e   WRITEGATTS_BINARY_FILE%COMMENT+BINARY_FILE_UTILITY :   n5          READGATTS_BINARY_FILE+BINARY_FILE_UTILITY >   6  @   e   READGATTS_BINARY_FILE%FID+BINARY_FILE_UTILITY G   N6  L   e   READGATTS_BINARY_FILE%WRITE_MODULE+BINARY_FILE_UTILITY E   6  L   e   READGATTS_BINARY_FILE%CREATED_ON+BINARY_FILE_UTILITY @   æ6  L   e   READGATTS_BINARY_FILE%TITLE+BINARY_FILE_UTILITY B   27  L   e   READGATTS_BINARY_FILE%HISTORY+BINARY_FILE_UTILITY B   ~7  L   e   READGATTS_BINARY_FILE%COMMENT+BINARY_FILE_UTILITY 8   Ê7  Z       MWWATERLUT_ASSOCIATED+MWWATERLUT_DEFINE =   $8  ]   e   MWWATERLUT_ASSOCIATED%SELF+MWWATERLUT_DEFINE 5   8  }       MWWATERLUT_INSPECT+MWWATERLUT_DEFINE E   þ8  @      MWWATERLUT_INSPECT%PRESENT+MWWATERLUT_DEFINE=PRESENT :   >9  ]   e   MWWATERLUT_INSPECT%SELF+MWWATERLUT_DEFINE ;   9  @   e   MWWATERLUT_INSPECT%PAUSE+MWWATERLUT_DEFINE 6   Û9  ö       MWWATERLUT_READFILE+MWWATERLUT_DEFINE F   Ñ:  @      MWWATERLUT_READFILE%PRESENT+MWWATERLUT_DEFINE=PRESENT @   ;  =      MWWATERLUT_READFILE%TRIM+MWWATERLUT_DEFINE=TRIM A   N;  ]   e   MWWATERLUT_READFILE%MWWATERLUT+MWWATERLUT_DEFINE ?   «;  L   e   MWWATERLUT_READFILE%FILENAME+MWWATERLUT_DEFINE ?   ÷;  @   e   MWWATERLUT_READFILE%NO_CLOSE+MWWATERLUT_DEFINE <   7<  @   e   MWWATERLUT_READFILE%QUIET+MWWATERLUT_DEFINE <   w<  L   e   MWWATERLUT_READFILE%TITLE+MWWATERLUT_DEFINE >   Ã<  L   e   MWWATERLUT_READFILE%HISTORY+MWWATERLUT_DEFINE >   =  L   e   MWWATERLUT_READFILE%COMMENT+MWWATERLUT_DEFINE <   [=  @   e   MWWATERLUT_READFILE%DEBUG+MWWATERLUT_DEFINE 7   =  ø       MWWATERLUT_WRITEFILE+MWWATERLUT_DEFINE G   >  @      MWWATERLUT_WRITEFILE%PRESENT+MWWATERLUT_DEFINE=PRESENT A   Ó>  =      MWWATERLUT_WRITEFILE%TRIM+MWWATERLUT_DEFINE=TRIM B   ?  ]   e   MWWATERLUT_WRITEFILE%MWWATERLUT+MWWATERLUT_DEFINE @   m?  L   e   MWWATERLUT_WRITEFILE%FILENAME+MWWATERLUT_DEFINE @   ¹?  @   e   MWWATERLUT_WRITEFILE%NO_CLOSE+MWWATERLUT_DEFINE =   ù?  @   e   MWWATERLUT_WRITEFILE%QUIET+MWWATERLUT_DEFINE =   9@  L   e   MWWATERLUT_WRITEFILE%TITLE+MWWATERLUT_DEFINE ?   @  L   e   MWWATERLUT_WRITEFILE%HISTORY+MWWATERLUT_DEFINE ?   Ñ@  L   e   MWWATERLUT_WRITEFILE%COMMENT+MWWATERLUT_DEFINE =   A  @   e   MWWATERLUT_WRITEFILE%DEBUG+MWWATERLUT_DEFINE "   ]A  Õ       MWWATERCOEFF_TYPE /   2B  ¤   a   MWWATERCOEFF_TYPE%IS_ALLOCATED *   ÖB  ¥   a   MWWATERCOEFF_TYPE%RELEASE *   {C  ¥   a   MWWATERCOEFF_TYPE%VERSION *    D  f   a   MWWATERCOEFF_TYPE%FCCOEFF *   D  f   a   MWWATERCOEFF_TYPE%FRCOEFF *   ìD  f   a   MWWATERCOEFF_TYPE%RCCOEFF *   RE  f   a   MWWATERCOEFF_TYPE%AZCOEFF +   ¸E  f   a   MWWATERCOEFF_TYPE%SSCCOEFF +   F  f   a   MWWATERCOEFF_TYPE%LSCCOEFF &   F  e   a   MWWATERCOEFF_TYPE%LUT (   éF  Z       MWWATERCOEFF_ASSOCIATED -   CG  _   a   MWWATERCOEFF_ASSOCIATED%SELF %   ¢G  R       MWWATERCOEFF_DESTROY *   ôG  _   a   MWWATERCOEFF_DESTROY%SELF $   SH  R       MWWATERCOEFF_CREATE )   ¥H  _   a   MWWATERCOEFF_CREATE%SELF %   I         MWWATERCOEFF_INSPECT -   I  @      MWWATERCOEFF_INSPECT%PRESENT *   ÃI  _   a   MWWATERCOEFF_INSPECT%SELF +   "J  @   a   MWWATERCOEFF_INSPECT%PAUSE *   bJ  Z       MWWATERCOEFF_VALIDRELEASE /   ¼J  _   a   MWWATERCOEFF_VALIDRELEASE%SELF "   K  Ï       MWWATERCOEFF_INFO +   êK  A      MWWATERCOEFF_INFO%LEN_TRIM &   +L  <      MWWATERCOEFF_INFO%LEN &   gL  <      MWWATERCOEFF_INFO%MIN (   £L  >      MWWATERCOEFF_INFO%ACHAR '   áL  _   a   MWWATERCOEFF_INFO%SELF '   @M  L   a   MWWATERCOEFF_INFO%INFO +   M  P       MWWATERCOEFF_DEFINEVERSION .   ÜM  L   a   MWWATERCOEFF_DEFINEVERSION%ID &   (N  ü       MWWATERCOEFF_READFILE +   $O  =      MWWATERCOEFF_READFILE%TRIM .   aO  @      MWWATERCOEFF_READFILE%PRESENT 3   ¡O  _   a   MWWATERCOEFF_READFILE%MWWATERCOEFF /    P  L   a   MWWATERCOEFF_READFILE%FILENAME /   LP  @   a   MWWATERCOEFF_READFILE%NO_CLOSE ,   P  @   a   MWWATERCOEFF_READFILE%QUIET ,   ÌP  L   a   MWWATERCOEFF_READFILE%TITLE .   Q  L   a   MWWATERCOEFF_READFILE%HISTORY .   dQ  L   a   MWWATERCOEFF_READFILE%COMMENT ,   °Q  @   a   MWWATERCOEFF_READFILE%DEBUG '   ðQ  þ       MWWATERCOEFF_WRITEFILE ,   îR  =      MWWATERCOEFF_WRITEFILE%TRIM /   +S  @      MWWATERCOEFF_WRITEFILE%PRESENT 4   kS  _   a   MWWATERCOEFF_WRITEFILE%MWWATERCOEFF 0   ÊS  L   a   MWWATERCOEFF_WRITEFILE%FILENAME 0   T  @   a   MWWATERCOEFF_WRITEFILE%NO_CLOSE -   VT  @   a   MWWATERCOEFF_WRITEFILE%QUIET -   T  L   a   MWWATERCOEFF_WRITEFILE%TITLE /   âT  L   a   MWWATERCOEFF_WRITEFILE%HISTORY /   .U  L   a   MWWATERCOEFF_WRITEFILE%COMMENT -   zU  @   a   MWWATERCOEFF_WRITEFILE%DEBUG 