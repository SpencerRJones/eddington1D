	  .0  V   k820309    ë          13.0        ±ÖDW                                                                                                           
       CRTM_MoleculeScatter.f90 CRTM_MOLECULESCATTER              CRTM_COMPUTE_MOLECULESCATTER CRTM_COMPUTE_MOLECULESCATTER_TL CRTM_COMPUTE_MOLECULESCATTER_AD                      @                              
       FP                      @                              
       SUCCESS FAILURE DISPLAY_MESSAGE                      @                              
       ZERO                      @                              
       CRTM_ATMOSPHERE_TYPE                      @                              
       CRTM_ATMOPTICS_TYPE                   !@               Á                '                   #IS_ALLOCATED    #MAX_LAYERS    #N_LAYERS 	   #N_ABSORBERS 
   #MAX_CLOUDS    #N_CLOUDS    #MAX_AEROSOLS    #N_AEROSOLS    #N_ADDED_LAYERS    #CLIMATOLOGY    #ABSORBER_ID    #ABSORBER_UNITS    #LEVEL_PRESSURE    #PRESSURE    #TEMPERATURE    #ABSORBER    #CLOUD    #AEROSOL !                $                                                                                                                                                  $                                                                                                                             0                 $                              	                                                                                               0                 $                              
                                                                                               0                 $                                                                                                                             0                 $                                                                                                                             0                 $                                                                                                                             0                 $                                                                                                                             0                 $                                           	                                                                                   0                 $                                   $       
                                                                                   6                $                                          (                             &                                                       $                                          p                             &                                                       $                                         ¸                 
            &                                                       $                                                          
            &                                                       $                                         H                
            &                                                       $                                                         
            &                   &                                                       $                                          ð      ð             #CRTM_CLOUD_TYPE              &                                                          @  !@              E                'ð                    #IS_ALLOCATED    #MAX_LAYERS    #N_LAYERS    #N_ADDED_LAYERS    #TYPE    #EFFECTIVE_RADIUS    #EFFECTIVE_VARIANCE    #WATER_CONTENT                  $                                                                                                                                                  $                                                                                                                             0                 $                                                                                                                             0                 $                                                                                                                             0                 $                                                                                                                              0                $                                                          
            &                                                       $                                         `                 
            &                                                       $                                          ¨                 
            &                                                       $                              !            8      ¨             #CRTM_AEROSOL_TYPE "             &                                                          @  !@              E           "     '¨                    #IS_ALLOCATED #   #MAX_LAYERS $   #N_LAYERS %   #N_ADDED_LAYERS &   #TYPE '   #EFFECTIVE_RADIUS (   #CONCENTRATION )                $                              #                                                                                                                    $                              $                                                                                               0                 $                              %                                                                                               0                 $                              &                                                                                               0                 $                              '                                                                                                0                $                             (                             
            &                                                       $                             )            `                 
            &                                                             !@               A           *     'È                   #IS_ALLOCATED +   #RELEASE ,   #N_LAYERS -   #N_LEGENDRE_TERMS .   #N_PHASE_ELEMENTS /   #MAX_LAYERS 0   #MAX_LEGENDRE_TERMS 1   #MAX_PHASE_ELEMENTS 2   #INCLUDE_SCATTERING 3   #LOFFSET 4   #SCATTERING_OPTICAL_DEPTH 5   #OPTICAL_DEPTH 6   #SINGLE_SCATTER_ALBEDO 7   #ASYMMETRY_FACTOR 8   #DELTA_TRUNCATION 9   #PHASE_COEFFICIENT :                $                              +                                                                                                                    $                              ,                                                                                               4                 $                              -                                                                                               0                 $                              .                                                                                               0                 $                              /                                                                                               0                 $                              0                                                                                               0                 $                              1                                                                                               0                 $                              2                                                                                               0                 $                              3             	                                                                      ÿÿÿÿÿÿÿÿ                         $                              4     $       
                                                                                   0                 $                             5     (         
                                                 
                                 0.0                $                             6            0                 
            &                                                       $                             7            x                 
            &                                                       $                             8            À                 
            &                                                       $                             9                            
            &                                                       $                             :            P                
            &                   &                   &                                                            !                            ;                                                                          !                            <                                                       0                 !                            =                                                      3#         @      !                           >                   #DISPLAY_MESSAGE%TRIM ?   #DISPLAY_MESSAGE%PRESENT @   #ROUTINE_NAME A   #MESSAGE B   #ERROR_STATE C   #MESSAGE_LOG D                 @                            ?     TRIM               @                            @     PRESENT           
   @                             A                    1           
   @                             B                    1           
   @                              C                     
  @                             D                    1                  !                           E     
                
                                 0.0%         @                                 F                           #WAVENUMBER G   #ATMOSPHERE H   #ATMOPTICS I   #MESSAGE_LOG J             
  @@                             G     
                
   @                              H                  #CRTM_ATMOSPHERE_TYPE              
D  @                              I     È              #CRTM_ATMOPTICS_TYPE *             
 @@                             J                    1 %         @                                 K                           #WAVENUMBER L   #ATMOSPHERE_TL M   #ATMOPTICS_TL N   #MESSAGE_LOG O             
  @@                             L     
                
   @                              M                  #CRTM_ATMOSPHERE_TYPE              
D  @                              N     È              #CRTM_ATMOPTICS_TYPE *             
 @@                             O                    1 %         @                                 P                           #WAVENUMBER Q   #ATMOPTICS_AD R   #ATMOSPHERE_AD S   #MESSAGE_LOG T             
  @@                             Q     
                
  @                              R     È              #CRTM_ATMOPTICS_TYPE *             
D  @                              S                   #CRTM_ATMOSPHERE_TYPE              
 @@                             T                    1        6      fn#fn *   Ö   m   b   uapp(CRTM_MOLECULESCATTER    C  C   J  TYPE_KINDS       `   J  MESSAGE_HANDLER     æ  E   J  CRTM_PARAMETERS '   +  U   J  CRTM_ATMOSPHERE_DEFINE &     T   J  CRTM_ATMOPTICS_DEFINE <   Ô  t      CRTM_ATMOSPHERE_TYPE+CRTM_ATMOSPHERE_DEFINE I   H  ¤   a   CRTM_ATMOSPHERE_TYPE%IS_ALLOCATED+CRTM_ATMOSPHERE_DEFINE G   ì  ¥   a   CRTM_ATMOSPHERE_TYPE%MAX_LAYERS+CRTM_ATMOSPHERE_DEFINE E     ¥   a   CRTM_ATMOSPHERE_TYPE%N_LAYERS+CRTM_ATMOSPHERE_DEFINE H   6  ¥   a   CRTM_ATMOSPHERE_TYPE%N_ABSORBERS+CRTM_ATMOSPHERE_DEFINE G   Û  ¥   a   CRTM_ATMOSPHERE_TYPE%MAX_CLOUDS+CRTM_ATMOSPHERE_DEFINE E     ¥   a   CRTM_ATMOSPHERE_TYPE%N_CLOUDS+CRTM_ATMOSPHERE_DEFINE I   %  ¥   a   CRTM_ATMOSPHERE_TYPE%MAX_AEROSOLS+CRTM_ATMOSPHERE_DEFINE G   Ê  ¥   a   CRTM_ATMOSPHERE_TYPE%N_AEROSOLS+CRTM_ATMOSPHERE_DEFINE K   o	  ¥   a   CRTM_ATMOSPHERE_TYPE%N_ADDED_LAYERS+CRTM_ATMOSPHERE_DEFINE H   
  ¥   a   CRTM_ATMOSPHERE_TYPE%CLIMATOLOGY+CRTM_ATMOSPHERE_DEFINE H   ¹
     a   CRTM_ATMOSPHERE_TYPE%ABSORBER_ID+CRTM_ATMOSPHERE_DEFINE K   M     a   CRTM_ATMOSPHERE_TYPE%ABSORBER_UNITS+CRTM_ATMOSPHERE_DEFINE K   á     a   CRTM_ATMOSPHERE_TYPE%LEVEL_PRESSURE+CRTM_ATMOSPHERE_DEFINE E   u     a   CRTM_ATMOSPHERE_TYPE%PRESSURE+CRTM_ATMOSPHERE_DEFINE H   	     a   CRTM_ATMOSPHERE_TYPE%TEMPERATURE+CRTM_ATMOSPHERE_DEFINE E     ¬   a   CRTM_ATMOSPHERE_TYPE%ABSORBER+CRTM_ATMOSPHERE_DEFINE B   I  ©   a   CRTM_ATMOSPHERE_TYPE%CLOUD+CRTM_ATMOSPHERE_DEFINE 2   ò  ß      CRTM_CLOUD_TYPE+CRTM_CLOUD_DEFINE ?   Ñ  ¤   a   CRTM_CLOUD_TYPE%IS_ALLOCATED+CRTM_CLOUD_DEFINE =   u  ¥   a   CRTM_CLOUD_TYPE%MAX_LAYERS+CRTM_CLOUD_DEFINE ;     ¥   a   CRTM_CLOUD_TYPE%N_LAYERS+CRTM_CLOUD_DEFINE A   ¿  ¥   a   CRTM_CLOUD_TYPE%N_ADDED_LAYERS+CRTM_CLOUD_DEFINE 7   d  ¥   a   CRTM_CLOUD_TYPE%TYPE+CRTM_CLOUD_DEFINE C   	     a   CRTM_CLOUD_TYPE%EFFECTIVE_RADIUS+CRTM_CLOUD_DEFINE E        a   CRTM_CLOUD_TYPE%EFFECTIVE_VARIANCE+CRTM_CLOUD_DEFINE @   1     a   CRTM_CLOUD_TYPE%WATER_CONTENT+CRTM_CLOUD_DEFINE D   Å  «   a   CRTM_ATMOSPHERE_TYPE%AEROSOL+CRTM_ATMOSPHERE_DEFINE 6   p  Ç      CRTM_AEROSOL_TYPE+CRTM_AEROSOL_DEFINE C   7  ¤   a   CRTM_AEROSOL_TYPE%IS_ALLOCATED+CRTM_AEROSOL_DEFINE A   Û  ¥   a   CRTM_AEROSOL_TYPE%MAX_LAYERS+CRTM_AEROSOL_DEFINE ?     ¥   a   CRTM_AEROSOL_TYPE%N_LAYERS+CRTM_AEROSOL_DEFINE E   %  ¥   a   CRTM_AEROSOL_TYPE%N_ADDED_LAYERS+CRTM_AEROSOL_DEFINE ;   Ê  ¥   a   CRTM_AEROSOL_TYPE%TYPE+CRTM_AEROSOL_DEFINE G   o     a   CRTM_AEROSOL_TYPE%EFFECTIVE_RADIUS+CRTM_AEROSOL_DEFINE D        a   CRTM_AEROSOL_TYPE%CONCENTRATION+CRTM_AEROSOL_DEFINE :           CRTM_ATMOPTICS_TYPE+CRTM_ATMOPTICS_DEFINE G   4  ¤   a   CRTM_ATMOPTICS_TYPE%IS_ALLOCATED+CRTM_ATMOPTICS_DEFINE B   Ø  ¥   a   CRTM_ATMOPTICS_TYPE%RELEASE+CRTM_ATMOPTICS_DEFINE C   }  ¥   a   CRTM_ATMOPTICS_TYPE%N_LAYERS+CRTM_ATMOPTICS_DEFINE K   "  ¥   a   CRTM_ATMOPTICS_TYPE%N_LEGENDRE_TERMS+CRTM_ATMOPTICS_DEFINE K   Ç  ¥   a   CRTM_ATMOPTICS_TYPE%N_PHASE_ELEMENTS+CRTM_ATMOPTICS_DEFINE E   l  ¥   a   CRTM_ATMOPTICS_TYPE%MAX_LAYERS+CRTM_ATMOPTICS_DEFINE M      ¥   a   CRTM_ATMOPTICS_TYPE%MAX_LEGENDRE_TERMS+CRTM_ATMOPTICS_DEFINE M   ¶   ¥   a   CRTM_ATMOPTICS_TYPE%MAX_PHASE_ELEMENTS+CRTM_ATMOPTICS_DEFINE M   [!  ¤   a   CRTM_ATMOPTICS_TYPE%INCLUDE_SCATTERING+CRTM_ATMOPTICS_DEFINE B   ÿ!  ¥   a   CRTM_ATMOPTICS_TYPE%LOFFSET+CRTM_ATMOPTICS_DEFINE S   ¤"  §   a   CRTM_ATMOPTICS_TYPE%SCATTERING_OPTICAL_DEPTH+CRTM_ATMOPTICS_DEFINE H   K#     a   CRTM_ATMOPTICS_TYPE%OPTICAL_DEPTH+CRTM_ATMOPTICS_DEFINE P   ß#     a   CRTM_ATMOPTICS_TYPE%SINGLE_SCATTER_ALBEDO+CRTM_ATMOPTICS_DEFINE K   s$     a   CRTM_ATMOPTICS_TYPE%ASYMMETRY_FACTOR+CRTM_ATMOPTICS_DEFINE K   %     a   CRTM_ATMOPTICS_TYPE%DELTA_TRUNCATION+CRTM_ATMOPTICS_DEFINE L   %  Ä   a   CRTM_ATMOPTICS_TYPE%PHASE_COEFFICIENT+CRTM_ATMOPTICS_DEFINE    _&  p       FP+TYPE_KINDS (   Ï&  q       SUCCESS+MESSAGE_HANDLER (   @'  q       FAILURE+MESSAGE_HANDLER 0   ±'  À       DISPLAY_MESSAGE+MESSAGE_HANDLER :   q(  =      DISPLAY_MESSAGE%TRIM+MESSAGE_HANDLER=TRIM @   ®(  @      DISPLAY_MESSAGE%PRESENT+MESSAGE_HANDLER=PRESENT =   î(  L   e   DISPLAY_MESSAGE%ROUTINE_NAME+MESSAGE_HANDLER 8   :)  L   e   DISPLAY_MESSAGE%MESSAGE+MESSAGE_HANDLER <   )  @   e   DISPLAY_MESSAGE%ERROR_STATE+MESSAGE_HANDLER <   Æ)  L   e   DISPLAY_MESSAGE%MESSAGE_LOG+MESSAGE_HANDLER %   *  s       ZERO+CRTM_PARAMETERS -   *         CRTM_COMPUTE_MOLECULESCATTER 8   +  @   a   CRTM_COMPUTE_MOLECULESCATTER%WAVENUMBER 8   U+  b   a   CRTM_COMPUTE_MOLECULESCATTER%ATMOSPHERE 7   ·+  a   a   CRTM_COMPUTE_MOLECULESCATTER%ATMOPTICS 9   ,  L   a   CRTM_COMPUTE_MOLECULESCATTER%MESSAGE_LOG 0   d,         CRTM_COMPUTE_MOLECULESCATTER_TL ;   ú,  @   a   CRTM_COMPUTE_MOLECULESCATTER_TL%WAVENUMBER >   :-  b   a   CRTM_COMPUTE_MOLECULESCATTER_TL%ATMOSPHERE_TL =   -  a   a   CRTM_COMPUTE_MOLECULESCATTER_TL%ATMOPTICS_TL <   ý-  L   a   CRTM_COMPUTE_MOLECULESCATTER_TL%MESSAGE_LOG 0   I.         CRTM_COMPUTE_MOLECULESCATTER_AD ;   ß.  @   a   CRTM_COMPUTE_MOLECULESCATTER_AD%WAVENUMBER =   /  a   a   CRTM_COMPUTE_MOLECULESCATTER_AD%ATMOPTICS_AD >   /  b   a   CRTM_COMPUTE_MOLECULESCATTER_AD%ATMOSPHERE_AD <   â/  L   a   CRTM_COMPUTE_MOLECULESCATTER_AD%MESSAGE_LOG 