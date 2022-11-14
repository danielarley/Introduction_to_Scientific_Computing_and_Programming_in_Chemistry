PROGRAM geom
IMPLICIT NONE 

! Define as variáveis e vetores que serão usadas ao longo do programa
REAL, DIMENSION(100) :: X, Y, Z
REAL, DIMENSION(3) :: ATOMO1, ATOMO2
REAL :: DISTANCIA, RAIO1, RAIO2, VDW_SOMA, PARAMETRO
CHARACTER(LEN=2), DIMENSION(100) :: SIMBOLOS
INTEGER :: I, NUM_ATOMOS, J, A
CHARACTER(LEN=2), DIMENSION(96) :: ELEM
REAL, DIMENSION(96) :: RAIOCOV
CHARACTER(LEN=100) :: ARQUIVO


ELEM=(/"H ","He",&                                                                      
"Li","Be","B ","C ","N ","O ","F ","Ne",&                                        
"Na","Mg","Al","Si","P ","S ","Cl","Ar",&                                        
"K ","Ca","Sc","Ti","V ","Cr","Mn","Fe","Co","Ni","Cu","Zn",&                    
"Ga","Ge","As","Se", "Br","Kr",&                                                 
"Rb","Sr","Y ","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd",&                    
"In","Sn","Sb","Te","I ","Xe",&                                                  
"Cs","Ba",&                                                                      
"La","Ce","Pr","Nd","Pm","Sm","Eu",&                                             
"Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu",&                                        
"Hf","Ta","W ","Re","Os","Ir","Pt","Au","Hg",&                                    
"Tl","Pb","Bi","Po","At","Rn",&                                                  
"Fr","Ra","Ac","Th","Pa","U ","Np","Pu","Am",&                                   
"Cm"/)
RAIOCOV=(/0.310,0.280,&                                                                    
1.280,0.960,0.840,0.760,0.710,0.660,0.570,0.580,&                                
1.660,1.410,1.210,1.110,1.070,1.050,1.020,1.060,&                                
2.030,1.760,1.700,1.600,1.530,1.390,1.500,1.420,1.380,1.240,1.320,1.220,&        
1.220,1.200,1.190,1.200,1.200,1.160,&                                            
2.200,1.950,1.900,1.750,1.640,1.540,1.470,1.460,1.420,1.390,1.450,1.440,&        
1.420,1.390,1.390,1.380,1.390,1.400,&                                            
2.440,2.150,&                                                                    
2.070,2.040,2.030,2.010,1.990,1.980,1.980,&                                      
1.960,1.940,1.920,1.920,1.890,1.900,1.870,1.870,&                                
1.750,1.700,1.620,1.510,1.440,1.410,1.360,1.360,1.320,&                          
1.450,1.460,1.480,1.400,1.500,1.500,&                                            
2.600,2.210,2.150,2.060,2.000,1.960,1.900,1.870,1.800,&                          
1.690/)

READ *, ARQUIVO

!Abre o arquivo .xyz passado como argumento
OPEN(1, FILE=ARQUIVO)
! Armazena o numero de átomos na variável NUM_ATOMS
READ(1,*) NUM_ATOMOS
! Lê a linha em branco após o numero de átomos
READ(1,*)

! Laço para armazenar os simbolos e as coordenadas de cada linha nas listas SIMBOLOS, X, Y e Z
DO I=1,NUM_ATOMOS
        READ(1,*) SIMBOLOS(I), X(I), Y(I), Z(I)
END DO

! Laços for aninhados para calcular todas distâncias possíveis
DO I=1,NUM_ATOMOS
        DO J=1,NUM_ATOMOS

! condicional que evita printar a distancia do atomo 1 para o atomo2 2 e do atomo 2 para o atomo 1 (computa as distânciasuma vez só)
               IF (I.LT.J) THEN
                        !armazena as coordenadas do átomo 1
                        ATOMO1=(/X(I),Y(I),Z(I)/)
                        !armazena as coordenadas do átomo 2
                        ATOMO2=(/X(J),Y(J),Z(J)/)
                        !calcula da distância de ligação pela formula de distância entre dois pontos
                        DISTANCIA=((ATOMO1(1)-ATOMO2(1))**2 + (ATOMO1(2)-ATOMO2(2))**2 + (ATOMO1(3)-ATOMO2(3))**2)**0.5
                        !Laço para calcular a soma dos raios de van der waals
                        DO A=1,96
                                IF (ELEM(A).EQ.SIMBOLOS(I)) THEN
                                        RAIO1=RAIOCOV(A)
                                END IF
                                IF (ELEM(A).EQ.SIMBOLOS(J)) THEN
                                        RAIO2=RAIOCOV(A)
                                END IF
                        END DO
                        VDW_SOMA=RAIO1+RAIO2
                        PARAMETRO=VDW_SOMA*1.3
                        !condicional que evita printar distancias maiores que 130%  da soma dos raios de van der waals
                        IF (DISTANCIA>0 .AND. DISTANCIA < PARAMETRO) THEN
                                !printa os resultados no terminal
                                PRINT *, SIMBOLOS(I), ' to ', SIMBOLOS(J), ':', DISTANCIA 
                
                        END IF
                END IF
        END DO
END DO     

END PROGRAM geom



