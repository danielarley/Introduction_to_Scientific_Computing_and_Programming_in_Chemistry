#!/usr/bin/perl
# boltzman.pl
#
# *********************************************************************************************
# *                                   Boltzman Program                                        *
# *********************************************************************************************
#
# https://github.com/danielarley/Computational_Chemistry_Scripts
#
# Esse script pode ser usado para o calculo das diferenças de energia entre diferentes conforma-
# -ções de uma mesma molécula a partir dos outputs gerados pelo cálculo de otimização de geometria
# e frequência gerados pelo programa Gausssian 16 (ou versão anterior). O script analisa qual dos 
# conformeros é mais estável e calcula as variações de energia eletrônica, energia livre e entalpia 
# com base nele. Após a isso o script faz o cálculo das populações de Boltzman para cada um dos 
# conformeros.
#
# Para executá-lo basta passar como argumento os diferentes arquivos.log que o script se encarrega
# do resto. Exemplo: 
#
# ./boltzman.pl conformero1.log conformero2.log conformero3.log conformero4.log (...)
#


### Início do Script ###

### Código para os dados de energia ###

# Cria listas vazias para armazenar os valores de energia
@ele_zpe=(); #Eletrônica + ZPE
@ele_thermal=(); # ELetrônica + térmica
@ele_enthalpies=(); #Entalpia
@ele_gibbs=(); # Energia livre


#Laço for para abrir cada arquivo e extrair os dados de energia (a lista @ARGV armazena todos argumentos passados)
foreach $arq_log (@ARGV){
	chomp($arq_log);
	open(ARQ,"$arq_log") or die "Não foi possível abrir $arq_log: $!\n";
	while($line=<ARQ>){
                if ($line !~/\bSum of electronic and zero-point Energies\b/){   # Busca o padrão e retorna verdadeiro onde ele não estiver
                }
                else{
                        $a=$line;                                               # Atribui a linha que contém o padrão à variável $a 
			@c=split(' ',$a);                                       # Fatia a linha em varias partes e monta uma lista @c
			$b=$c[6];						# Atribui o 6º elemento de @c a varável b
			push (@ele_zpe,$b);                                     # Adiciona o valor de $b (que é a energia + zpe) a lista @ele_zpe
                }
		if ($line !~/\bSum of electronic and thermal Energies\b/){
                }
                else{
                        $a=$line;
                        @c=split(' ',$a);
                        $b=$c[6];
                        push (@ele_thermal,$b);
                }

		if ($line !~/\bSum of electronic and thermal Enthalpies\b/){
                }
                else{
                        $a=$line;
                        @c=split(' ',$a);
                        $b=$c[6];
                        push (@ele_enthalpies,$b);
                }

		if ($line !~/\bSum of electronic and thermal Free Energies\b/){
                }
                else{
                        $a=$line;
                        @c=split(' ',$a);
                        $b=$c[7];
                        push (@ele_gibbs,$b);
                }
        }
	close(ARQ);
}


# Cria listas vazias para armazenar os valores dos deltas das energias
@delta_zpe=();
@delta_thermal=();
@delta_enthalpies=();
@delta_gibbs=();

# Determina o menor valor da energia eletrônica + zpe para tomar como referencial
$menor_zpe=0;
foreach $i (@ele_zpe){
	if ($i<$menor_zpe){
		$menor_zpe=$i;
	}
}
# Calcula o delta_zpe e adiciona na lista @delta_zpe
foreach $i (@ele_zpe){
	$x=(($i-$menor_zpe)*627.5);
	push (@delta_zpe,$x);
}


# Determina o menor valor da energia eletrônica + térmica para tomar como referencial
$menor_thermal=0;
foreach $i (@ele_thermal){
        if ($i<$menor_thermal){
                $menor_thermal=$i;
        }
}
# Calcula o delta da energia térmica e adiciona na lista @delta_thermal
foreach $i (@ele_thermal){
        $x=(($i-$menor_thermal)*627.5);
        push (@delta_thermal,$x);
}


# Determina o menor valor da entalpia para tomar como referencial
$menor_enthalpies=0;
foreach $i (@ele_enthalpies){
        if ($i<$menor_enthalpies){
                $menor_enthalpies=$i;
        }
}
# Calcula o delta das entalpias e adiciona na lista @delta_enthalpies
foreach $i (@ele_enthalpies){
        $x=(($i-$menor_enthalpies)*627.5);
        push (@delta_enthalpies,$x);
}


# Determina o menor valor da energia livre para tomar como referencial
$menor_gibbs=0;
foreach $i (@ele_gibbs){
        if ($i<$menor_gibbs){
                $menor_gibbs=$i;
        }
}
# Calcula o delta da energia livre e adiciona na lista @delta_gibbs
foreach $i (@ele_gibbs){
        $x=(($i-$menor_gibbs)*627.5);
        push (@delta_gibbs,$x);
}



### Código para cálculo das populações ###

@zpe_boltz=();
@thermal_bolz=();
@enthalpies_boltz=();
@gibbs_boltz=();

# Calculando o valor das razões da população de cada estado pela população do estado mais provável
foreach $i (@delta_zpe){
	$nx_ne=2.71828182845**(-($i/(298.15*1.998557)));
	push (@zpe_boltz,$nx_ne);
}
# Somando todas razões calculadas no passo anterior
$soma_nx_ne=0;
foreach $i (@zpe_boltz){
	$soma_nx_ne=$soma_nx_ne+$i;
}

$ne=1/$soma_nx_ne;  #Obtém o valor da população do estado mais provável

#calculando a população de cada um dos estados restantes
@pop_zpe=();
foreach $i (@zpe_boltz){
	$nx=$i*$ne*100;
	push (@pop_zpe,$nx);
}

# Repetindo os passos anteriores usando a energia térmica
foreach $i (@delta_thermal){
	$nx_ne=2.71828182845**(-($i/(298.15*1.998557)));
	push (@thermal_boltz,$nx_ne);
}

$soma_nx_ne=0;
foreach $i (@thermal_boltz){
	$soma_nx_ne=$soma_nx_ne+$i;
}

$ne=1/$soma_nx_ne;

@pop_thermal=();
foreach $i (@thermal_boltz){
	$nx=$i*$ne*100;
	push (@pop_thermal,$nx);
}

# Repetindo os passos anteriores usando as entalpias
foreach $i (@delta_enthalpies){
	$nx_ne=2.71828182845**(-($i/(298.15*1.998557)));
	push (@enthalpies_boltz,$nx_ne);
}

$soma_nx_ne=0;
foreach $i (@enthalpies_boltz){
	$soma_nx_ne=$soma_nx_ne+$i;
}

$ne=1/$soma_nx_ne;

@pop_enthalpies=();
foreach $i (@enthalpies_boltz){
	$nx=$i*$ne*100;
	push (@pop_enthalpies,$nx);
}

# Repetindo os passos anteriores usando a energia livre
foreach $i (@delta_gibbs){
	$nx_ne=2.71828182845**(-($i/(298.15*1.998557)));
	push (@gibbs_boltz,$nx_ne);
}

$soma_nx_ne=0;
foreach $i (@gibbs_boltz){
	$soma_nx_ne=$soma_nx_ne+$i;
}

$ne=1/$soma_nx_ne;

@pop_gibbs=();
foreach $i (@gibbs_boltz){
	$nx=$i*$ne*100;
	push (@pop_gibbs,$nx);
}


### escrevendo resultados na tela ###

$len=@delta_zpe;		 #Pega o numero de elementos da lista delta_zpe (que é o mesmo numero para as outras listas)

print "***********************************************************************************************************\n*		 				 Boltzman Program		      		          *\n*						     (v 1.0)					          *\n*					           Daniel Arley	 			                  *\n***********************************************************************************************************\n";

# Escreve na tela os resultados das variações de energia dos conformetos em relação ao conformero mais estável
print "DELTA DE ENERGIA (kcal/mol)\n";
print "\t\t     D(E+ZPE)\t\t    D(E+Thermal)\t\tDH\t\t\tDG\n";
for ($i=0;$i<$len;++$i){
	print "$ARGV[$i]\t $delta_zpe[$i]\t $delta_thermal[$i]\t $delta_enthalpies[$i]\t $delta_gibbs[$i]\n";
}
print "\n";

# Escreve na tela as populações de cada conformero com base no tipo de energia usada
print "POPULAÇÃO DAS CONFORMAÇÕES\n";
print "\t\t     D(E+ZPE)\t\t    D(E+Thermal)\t\tDH\t\t\tDG\n";
for ($i=0;$i<$len;++$i){
	print "$ARGV[$i]\t $pop_zpe[$i]\t $pop_thermal[$i]\t $pop_enthalpies[$i]\t $pop_gibbs[$i]\n";
}

