#!/usr/bin/perl

# Script para cálculo de distâncias de ligação
$arq_xyz=$ARGV[0]; #Define a variável arq_xyz como sendo o primeiro argumento passado na linha de comando
chomp($arq_xyz);
open(ARQ, "<$arq_xyz");
$num_atoms=<ARQ>; #A primeira linha de um arquivo.xyz é sempre o número de átomos, armazena esse valor na variavel num_atoms
$remove=<ARQ>; #armazena a linha em branco entre o numero de átomos e as coordenadas na variavel remove (cada linha que é lida uma vez não é lida novamente)


#cria listas vazias para armazenas os valores
@symbols=();
@x=();
@y=();
@z=();

#laço while para ler as coordenadas do arquivo.xyz
while ($line=<ARQ>){
	@atom_data=split(' ',$line); #cada linha é transformada em uma lista do tipo (simbolo atomico, coordenada x, coordenada y, coordenada z)
        push(@symbols,$atom_data[0]); #adiciona os simbolos dos elementos na lista symbols
        push(@x,$atom_data[1]); #adiciona as coordenadas x na lista @x
        push(@y,$atom_data[2]); #adiciona as coordenadas y na lista @y
        push(@z,$atom_data[3]); #adiciona as coordenadas z na lista @z
}

# Laços for aninhados para calcular todas distâncias possíveis
for ($i=0;$i<$num_atoms;++$i){
	for ($j=0;$j<$num_atoms;++$j){
                if ($i<$j){		# condicional que evita printar a distancia do atomo 1 para o atomo2 2 e do atomo 2 para o atomo 1 (computa as distânciasuma vez só)
			@atom1=($x[$i],$y[$i],$z[$i]); #armazena as coordenadas do átomo 1
			@atom2=($x[$j],$y[$j],$z[$j]); #armazena as coordenadas do átomo 2
			$x_distance=$atom1[0]-$atom2[0];
			$y_distance=$atom1[1]-$atom2[1];
			$z_distance=$atom1[2]-$atom2[2];
			$bond_lenght=($x_distance**2 + $y_distance**2 + $z_distance**2)**0.5; #calcula da distância de ligação pela formula de distância entre dois pontos
			if ($bond_lenght>0 and $bond_lenght<1.5){             #condicional que evita printar distancias maiores que 1.5 (caso desejado pode trocar para um valor de corte maior)
				print "$symbols[$i] to $symbols[$j] : $bond_lenght\n"; #printa os resultados no terminal
			}
		}
	}


}

close(ARQ);

