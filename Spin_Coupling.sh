#!/bin/bash
# Spin_Coupling.sh
#
# *********************************************************************************************
# *                          Spin-spin Coupling constants Extraction                          *
# *********************************************************************************************
#
# https://github.com/danielarley/Computational_Chemistry_Scripts
# 
# Esse script pode ser usado para extrair os valores da constantes de acoplamento spin-spin
# e das suas 4 contribuições aditivas geradas em um calculo de RMN do software de química
# computacional Gaussian16 (ou outra versão anterior do programa). Para mais informações 
# referentes a esse tipo de calculo acessar a página oficial do gaussian: "https://gaussian.com/".
#
# Para utilizar esse script é preciso passar no mínimo 4 argumentos. O primeiro argumento 
# sempre será o nome do arquivo que contém os dados, o segundo argumento é o tipo de constante ou
# contribuição desejada, podendo ser escolhida os tipos: KFC, JFC, KSD, JSD, KDSO, JDSO, KPSO, JDSO, 
# KTOTAL, JTOTAL. O terceiro e o quarto argumento (assim como os argumentos posteriores, caso seja 
# passado mais de 4 argumentos na linha de comando) devem ser números inteiros menores que a quantidade
# total de átomos do sistema. O script sempre irá retornar os valores das constantes solicitadas  
# entre o átomo passado no terceiro argumento com os átomos passados nos argumentos posteriores,
# assim como também irá retornar a média desses valores.
#
# Exemplo de execução do script:
#
# ./Spin_Coupling.sh arquivo.log JDSO 14 19 22 1 7
#
# Caso a linha acima seja executada o script irá acessar o arquivo "arquivo_log" e irá extrair 
# os valores da contiruição JDSO (Diamagnético Spin-Órbita) para os pares 14 e 19, 14 e 22, 14 
# e 1 e para 14 e 7. 



### Ínicio do script ###

# Analisa se a quantidade de argumentos passada é suficiente para o funcionamento do script
if [ $# -lt 4 ]; then
	echo "Digite pelo menos 4 argumentos para o script funcionar."
	echo "
	Para utilizar esse script é preciso passar no mínimo 4 argumentos. O primeiro argumento 
sempre será o nome do arquivo que contém os dados, o segundo argumento é o tipo de constante ou
contribuição desejada, podendo ser escolhida os tipos: KFC, JFC, KSD, JSD, KDSO, JDSO, KPSO, JDSO, 
KTOTAL, JTOTAL. O terceiro e o quarto argumento (assim como os argumentos posteriores, caso seja 
passado mais de 4 argumentos na linha de comando) devem ser números inteiros menores que a quantidade
total de átomos do sistema. O script sempre irá retornar os valores das constantes solicitadas  
entre o átomo passado no terceiro argumento com os átomos passados nos argumentos posteriores,
assim como também irá retornar a média desses valores.

 Exemplo de execução do script:

 ./Spin_Coupling.sh arquivo.log JDSO 14 19 22 1 7"
	exit
fi

# Filtrando a matriz de dados de acordo a constante escolhida

if [ "$2" = "KFC" ]; then
        sed -n '/Fermi Contact (FC) contribution to K (Hz):/,/Fermi Contact (FC) contribution to J (Hz):/p' $1 | sed -n '$!p' > KFC.tmp 
# O primeiro sed filtra os resultados apenas da constante desejada e o segundo sed remove a ultima linha que conrresponde ao nome da próxima 
# constante, no final o resultado do comando é escrito no arquivo KFC.tmp

elif [ "$2" = "JFC" ]; then
        sed -n '/Fermi Contact (FC) contribution to J (Hz):/,/Spin-dipolar (SD) contribution to K (Hz):/p' $1 | sed -n '$!p' > JFC.tmp

elif [ "$2" = "KSD" ]; then
        sed -n '/Spin-dipolar (SD) contribution to K (Hz):/,/Spin-dipolar (SD) contribution to J (Hz):/p' $1 | sed -n '$!p' > KSD.tmp

elif [ "$2" = "JSD" ]; then
        sed -n '/Spin-dipolar (SD) contribution to J (Hz):/,/Paramagnetic spin-orbit (PSO) contribution to K (Hz):/p' $1 | sed -n '$!p' > JSD.tmp

elif [ "$2" = "KPSO" ]; then
        sed -n '/Paramagnetic spin-orbit (PSO) contribution to K (Hz):/,/Paramagnetic spin-orbit (PSO) contribution to J (Hz):/p' $1 | sed -n '$!p' > KPSO.tmp

elif [ "$2" = "JPSO" ]; then
        sed -n '/Paramagnetic spin-orbit (PSO) contribution to J (Hz):/,/Diamagnetic spin-orbit (DSO) contribution to K (Hz):/p' $1 | sed -n '$!p' > JPSO.tmp

elif [ "$2" = "KDSO" ]; then
        sed -n '/Diamagnetic spin-orbit (DSO) contribution to K (Hz):/,/Diamagnetic spin-orbit (DSO) contribution to J (Hz):/p' $1 | sed -n '$!p' > KDSO.tmp

elif [ "$2" = "JDSO" ]; then
        sed -n '/Diamagnetic spin-orbit (DSO) contribution to J (Hz):/,/Total nuclear spin-spin coupling K (Hz):/p' $1 | sed -n '$!p' > JDSO.tmp

elif [ "$2" = "KTOTAL" ]; then
        sed -n '/Total nuclear spin-spin coupling K (Hz):/,/Total nuclear spin-spin coupling J (Hz):/p' $1 | sed -n '$!p' > KTOTAL.tmp

elif [ "$2" = "JTOTAL" ]; then
        sed -n '/Total nuclear spin-spin coupling J (Hz):/,/End of Minotr Frequency-dependent properties file/p' $1 | sed -n '$!p' > JTOTAL.tmp

else
	echo "Tipo de constante não permitida, tente novamente com algum dos tipos (KFC, JFC, KSD, JSD, KDSO, JDSO, KPSO, JDSO, KTOTAL, JTOTAL)"
	exit
fi


# Obtendo o número de átomos do sistema

# O primeiro sed filtra o padrão passado (o qual contém as coordenadas de cada átomo) e redireciona para o segundo sed que evita que o padrão passado 
# (que ocorre várias vezes ao longo do arquivo) sejá lido mais de uma vez, a saida dos comandos sed são redirecionadas para o comando wc -l que conta 
# a quantidade de linhas, por fim subtraimos 7 que são as linhas que não correspondem aos átomos e obtemos o numero de átomos desejado.
num_atoms=$[`sed -ne '/Standard/,/Rotational/p' -e '/Rotational/q' $1 | wc -l` -7]
a=$num_atoms # renomeando a variável para facilitar as expressões seguintes


# Criando lista apenas com os parametros passados da posição 4 adiante
parametros=()
for i in $@ # itera sobre os parametros passados
do
        parametros+=($i) #adiciona cada um dos parametros a lista  "parametros"
done

unset parametros[0] #Remove o primeiro parâmetro da lista
unset parametros[1] #Remove o segundo parâmetro da lista
unset parametros[2] #Remove o terceiro parâmetro da lista



# Análise se algum dos parametros digitados é maior que o número de átomos
for i in ${parametros[*]}
do
        if [ $3 -gt $num_atoms -o $i -gt $num_atoms ]; then
                echo "O sistema contém $num_atoms átomos, você digitou um número maior que esse.Tente novamente com números entre 1 e $num_atoms"
                rm *.tmp
		exit
        fi
done



# Laço for para obter os valores das constantes para todos os pares solicitados
dados=() # Lista que serão armazenadas os dados para futuro cálculo de média
for i in ${parametros[*]}
do	
# Condicional if para analisar qual dos números do par é menor e fazer a seleção da matriz de dados com base nisso
	if [ $3 -lt $i ]; then
		atomo1=$3
		atomo2=$i
	else
		atomo1=$i
		atomo2=$3
	fi

	x=$[$atomo1%5] # Calcula o resto da divisão entre o numero do átomo passado e o número 5

# Testa se o resto da divisão é igual ou diferente de zero e define em qual matriz, qual coluna
# e qual a linha que os dados serão extraídos
	if [ $x -ne 0 ]; then
		num_matriz=$[($atomo1 / 5) + 1]
		num_coluna=$[$x+2]
		num_linha=$[$atomo2-($num_matriz-1)*5]

	else
		num_matriz=$(($atomo1 / 5))
		num_coluna=7
		num_linha=$[$atomo2-($num_matriz-1)*5]
	fi

	c=$num_coluna # renomeia a variável para facilitar a escrita das expressões seguintes
	m=$num_matriz # renomeia a variável para facilitar a escrita das expressões seguintes
	l=$num_linha  # renomeia a variável para facilitar a escrita das expressões seguintes

# Laço for para obter valor necessário para o cálculo do número de linhas que se deve remover de cada matriz para filtragem dos resultados 
	soma=0
	for ((j=2;j<=$m;j++))
	{
        	soma=$[$soma+$[($m-$j)*5]]
	}

# Calculo do número de linhas que se deve remover e printar da matriz gerada no primeiro condicional para obter a apenas a matriz com os resultados das
# constantes para os pares solicitados
	num_de_linhas_para_remover=$[$m+1+($m-1)*$a-$soma]
	nlr=$num_de_linhas_para_remover # renomeia a variável para facilitar a escrita das expressões seguintes

	num_linhas_para_printar=$[$a-($m-1)*5]
	nlp=$num_linhas_para_printar    # renomeia a variável para facilitar a escrita das expressões seguintes

# Conjunto de seds para filtrar apenas o valor do par desejado
	sed -n "1,$nlr!p" $2.tmp | sed "1,$nlp!d" > $2+2.tmp # Remove e printa um numero preciso de linhas e armazena o resultado contendo apenas a matriz com
                                                             # os dados para o par solicitado em outro arquivo temporário.
	sed -i 's/     / /g' $2+2.tmp                        # Remove excesso de espaços
	sed -i 's/  / /g' $2+2.tmp                           # Remove execesso de espaços
	cut -f$c -d " " $2+2.tmp > $2+3.tmp                  # O comando cut filtra apenas a coluna com os dados para o átomo1 e escreve um terceiro arquivo temporário
	sed -i 's/D+00/*1/g' $2+3.tmp                                 # Sustitui as letras D+00 por *1
	sed -i 's/D+01/*10/g' $2+3.tmp                                # Sustitui as letras D+01 por *10
	sed -i 's/D+02/*100/g' $2+3.tmp                               # Sustitui as letras D+02 por *100
	sed -i 's/D+03/*1000/g' $2+3.tmp                              # Sustitui as letras D+03 por *1000
	sed -i 's/D+04/*10000/g' $2+3.tmp                             # Sustitui as letras D+04 por *10000
	sed -i 's/D-01/*0.1/g' $2+3.tmp                               # Sustitui as letras D-01 por *0.1
	sed -i 's/D-02/*0.01/g' $2+3.tmp                              # Sustitui as letras D-02 por *0.01
	sed -i 's/D-03/*0.001/g' $2+3.tmp                             # Sustitui as letras D-03 por *0.001
	sed -i 's/D-04/*0.0001/g' $2+3.tmp                            # Sustitui as letras D-04 por *0.0001
	sed -i 's/D-05/*0.00001/g' $2+3.tmp                           # Sustitui as letras D-05 por *0.00001
	resultado=`sed "$num_linha!d" $2+3.tmp`                       # Armazena a constante solicitada na variável resultado
	resultado_final=`bc <<EOF                                     # Transforma o resultado que está em notação cientifica em um numero decimal
scale=6
$resultado
EOF
`
	dados+=($resultado_final)                                           # Armazena todos resultados na lista dados para posterior cálculo de média
	echo "A constante $2 para os átomos $atomo1 e $atomo2 é $resultado_final"   #Printa na tela o resultado para cada par
done


# Loop for para cálculo de média
j=0 # j começa em 0 e será adicionado em cada passagem do loop um dos resultados obtidos nos passos anteriores
contador=0 # Contador começa em 0 e será adicionado uma unidade em cada passagem do loop
for i in ${dados[*]}
do
	j=`echo "scale=6;$j+$i" | bc`
	contador=$[$contador+1]
done

echo "A média entre as constantes é `echo "scale=6; $j/$contador" | bc`"

rm *.tmp  # Apaga todos arquivos temporários gerados

