# ProjetoJuno
SImulação da missão Juno da NASA. Testando as condições para a 'Deep Space Maneuver'.

Simulation of the NASA mission Juno, testing the conditions for the Deep Space Maneuver.

##

#### Introdução:

A sonda Juno é uma sonda espacial enviada pela NASA para coletar informações sobre o Planeta Júpiter. Foi lançada no dia 5 de agosto de 2011 e chegou a Júpiter no dia 5 de Julho de 2016. Um dos fatos mais impressionantes, e mais importantes para nós, nessa missão, foi seu plano de voo extremamente complexo. Esse plano de voo consistia em, após ser lançada, realizar uma trajetória de dois anos orbitando o sol e, após essa orbita, passar perto da Terra para ganhar mais velocidade com sua atração gravitacional e com isso fazer uma segunda órbita bem maior, que se encontra com Júpiter depois de mais três anos de viagem.

Esse ”estilingue” feito com a Terra é extremamente complexo por envolver diversas variáveis além da escolha adequada de parâmetros iniciais de lançamento, que afetam toda a trajetória. “This is the hardest thing NASA has ever done, that’s my claim.” [Isso é a coisa mais difı́cil que a NASA já fez, essa é minha afirmação] Scott Bolton, lider do grupo de pesquisa na missão Juno
O que almejamos neste projeto é simular essa manobra espacial e entender o quão otimizada ela é em relação a um lançamento direto, sem estilingue gravitacional. 

![Trajetoria_Juno](https://github.com/pemariano/ProjetoJuno/assets/85647121/9c06f6c3-1082-4ab5-ac12-95262b234289)

Trajetória prevista da sonda Juno.

##

#### Manual do usuário:

Primeiro se deve usar o script ProjetoJuno_Calculo.m para achar condições iniciais adequadas, e depois o script ProjetoJuno.m para ver o resultado de tais condições iniciais. 

O script ProjetoJunoCalculo tem como função otimizar a órbita em função de duas condições iniciais, o módulo da velocidade de lançamento e o ângulo de lançamento da sonda em relação ao eixo x. Podemos variar ambas as condições simultaneamente ou apenas uma. A otimização da órbita é feita tentando otimizar dois parâmetros, a distância mı́nima de Juno com a Terra no sobrevoo, que queremos minimizar, e a distância máxima entre Juno e o Sol após o sobrevoo, que queremos maximizar.

![image](https://github.com/pemariano/ProjetoJuno/assets/85647121/c84994a2-7ea1-41df-a501-2bf7e916c47e)

Plots do programa ProjetoJuno_Calculo.m com cálculo das distâncias de acordo com variadas condições iniciais. Aqui podemos ver que um determinado par de valores alcança uma distância muito maior do que o resto!

##

#### Resultado:

Por fim após inúmeras simulações conseguimos identificar condições iniciais precisas que realizam com sucesso a manobra gravitacional e a sonda Juno chega até a órbita de Júpiter! 

Uma órbita direta, sem manobra gravitacional, precisaria de uma velocidade inicial da sonda Juno de 30.000 km/h para chegar a órbita de Júpiter. A partir da simulação com a manobra conseguimos chegar na órbita de Júpiter com uma velocidade inicial da sonda Juno igual a 18450.40 km/h. Ou seja, conseguimos diminuir em um terço a velocidade inicial necessária da sonda Juno para que ela alcance Júpiter. Isso significa um foguete mais leve, mais barato e mais fácil de fazer.

![estilingue](https://github.com/pemariano/ProjetoJuno/assets/85647121/ba882a13-7b5e-40b8-873c-557b7f693c1a)

Órbita do planeta Terra (azul) e de Júpiter (marrom) e da sonda Juno (laranja). 

##

#### Arquivos:

- **PedroMariano_Projeto1_Relatorio.pdf** contém a descrição detalhada do projeto contando com manual do usuário.

- **ProjetoJuno.m** principal para cálculo dos planetas e da sonda. Esse cálculo é feito através do Método de Runge-Kutta de ordem 2, que consiste em atualizar as posições e forças de forma discreta a cada intervalo *dt*. Produz o gráfico e a animação finais. 

- **ProjetoJuno_Calculo.m** secundário para cálculo de valores ideais de lançamento da sonda Juno para se chegar a Júpiter. 
  - Cálculo da menor distancia entre Juno e a Terra depois da primeira volta em funcao do angulo de lançamento e da velocidade inicial. Buscando que a sonda Juno tenha a menor distância possível da Terra para que a manobra seja mais efetiva.
  - Cálculo da maior distancia de Juno em relacao ao Sol de 2.5 anos do lançamento em diante em funcao do angulo de lançamento e da velocidade inicial. Buscando que a sonda Juno vá o mais longe possível depois da manobra com a Terra.
