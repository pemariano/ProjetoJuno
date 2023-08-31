# ProjetoJuno
Simulation of the NASA mission Juno, testing the conditions for the Deep Space Maneuver

##

Arquivos:

*PedroMariano_Projeto1_Relatorio.pdf* contém a descrição detalhada do projeto contando com manual do usuário.

*ProjetoJuno.m* principal para cálculo dos planetas e da sonda. Esse cálculo é feito através do Método de Runge-Kutta de ordem 2, que consiste em atualizar as posições e forças de forma discreta a cada intervalo *dt*. Produz o gráfico e a animação finais. 

*ProjetoJuno_Calculo.m* secundário para cálculo de valores ideais de lançamento da sonda Juno para se chegar a Júpiter. 
- Cálculo da menor distancia entre Juno e a Terra depois da primeira volta em funcao do angulo de lançamento e da velocidade inicial. Buscando que a sonda Juno tenha a menor distância possível da Terra para que a manobra seja mais efetiva.
- Cálculo da maior distancia de Juno em relacao ao Sol de 2.5 anos do lançamento em diante em funcao do angulo de lançamento e da velocidade inicial. Buscando que a sonda Juno vá o mais longe possível depois da manobra com a Terra.
