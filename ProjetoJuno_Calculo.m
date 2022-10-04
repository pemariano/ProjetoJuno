#**********************************************************
#                                                         *
#       CÁLCULO DE DISTÂNCIA JUNO-TERRA E JUNO-SOL        *
#                                                         *
#**********************************************************

#Programa para calculo da menor distancia entre Juno e a Terra depois da primeira 
#volta em funcao do angulo de lançamento e da velocidade inicial.

#Programa para calculo da maior distancia de Juno em relacao ao Sol de 2.5 anos em 
#diante em funcao do angulo de lançamento e da velocidade inicial.

#Utilizado para saber quais as melhores condicoes iniciais de lançamento para se 
#chegar a Jupiter.

#Pode-se alterar as condicoes inciais do problema, das linhas 44 a 106 e tambem
#as linhas 308 a 310 que constituem parametros para analise dos dados.


#********************INICIO DO PROGRAMA: ***********************************************


tempo1 = time(); #para calcular o tempo que o programa demora para rodar
#tempo para 1 condicao inicial: 9.038 ~ 11.769 segundos com o passo temporal de 1 dia.

#Funcao que calcula a forca que astro1 faz em astro2 na coordenada x ou y no tempo t, 
#recebe o indice do tempo, a massa, a posição dos astros e uma string com a coordenada.
function [f] = Fg(astro1,astro2,t,coord,massa)
  if coord == "x"
    i = 1;
  elseif coord == "y"
    i = 2;
  elseif coord == "xint"
    i = 3;
  elseif coord == "yint"
    i = 4;
  endif
  f = -massa*39.434/((astro1(t,i)-astro2(t,i))^2+(astro1(t,(i+(-1)^(i+1)))-...
      astro2(t,(i+(-1)^(i+1))))^2)^(3/2) * (astro2(t,i)-astro1(t,i));
endfunction



#***********PARTE QUE PODE SER ALTERADA****************************************************

#Parametros iniciais do problema:
G = 6.674184 * 10^-11 * (149597870700)^-3 * (1.989*10^30)^1 * (8760*3600)^2; #ua^3 * Ms^-1 * ano^-2
dt = 1/365;          #anos
tfinal = 6;              #anos
t = 0:dt:tfinal;         #vetor dos tempos
Msol = 1;                #massas solares
Mterra = 3.003*10^-6;    #massas solares
Mjupiter = 9.546*10^-4;  #massas solares

#Posição a cada t: Px=(:,1), Py=(:,2), Pxint=(:,3), Pyint=(:,4).
Sol = zeros(length(t),4);
Terra = zeros(length(t),4);
Jupiter = zeros(length(t),4);
Juno = zeros(length(t),4);
#Velocidade a cada t: Vx=(:,1), Vy=(:,2), Vxint=(:,3), Vyint=(:,4).
Vsol = zeros(length(t),4);
Vterra = zeros(length(t),4);
Vjupiter = zeros(length(t),4);
Vjuno = zeros(length(t),4);

#Condicoes iniciais do nosso problema:
#Posições (ua)
#x inicial da Terra:
Terra(1,1) = 1;
#y inicial da Terra:
Terra(1,2) = 0;
#x inicial do Sol:
Sol(1,1) = 0;
#y inicial do Sol:
Sol(1,2) = 0;
#x inicial de Jupiter:
Jupiter(1,1) = 1.1021;
#y inicial de Jupiter:
Jupiter(1,2) = 5.0938;
#x inicial de Juno:
Juno(1,1) = 1 + 0*4*10^-5;
#y inicial de Juno:
Juno(1,2) = 4*10^-5;

#Velocidades (ua/ano) 
#Vx inicial da Terra:
Vterra(1,1) = 0;
#Vy inicial da Terra:
Vterra(1,2) = 109040 * 8760/149597871; 
#Vx inicial do Sol:
Vsol(1,1) = 0;
#Vy inicial do Sol:
Vsol(1,2) = 0;
#Vx inicial de Jupiter:
Vjupiter(1,1) = -2.6883;
#Vy inicial de Jupiter:
Vjupiter(1,2) = 0.58625;

#Condicoes iniciais da sonda Juno:
Mod_Vel_Juno = 1.0;               #Módulo da velocidade inicial da Juno (ua/ano)
angulo = 1.5;                     #angulação em relação ao eixo x do lançamento
passo_angulo = 0.01;              #passo de aumento do angulo a cada simulacao (rad)
num_simulacoes_ang = 10;          #numero de simulacoes do angulo a serem feitas 
passo_velocidade = 0.01;          #passo de aumento da velocidade a cada simulacao (ua)
num_simulacoes_vel = 10;          #numero de simulacoes da velocidade a serem feitas 

#***********FIM****************************************************************************

#indices do intervalo de tempo no qual Juno acaba a primeira volta:
cont = 1;
while t(cont) < 2
  cont = cont + 1;
endwhile
t_inicial = idivide(cont - cont/5,1);
t_final = idivide(cont + cont/5,1);

#Duracao aproximada do programa:
#Caso veja que va demorar muito tempo pode parar o programa e escolher um
#numero menor de passos. Lembrando que o tempo para cada cond. inicial é de
#aproximadamente 10s com o passo temporal de 1 dia.
segundos = num_simulacoes_ang * num_simulacoes_vel * 10/(365*dt);
disp ('Tempo aproximado para calcular todas as condicoes iniciais:')
disp('(hh:mm:ss)')
disp(datestr(segundos/(24*60*60), 'HH:MM:SS'))

#Matrizes que receberao a distancia para cada config. inicial:
aprox_JT = zeros(num_simulacoes_ang,num_simulacoes_vel);
dist_JS = zeros(num_simulacoes_ang,num_simulacoes_vel);
#Vetores dos valores de angulo e velocidade:
angulos = angulo:passo_angulo:angulo + passo_angulo * (num_simulacoes_ang-1);
velocidades = Mod_Vel_Juno:passo_velocidade:Mod_Vel_Juno + passo_velocidade * (num_simulacoes_vel-1);

Mod_Vel_Juno = Mod_Vel_Juno - passo_velocidade;
angulo = angulo - passo_angulo;

for vel=1:num_simulacoes_vel
  Mod_Vel_Juno = Mod_Vel_Juno + passo_velocidade;
  
    for ang=1:num_simulacoes_ang
      angulo = angulo + passo_angulo;
      #Vx inicial de Juno:
      Vjuno(1,1) = Mod_Vel_Juno * cos(angulo) + Vterra(1,1); #velocidade de lançamento x mais a da Terra em x
      #Vy inicial de Juno:
      Vjuno(1,2) = Mod_Vel_Juno * sin(angulo) + Vterra(1,2); #velocidade de lançamento y mais a da Terra em y
   
      
      for i=1:length(t)-1;
        
        #k1 Terra:
        k1_Vt_x= Fg(Sol,Terra,i,"x",Msol) * dt + Fg(Jupiter,Terra,i,"x",Mjupiter) * dt;
        k1_Terra_x=Vterra(i,1)*dt; 
        k1_Vt_y= Fg(Sol,Terra,i,"y",Msol) * dt + Fg(Jupiter,Terra,i,"y",Mjupiter) * dt;
        k1_Terra_y=Vterra(i,2)*dt; 
        #k1 Sol:
        k1_Vs_x= -Fg(Sol,Terra,i,"x",Mterra) * dt + Fg(Jupiter,Sol,i,"x",Mjupiter) * dt;
        k1_Sol_x=Vsol(i,1)*dt; 
        k1_Vs_y= -Fg(Sol,Terra,i,"y",Mterra) * dt + Fg(Jupiter,Sol,i,"y",Mjupiter) * dt;
        k1_Sol_y=Vsol(i,2)*dt; 
        #k1 Jupiter:
        k1_Vjp_x= Fg(Sol,Jupiter,i,"x",Msol) * dt + Fg(Terra,Jupiter,i,"x",Mterra) * dt;
        k1_Jupiter_x=Vjupiter(i,1)*dt; 
        k1_Vjp_y= Fg(Sol,Jupiter,i,"y",Msol) * dt + Fg(Terra,Jupiter,i,"y",Mterra) * dt;
        k1_Jupiter_y=Vjupiter(i,2)*dt; 
        #k1 Juno:
        k1_Vjn_x= Fg(Sol,Juno,i,"x",Msol) * dt + Fg(Terra,Juno,i,"x",Mterra) * dt + Fg(Jupiter,Juno,i,"x",Mjupiter) * dt;
        k1_Juno_x=Vjuno(i,1)*dt; 
        k1_Vjn_y= Fg(Sol,Juno,i,"y",Msol) * dt + Fg(Terra,Juno,i,"y",Mterra) * dt + Fg(Jupiter,Juno,i,"y",Mjupiter) * dt;
        k1_Juno_y=Vjuno(i,2)*dt; 
        
        
        
        #intermediario Terra:
        Vterra(i,3)=Vterra(i,1) + k1_Vt_x/2;
        Terra(i,3)=Terra(i,1) + k1_Terra_x/2;
        Vterra(i,4)=Vterra(i,2) + k1_Vt_y/2;
        Terra(i,4)=Terra(i,2) + k1_Terra_y/2;
        #intermediario Sol:
        Vsol(i,3)=Vsol(i,1) + k1_Vs_x/2;
        Sol(i,3)=Sol(i,1) + k1_Sol_x/2;
        Vsol(i,4)=Vsol(i,2) + k1_Vs_y/2;
        Sol(i,4)=Sol(i,2) + k1_Sol_y/2;
        #intermediario Jupiter:
        Vjupiter(i,3)=Vjupiter(i,1) + k1_Vjp_x/2;
        Jupiter(i,3)=Jupiter(i,1) + k1_Jupiter_x/2;
        Vjupiter(i,4)=Vjupiter(i,2) + k1_Vjp_y/2;
        Jupiter(i,4)=Jupiter(i,2) + k1_Jupiter_y/2;
        #intermediario Juno:
        Vjuno(i,3)=Vjuno(i,1) + k1_Vjn_x/2;
        Juno(i,3)=Juno(i,1) + k1_Juno_x/2;
        Vjuno(i,4)=Vjuno(i,2) + k1_Vjn_y/2;
        Juno(i,4)=Juno(i,2) + k1_Juno_y/2;
        
        
        
        #k1 Terra:
        k2_Vt_x= Fg(Sol,Terra,i,"xint",Msol) * dt + Fg(Jupiter,Terra,i,"xint",Mjupiter) * dt;
        k2_Terra_x=Vterra(i,3)*dt; 
        k2_Vt_y= Fg(Sol,Terra,i,"yint",Msol) * dt + Fg(Jupiter,Terra,i,"yint",Mjupiter) * dt;
        k2_Terra_y=Vterra(i,4)*dt; 
        #k1 Sol:
        k2_Vs_x= -Fg(Sol,Terra,i,"xint",Mterra) * dt + Fg(Jupiter,Sol,i,"xint",Mjupiter) * dt;
        k2_Sol_x=Vsol(i,3)*dt; 
        k2_Vs_y= -Fg(Sol,Terra,i,"yint",Mterra) * dt + Fg(Jupiter,Sol,i,"yint",Mjupiter) * dt;
        k2_Sol_y=Vsol(i,4)*dt; 
        #k1 Jupiter:
        k2_Vjp_x= Fg(Sol,Jupiter,i,"xint",Msol) * dt + Fg(Terra,Jupiter,i,"xint",Mterra) * dt;
        k2_Jupiter_x=Vjupiter(i,3)*dt; 
        k2_Vjp_y= Fg(Sol,Jupiter,i,"yint",Msol) * dt + Fg(Terra,Jupiter,i,"yint",Mterra) * dt;
        k2_Jupiter_y=Vjupiter(i,4)*dt; 
        #k1 Juno:
        k2_Vjn_x= Fg(Sol,Juno,i,"xint",Msol) * dt + Fg(Terra,Juno,i,"xint",Mterra) * dt + Fg(Jupiter,Juno,i,"xint",Mjupiter) * dt;
        k2_Juno_x=Vjuno(i,3)*dt; 
        k2_Vjn_y= Fg(Sol,Juno,i,"yint",Msol) * dt + Fg(Terra,Juno,i,"yint",Mterra) * dt + Fg(Jupiter,Juno,i,"yint",Mjupiter) * dt;
        k2_Juno_y=Vjuno(i,4)*dt; 
       
        
        
        #Atualizacao Terra:
        Vterra(i+1,1)=Vterra(i,1) + k2_Vt_x;
        Terra(i+1,1)=Terra(i,1) + k2_Terra_x;
        Vterra(i+1,2)=Vterra(i,2) + k2_Vt_y;
        Terra(i+1,2)=Terra(i,2) + k2_Terra_y;
        #Atualizacao Sol:
        VSol(i+1,1)=Vsol(i,1) + k2_Vs_x;
        Sol(i+1,1)=Sol(i,1) + k2_Sol_x;
        Vsol(i+1,2)=Vsol(i,2) + k2_Vs_y;
        Sol(i+1,2)=Sol(i,2) + k2_Sol_y;
        #Atualizacao Jupiter:
        Vjupiter(i+1,1)=Vjupiter(i,1) + k2_Vjp_x;
        Jupiter(i+1,1)=Jupiter(i,1) + k2_Jupiter_x;
        Vjupiter(i+1,2)=Vjupiter(i,2) + k2_Vjp_y;
        Jupiter(i+1,2)=Jupiter(i,2) + k2_Jupiter_y;
        #Atualizacao Juno:
        Vjuno(i+1,1)=Vjuno(i,1) + k2_Vjn_x;
        Juno(i+1,1)=Juno(i,1) + k2_Juno_x;
        Vjuno(i+1,2)=Vjuno(i,2) + k2_Vjn_y;
        Juno(i+1,2)=Juno(i,2) + k2_Juno_y;
      
      
    endfor
    
    #Distancia^2 dada uma velocidade e um angulo: (ang,vel)
    aprox_JT(ang,vel) = sqrt(min((Terra(t_inicial:t_final,1)-Juno(t_inicial:t_final,1)).^2 + (Juno(t_inicial:t_final,2)-Terra(t_inicial:t_final,2)).^2)); #calcula a distancia minima apos a primeira volta
    dist_JS(ang,vel) = sqrt(max((Juno(t_final:length(t),1)-Sol(t_final:length(t),1)).^2 + (Sol(t_final:length(t),2)-Juno(t_final:length(t),2)).^2)); #calcula a distância maxima apos a primeira volta
    
  endfor
endfor  

#********************OPCAO DE SALVAR OS DADOS: ******************************

#salvar 1
fid = fopen('dataAprox_JT.txt', 'w+');
for i=1:size(aprox_JT, 1)
    fprintf(fid, '%f ', aprox_JT(i,:));
    fprintf(fid, '\n');
end
fclose(fid);
csvwrite('dataAprox_JT.txt', aprox_JT);

#salvar 2
fid = fopen('dataDist_JS.txt', 'w+');
for i=1:size(dist_JS, 1)
    fprintf(fid, '%f ', dist_JS(i,:));
    fprintf(fid, '\n');
end
fclose(fid);
csvwrite('dataDist_JS.txt', dist_JS);

#********************GRAFICOS: **********************************************

if num_simulacoes_ang>1 && num_simulacoes_vel>1
  figure(1)
  subplot(1,2,1) #Plot 3D aproximacao com a Terra
  mesh(velocidades,angulos,aprox_JT);
  colorbar;
  hold on;
  title('Distancia minima entre Juno e a Terra apos a primeira volta','FontSize',16);  # coloca um titulo   
  ylabel('Angulo inicial (rad)','FontSize',15,'fontweight','bold');                    # coloca uma legenda no eixo y   
  xlabel('Velocidade inicial (ua/ano)','FontSize',15,'fontweight','bold');             # coloca uma legenda no eixo x      
  zlabel('Distancia (ua)','FontSize',15,'fontweight','bold');                          # coloca uma legenda no eixo z       
  zlim([0 5]);
  set(gca,'zdir','reverse');                          #inverte o eixo z  

  subplot(1,2,2) #Curvas de nivel aproximacao com a Terra
  contour(velocidades,angulos,aprox_JT);
  colorbar;
  hold on;
  title('Distancia minima (ua) entre Juno e a Terra apos a primeira volta','FontSize',16);  # coloca um titulo   
  ylabel('Angulo inicial (rad)','FontSize',15,'fontweight','bold');                         # coloca uma legenda no eixo y   
  xlabel('Velocidade inicial (ua/ano)','FontSize',15,'fontweight','bold');                  # coloca uma legenda no eixo x       
  print Distancia_minima_3D_apos_1_volta.pdf                                                # opção de salvar o plot   


  figure(2)
  subplot(1,2,1) #Plot 3D distancia do Sol
  mesh(velocidades,angulos,dist_JS);
  colorbar;
  hold on;
  title('Distancia maxima entre Juno e o Sol apos 2.5 anos','FontSize',16);  # coloca um titulo   
  ylabel('Angulo inicial (rad)','FontSize',15,'fontweight','bold');          # coloca uma legenda no eixo y   
  xlabel('Velocidade inicial (ua/ano)','FontSize',15,'fontweight','bold');   # coloca uma legenda no eixo x      
  zlabel('Distancia (ua)','FontSize',15,'fontweight','bold');                # coloca uma legenda no eixo z       

  subplot(1,2,2) #Curvas de nivel distancia do Sol
  contour(velocidades,angulos,dist_JS);
  colorbar;
  hold on;
  title('Distancia maxima (ua) entre Juno e o Sol apos 2.5 anos','FontSize',16);  # coloca um titulo   
  ylabel('Angulo inicial (rad)','FontSize',15,'fontweight','bold');               # coloca uma legenda no eixo y   
  xlabel('Velocidade inicial (ua/ano)','FontSize',15,'fontweight','bold');        # coloca uma legenda no eixo x     
  print Distancia_Maxima_3D_apos_2_5_anos.pdf                                     # opção de salvar o plot
else
  if  num_simulacoes_ang == 1
    figure(1)
    plot(velocidades,aprox_JT(1,:)) #plot 2D, aproximacao com a Terra x velocidade
    hold on;
    title('Distancia minima (ua) entre Juno e Terra apos a primeira volta','FontSize',18);  # coloca um titulo   
    ylabel('Distancia (ua)','FontSize',15,'fontweight','bold');                             # coloca uma legenda no eixo y   
    xlabel('Velocidade inicial (ua/ano)','FontSize',15,'fontweight','bold');                # coloca uma legenda no eixo x 
    
    figure(2)
    plot(velocidades,dist_JS(1,:)) #plot 2D, distancia do Sol x velocidade
    hold on;
    title('Distancia maxima (ua) entre Juno e o Sol apos 2.5 anos','FontSize',18);  # coloca um titulo   
    ylabel('Distancia (ua)','FontSize',15,'fontweight','bold');                     # coloca uma legenda no eixo y   
    xlabel('Velocidade inicial (ua/ano)','FontSize',15,'fontweight','bold');        # coloca uma legenda no eixo x
  elseif num_simulacoes_vel == 1
    figure(3)
    plot(angulos,aprox_JT(:,1)) #plot 2D, aproximacao com a Terra x angulo
    hold on;
    title('Distancia minima (ua) entre Juno e Terra apos a primeira volta','FontSize',18);  # coloca um titulo   
    ylabel('Distancia (ua)','FontSize',15,'fontweight','bold');                             # coloca uma legenda no eixo y   
    xlabel('Angulo inicial (rad)','FontSize',15,'fontweight','bold');                       # coloca uma legenda no eixo x
  
    figure(4)
    plot(angulos,dist_JS(:,1)) #plot 2D, distancia do Sol x angulo
    hold on;
    title('Distancia maxima (ua) entre Juno e o Sol apos 2.5 anos','FontSize',18);  # coloca um titulo   
    ylabel('Distancia (ua)','FontSize',15,'fontweight','bold');                     # coloca uma legenda no eixo y   
    xlabel('Angulo inicial (rad)','FontSize',15,'fontweight','bold');               # coloca uma legenda no eixo x
  endif

endif

#********************ANALISE: ***********************************************
if num_simulacoes_ang>1 && num_simulacoes_vel>1
  #Condicao de valores maximos e minimos nas matrizes:
  #PODEM SER ALTERADOS:
  aprox_JT_B = aprox_JT < 0.01;       #aproximacao (ua)
  dist_JS_B = dist_JS > 4.9;            #distancia (ua)
  #ACABA ALTERACAO

  A = dist_JS_B .* aprox_JT_B;          #junta as duas condicoes, matriz de 1's e 0's
  valores_possiveis = zeros(2,length(find(A==1))); #angulos=(:,1), velocidades=(:,2)
  n = 1;
  disp('')
  disp ('Numero de valores que se enquadram na solicitacao:')
  disp (length(find(A==1)))

  if length(find(A==1))>0 #se tiver um valor que se adeque as condicoes
    for i=1:length(A(:,1))
      for j=1:length(A(1,:)) 
        if A(i,j) == 1
          valores_possiveis(n,1) = angulos(i); #angulo
          valores_possiveis(n,2) = velocidades(j); #velocidade
          n = n+1;
        endif
      endfor
    endfor                                                 
  endif

  #********************MAIS GRAFICOS: ******************************************

  #grafico dos valores possiveis de angulo e velocidade inicial
  if length(find(A==1))>0 #se tiver um valor que se adeque as condicoes
    figure(5)
    plot(valores_possiveis(:,2),valores_possiveis(:,1),'ro');hold on;grid on;
    title('Valores possiveis','FontSize',20);  # coloca um titulo   
    ylabel('Angulo inicial (rad)','FontSize',16,'fontweight','bold');          # coloca uma legenda no eixo x   
    xlabel('Velocidade inicial (ua/ano)','FontSize',16,'fontweight','bold');   # coloca uma legenda no eixo y
  endif
endif
#agora podemos usar esses valores na funcao principal ProjetoJuno_.m ! :)
#tambem podemos olhar os pontos onde ficamos mais proximos da terra ou mais longe
#do sol e obter uma definicao maior desses pontos usando um passo menor
#ao redor desses pontos.

#********************MARCA O TEMPO TOTAL DO PROGRAMA: ************************

tempo2 = time(); #tempo que o programa roda
#tempo que o programa demorou:
disp('')
disp ('Tempo que o programa demorou:')
disp('(hh:mm:ss)')
disp(datestr((tempo2-tempo1)/(24*60*60), 'HH:MM:SS'))

