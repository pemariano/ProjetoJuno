#**********************************************************
#                                                         *
#                  PROGRAMA PRINCIPAL                     *
#                                                         *
#**********************************************************

#Programa para calculo da trajetória da sonda Juno.

#***********PARTE QUE PODE SER ALTERADA****************************************************

#Parâmetros iniciais do problema:
G = 6.674184 * 10^-11 * (149597870700)^-3 * (1.989*10^30)^1 * (8760*3600)^2; #ua^3 * Ms^-1 * ano^-2
dt = 1/365;              #anos
tfinal = 4.05;           #anos
t = 0:dt:tfinal;       
Msol = 1;                #massas solares71500
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

#Condições iniciais do nosso problema:
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
Juno(1,1) = 1 - 0.01;
#y inicial de Juno:
Juno(1,2) = 0.01;

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
Mod_Vel_Juno = 1.080400;                                #Modulo da velocidade inicial da Juno (ua/ano)
angulo = 1.540;                                         #angulo de lancamento em relacao ao eixo x
Vjuno(1,1) = Mod_Vel_Juno * cos(angulo) + Vterra(1,1);  #Vx inicial de Juno
Vjuno(1,2) = Mod_Vel_Juno * sin(angulo) + Vterra(1,2);  #Vy inicial de Juno

#***********FIM****************************************************************************

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



#Método de Runge-Kutta de ordem 2

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
  
  
  
  #intermediário Terra:
  Vterra(i,3)=Vterra(i,1) + k1_Vt_x/2;
  Terra(i,3)=Terra(i,1) + k1_Terra_x/2;
  Vterra(i,4)=Vterra(i,2) + k1_Vt_y/2;
  Terra(i,4)=Terra(i,2) + k1_Terra_y/2;
  #intermediário Sol:
  Vsol(i,3)=Vsol(i,1) + k1_Vs_x/2;
  Sol(i,3)=Sol(i,1) + k1_Sol_x/2;
  Vsol(i,4)=Vsol(i,2) + k1_Vs_y/2;
  Sol(i,4)=Sol(i,2) + k1_Sol_y/2;
  #intermediário Jupiter:
  Vjupiter(i,3)=Vjupiter(i,1) + k1_Vjp_x/2;
  Jupiter(i,3)=Jupiter(i,1) + k1_Jupiter_x/2;
  Vjupiter(i,4)=Vjupiter(i,2) + k1_Vjp_y/2;
  Jupiter(i,4)=Jupiter(i,2) + k1_Jupiter_y/2;
  #intermediário Juno:
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
 
  
  
  #Atualização Terra:
  Vterra(i+1,1)=Vterra(i,1) + k2_Vt_x;
  Terra(i+1,1)=Terra(i,1) + k2_Terra_x;
  Vterra(i+1,2)=Vterra(i,2) + k2_Vt_y;
  Terra(i+1,2)=Terra(i,2) + k2_Terra_y;
  #Atualização Sol:
  VSol(i+1,1)=Vsol(i,1) + k2_Vs_x;
  Sol(i+1,1)=Sol(i,1) + k2_Sol_x;
  Vsol(i+1,2)=Vsol(i,2) + k2_Vs_y;
  Sol(i+1,2)=Sol(i,2) + k2_Sol_y;
  #Atualização Jupiter:
  Vjupiter(i+1,1)=Vjupiter(i,1) + k2_Vjp_x;
  Jupiter(i+1,1)=Jupiter(i,1) + k2_Jupiter_x;
  Vjupiter(i+1,2)=Vjupiter(i,2) + k2_Vjp_y;
  Jupiter(i+1,2)=Jupiter(i,2) + k2_Jupiter_y;
  #Atualização Juno:
  Vjuno(i+1,1)=Vjuno(i,1) + k2_Vjn_x;
  Juno(i+1,1)=Juno(i,1) + k2_Juno_x;
  Vjuno(i+1,2)=Vjuno(i,2) + k2_Vjn_y;
  Juno(i+1,2)=Juno(i,2) + k2_Juno_y;

endfor

#Grafico da orbita
figure(1)
plot(Terra(:,1),Terra(:,2),Jupiter(:,1),Jupiter(:,2),Juno(:,1),Juno(:,2));
hold on;
grid on;
plot(Sol(1,1),Sol(1,2),'ko','LineWidth',2,'MarkerEdgeColor','k',...
      'MarkerfaceColor','r','MarkerSize',30);
str = strcat('Orbita da Terra, de Jupiter e da sonda Juno');  # gera uma string
title(str,'FontSize',18);               # coloca um título
xlabel('x(ua)','FontSize',18);          # coloca uma legenda no eixo x
ylabel('y(ua)','FontSize',18);          # coloca uma legenda no eixo y
xlim([-6 6]);
ylim([-6 6]);


#Fazer a animacao:
#Animacao não esta em escala.
figure1 = figure(2);
set(figure1,'color','white');
winsize = get(figure1,'Position');
winsize(1:2) = [0 0];
clear M;
count=1;
p1=plot(Terra(1,1),Terra(1,2),'ko','LineWidth',2,'MarkerEdgeColor','k',...
      'MarkerfaceColor','b','MarkerSize',8,...
      Sol(1,1),Sol(1,2),'ko','LineWidth',2,'MarkerEdgeColor','k',...
      'MarkerfaceColor','r','MarkerSize',30,...
      Jupiter(1,1),Jupiter(1,2),'ko','LineWidth',2,'MarkerEdgeColor','k',...
      'MarkerfaceColor',[1.0 0.7 0.0],'MarkerSize',18,...
      Juno(1,1),Juno(1,2),'ko','LineWidth',2,'MarkerEdgeColor','k',...
      'MarkerfaceColor','k','MarkerSize',1);
str = strcat('Orbita da Terra, de Jupiter e da sonda Juno');        # gera uma string
title(str,'FontSize',18);               # coloca um título
xlabel('x(ua)','FontSize',18);          # coloca uma legenda no eixo x
ylabel('y(ua)','FontSize',18);          # coloca uma legenda no eixo y
legenda = strcat('tempo = ',num2str(0),'anos',num2str(0),'dias');   # gera uma string
text(1,5.5,legenda,'FontSize',16);      # coloca uma legenda
ylim([-6 6]);                           # mantem a escala do gráfico constante
xlim([-6 6]);                           # mantem a escala do gráfico constante
grid on;                                # coloca o grid
pause(10^-6);                           # tempo entre cada imagem
M(count)=getframe(figure1);             # armazena a frame em 'M' - Octave
delete(p1);                             # deleta o gráfico (já armazenou em 'M')
count=count+1;                          # atualiza o contador  
for i=1:20:length(t)                    #coluna = tempo
  p1=plot(Terra(i,1),Terra(i,2),'ko','LineWidth',2,'MarkerEdgeColor','k',...
      'MarkerfaceColor','b','MarkerSize',8,...
      Sol(i,1),Sol(i,2),'ko','LineWidth',2,'MarkerEdgeColor','k',...
      'MarkerfaceColor','r','MarkerSize',30,...
      Jupiter(i,1),Jupiter(i,2),'ko','LineWidth',2,'MarkerEdgeColor','k',...
      'MarkerfaceColor',[1.0 0.7 0.0],'MarkerSize',18,...
      Juno(i,1),Juno(i,2),'ko','LineWidth',2,'MarkerEdgeColor','k',...
      'MarkerfaceColor','k','MarkerSize',2,...
      Terra(1:i,1),Terra(1:i,2),Jupiter(1:i,1),Jupiter(1:i,2),Juno(1:i,1),Juno(1:i,2));
  str = strcat('Orbita da Terra, de Jupiter e da sonda Juno');      # gera uma string
  title(str,'FontSize',18);             # coloca um título
  xlabel('x(ua)','FontSize',18);        # coloca uma legenda no eixo x
  ylabel('y(ua)','FontSize',18);        # coloca uma legenda no eixo y
  legenda = strcat('tempo = ',num2str(idivide(t(i),1,'fix')),'anos',...
      num2str(rem(t(i),1)*365),'dias');  # gera uma string
  text(1,5.5,legenda,'FontSize',16);    # coloca uma legenda
  ylim([-6 6]);                         # mantem a escala do gráfico constante
  xlim([-6 6]);                         # mantem a escala do gráfico constante
  grid on;                              # coloca o grid
  pause(10^-6);                         # tempo entre cada imagem
  M(count)=getframe(figure1);           # armazena a frame em 'M' - Octave
  delete(p1);                           # deleta o gráfico (já armazenou em 'M')
  count=count+1;                        # atualiza o contador  
endfor

