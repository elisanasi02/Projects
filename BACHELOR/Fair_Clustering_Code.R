rm(list=ls())
dataset=Dataset_Massachusetts
dataset2=Dataset_Massachusetts
rm(list="Dataset_Massachusetts")
#1816 scuole e 16 variabili

names(dataset)
# Installo i pacchetti necessari
library(fpc)
library(cluster)
library(rgl)
library(rattle)
library(mclust)
library(corrplot)
library(factoextra)
install.packages("ggmap")
library(ggmap)
library(ggplot2)
library(httr)
library(jsonlite)
library(devtools)
devtools::install_github("HristoInouzhe/AttractionRepulsionClustering",force=TRUE)
library(AttractionRepulsionClustering)
install.packages("geodist")
library(geodist)


# Per usare il pacchetto ggmap e poi la funzione geocode è necessario abilitare i
# servizi google creando una API key e abilitando i servizi di mapping e geocoding

# La mia API key privata
register_google(key = "***",write=TRUE)

# Per assicurarmi che le scuole vengano ricercate nello stato del Massachusetts,
# aggiungo nel nome lo stato.
nomi= paste(dataset$SCH_NAME,"Massachusetts",sep="")

# con geocode ricavo latitudine e longitudine delle scuole
coordinate_scuole=geocode(nomi,source="google")
summary(coordinate_scuole)
coordinate_scuole=as.data.frame(coordinate_scuole)

# creo una copia delle coordinate per fare delle prove
prova=coordinate_scuole
prova=as.data.frame(prova)
row.names(prova)=dataset$SCH_NAME

prova[which(prova[,1]< (-73.5)),] #queste sono le 2 scuole localizzate troppo a ovest
prova=prova[-which(prova[,1]< (-73.5)),]
prova[which(prova[,1]> (-69.9)),] #0
prova[which(prova[,2]> (43)),]     #0
prova[which(prova[,2]< (41.2)),]   #0

summary(prova)
prova = prova[!is.na(prova[,1]),]

#FINALE: 1807 scuole localizzate correttamente nel Massachusetts
dataset=dataset[dataset$SCH_NAME %in% row.names(prova),]
rm(list="coordinate_scuole")
coordinate_scuole=prova
rm(list="prova")
# Rappresentazione geografica
ggmap(get_map(c(left = -73.5, bottom = 41.2, right = -69.9, top = 43),source = "stamen"))+
  geom_point(data = coordinate_scuole, aes(x = lon, y = lat), color = "blue", size = 1,na.rm =F)

# Proporzione delle varie etnie nel dataset
p_tot_hispanic= sum(dataset$Hispanic)/sum(dataset$Total)
p_tot_NAmerica= sum(dataset$Namerican)/sum(dataset$Total)
p_tot_Asian= sum(dataset$Asian)/sum(dataset$Total)
p_tot_HP= sum(dataset$`Hawaiian/Pacific`)/sum(dataset$Total)
p_tot_black= sum(dataset$Black)/sum(dataset$Total)
p_tot_white= sum(dataset$White)/sum(dataset$Total)

p_tot= c(p_tot_hispanic,p_tot_NAmerica,p_tot_Asian,p_tot_HP,p_tot_black,p_tot_white)
sum(p_tot)

# Matrice delle distanze tra le scuole
distanze=geodist(coordinate_scuole,paired = T,measure="geodesic")
distanze=as.matrix(distanze)
row.names(distanze)=dataset$SCH_NAME
colnames(distanze)=dataset$SCH_NAME
dim(distanze) #matrice 1807 X 1807

# Controlli sulla matrice delle distanze
sum(distanze-t(distanze))
distanze[100,100]
min(distanze)

distanze=as.dist(distanze) #hclust vuole un oggetto dist come input

# Analisi dei gruppi classica
# LEGAME COMPLETO
distanze=as.dist(distanze)
completo=hclust(distanze)
plot(completo,main="Complete method",labels=F) 


# CONSIDERO LE PARTIZIONI CON 2 FINO A 15 GRUPPI E CALCOLO SI E II PER OGNUNA DI ESSE
II_comp=c() #vettore in cui raccolgo Imbalance Index
SI_comp=c() #vettore in cui raccolgo Silhouette Index
for (num_clu in 2:15){
  cluster <- cutree(completo,k=num_clu)
  SI_comp=c(SI_comp,summary(silhouette(cluster,dist = distanze))$avg.width)
  prop=e2DistProp(cluster,dataset[,4:9],p_tot)
  addendi=c()
  for(i in 1:num_clu){
    addendi=c(addendi,sqrt(sum((prop$clusterProportions[i,]-p_tot)^2)))
  }
  II_comp=c(II_comp,sum(addendi)*(1/num_clu))
  
}
II_comp 
SI_comp

#LEGAME SINGOLO
singolo= hclust(distanze,method = "single")
plot(singolo, main='Single method',labels = F) 


#LEGAME MEDIO
medio=hclust(distanze,method = "average")
plot(medio,main="average method",labels = F)

II_medio=c() #vettore in cui raccolgo Imbalance Index
SI_medio=c() #vettore in cui raccolgo Silhouette Index
for (num_clu in 2:15){
  cluster <- cutree(medio,k=num_clu)
  SI_medio=c(SI_medio,summary(silhouette(cluster,dist = distanze))$avg.width)
  prop=e2DistProp(cluster,dataset[,4:9],p_tot)
  addendi=c()
  for(i in 1:num_clu){
    addendi=c(addendi,sqrt(sum((prop$clusterProportions[i,]-p_tot)^2)))
  }
  II_medio=c(II_medio,sum(addendi)*(1/num_clu))
}
II_medio
SI_medio

#Rappresentazione degli Imbalance Index dei metodi COMPLETO E MEDIO
plot(2:15,II_comp,type = "l",col="blue",ylim = c(0,0.7),main="Imbalance Index e 
     Silhouette Index", xlab = "Numero di cluster",ylab="")
lines(2:15,II_medio,col='red')
legend("topright",legend = c("II completo","II medio","SI completo","SI medio"),
       col=c("blue","red","blue","red"),lty = c(1,1,2,2),lwd=1)
lines(2:15,SI_comp,type = "l",lty=2,col="blue")
lines(2:15,SI_medio,type = "l",lty=2,col="red")

# Scelta1: metodo del legame medio con 6 gruppi
II_medio[5] #0.12
SI_medio[5] #0.41

scelta1= cutree(medio,k=6)
coordinate_scelta1= cbind(coordinate_scuole,scelta1)
coordinate_1_1= coordinate_scelta1[which(coordinate_scelta1[3]==1),c(1,2)]
dim(coordinate_1_1) #1311 gruppo 1
coordinate_2_1= coordinate_scelta1[which(coordinate_scelta1[3]==2),c(1,2)]
dim(coordinate_2_1) #38
coordinate_3_1= coordinate_scelta1[which(coordinate_scelta1[3]==3),c(1,2)]
dim(coordinate_3_1) #143
coordinate_4_1= coordinate_scelta1[which(coordinate_scelta1[3]==4),c(1,2)]
dim(coordinate_4_1) #186
coordinate_5_1= coordinate_scelta1[which(coordinate_scelta1[3]==5),c(1,2)]
dim(coordinate_5_1) #95
coordinate_6_1= coordinate_scelta1[which(coordinate_scelta1[3]==6),c(1,2)]
dim(coordinate_6_1) #34

# Rappresentazione geografica
ggmap(get_map(c(left = -73.5, bottom = 41.2, right = -69.9, top = 43),source = "stamen"))+
  geom_point(data = coordinate_1_1, aes(x = lon, y = lat), color ="blue")+
  geom_point(data = coordinate_2_1, aes(x = lon, y = lat), color ="red")+
  geom_point(data = coordinate_3_1, aes(x = lon, y = lat), color ="yellow")+
  geom_point(data = coordinate_4_1, aes(x = lon, y = lat), color ="green")+
  geom_point(data = coordinate_5_1, aes(x = lon, y = lat), color ="purple")+
  geom_point(data = coordinate_6_1, aes(x = lon, y = lat), color ="orange")

# Scelta2: metodo del legame completo con 10 gruppi
II_completo[9]
SI_completo[9]

scelta2= cutree(completo,k=10)
coordinate_scelta2= cbind(coordinate_scuole,scelta2)
coordinate_1_2= coordinate_scelta2[which(coordinate_scelta2[3]==1),c(1,2)]
dim(coordinate_1_2) #519
coordinate_2_2= coordinate_scelta2[which(coordinate_scelta2[3]==2),c(1,2)]
dim(coordinate_2_2) #221
coordinate_3_2= coordinate_scelta2[which(coordinate_scelta2[3]==3),c(1,2)]
dim(coordinate_3_2) #408
coordinate_4_2= coordinate_scelta2[which(coordinate_scelta2[3]==4),c(1,2)]
dim(coordinate_4_2) #191
coordinate_5_2= coordinate_scelta2[which(coordinate_scelta2[3]==5),c(1,2)]
dim(coordinate_5_2) #146
coordinate_6_2= coordinate_scelta2[which(coordinate_scelta2[3]==6),c(1,2)]
dim(coordinate_6_2) #55
coordinate_7_2= coordinate_scelta2[which(coordinate_scelta2[3]==7),c(1,2)]
dim(coordinate_7_2) #137
coordinate_8_2= coordinate_scelta2[which(coordinate_scelta2[3]==8),c(1,2)]
dim(coordinate_8_2) #72
coordinate_9_2= coordinate_scelta2[which(coordinate_scelta2[3]==9),c(1,2)]
dim(coordinate_9_2) #34
coordinate_10_2= coordinate_scelta2[which(coordinate_scelta2[3]==10),c(1,2)]
dim(coordinate_10_2) #24

# Rappresentazione geografica
ggmap(get_map(c(left = -73.5, bottom = 41.2, right = -69.9, top = 43),source = "stamen"))+
  geom_point(data = coordinate_1_2, aes(x = lon, y = lat), color ="blue" , size = 1,na.rm =F)+
  geom_point(data = coordinate_2_2, aes(x = lon, y = lat), color ="green" , size = 1,na.rm =F)+
  geom_point(data = coordinate_3_2, aes(x = lon, y = lat), color ="purple" , size = 1,na.rm =F)+
  geom_point(data = coordinate_4_2, aes(x = lon, y = lat), color ="yellow" , size = 1,na.rm =F)+
  geom_point(data = coordinate_5_2, aes(x = lon, y = lat), color ="red" , size = 1,na.rm =F)+
  geom_point(data = coordinate_6_2, aes(x = lon, y = lat), color ="orange" , size = 1,na.rm =F)+
  geom_point(data = coordinate_7_2, aes(x = lon, y = lat), color ="cyan" , size = 1,na.rm =F)+
  geom_point(data = coordinate_8_2, aes(x = lon, y = lat), color ="black" , size = 1,na.rm =F)+
  geom_point(data = coordinate_9_2, aes(x = lon, y = lat), color ="darkgreen" , size = 1,na.rm =F)+
  geom_point(data = coordinate_10_2, aes(x = lon, y = lat), color ="violet" , size = 1,na.rm =F)


# Attraction-Repulsion clustering
# SCELTA 1: medio con 6 gruppi
dataset=as.data.frame(dataset)
# Matrice delle distanze euclidee tra gli attributi sensibili
CD=daisy(dataset[,4:9])
CD_z=scale(dataset[,4:9])

# Delta 1
# Creo i parametri necessari
v= seq(from=0.1,to=10, by=0.1)
U= matrix(data=0,ncol=6,nrow=6)
U= as.matrix(U)
V= matrix(data = c(1,-1,-1,-1,-1,-1,
                   -1, 1,-1,-1,-1,-1,
                   -1,-1, 1,-1,-1,-1,
                   -1,-1,-1, 1,-1,-1,
                   -1,-1,-1,-1, 1,-1,
                   -1,-1,-1,-1,-1, 1),ncol=6,nrow=6,byrow = T) 
#matrice ovvia in cui cerco di mixare tutte le etnie
V= as.matrix(V)

# Minimizzo l'Imbalance Index
imb_opt=Inf
sol_opt=NA
for(valore in v){
  params = list(U,valore,V)
  funzione_D1=fairDistanceCalc(params,type = "additive1",distances=as.matrix(distanze),
                               chargeMatrix = as.matrix(dataset[,4:9]),Ks=6,totalPropMatrix = p_tot,
                               chargeDistance = as.matrix(CD)) 
  clus_prop_D1= funzione_D1$averageHclust[[1]]$clusterProportions
  addendi=c()
  for(riga in 1:6){
    addendi=c(addendi,sqrt(sum((clus_prop_D1[riga,]-p_tot)^2)))
  }
  imb= sum(addendi)*(1/6)
  if (imb<imb_opt){
    imb_opt=imb
    sol_opt=valore
  }
}
# Con Delta1 si è ottenuta la seguente soluzione ottima
SOL_D1_1=sol_opt #v=7.2
IMB_D1_1=imb_opt #0.1344883

funzione_D1_1=fairDistanceCalc(list(U,SOL_D1-1,V),type = "additive1",distances=as.matrix(distanze),
                             chargeMatrix = as.matrix(dataset[,4:9]),Ks=6,totalPropMatrix = p_tot,
                             chargeDistance = as.matrix(CD)) 
scelta1_D1=funzione_D1_1$clusterAverageHclust
SI_D1_1= summary(silhouette(scelta1_D1[[1]],dist=distanze))$avg.width #0.43
coordinate_scelta1_D1= cbind(coordinate_scuole,scelta1_D1)

coordinate_1_D1_1= coordinate_scelta1_D1[which(coordinate_scelta1_D1[3]==1),c(1,2)]
dim(coordinate_1_D1_1) #1311 gruppo 1
coordinate_2_D1_1= coordinate_scelta1_D1[which(coordinate_scelta1_D1[3]==2),c(1,2)]
dim(coordinate_2_D1_1) #38
coordinate_3_D1_1= coordinate_scelta1_D1[which(coordinate_scelta1_D1[3]==3),c(1,2)]
dim(coordinate_3_D1_1) #164
coordinate_4_D1_1= coordinate_scelta1_D1[which(coordinate_scelta1_D1[3]==4),c(1,2)]
dim(coordinate_4_D1_1) #165
coordinate_5_D1_1= coordinate_scelta1_D1[which(coordinate_scelta1_D1[3]==5),c(1,2)]
dim(coordinate_5_D1_1) #95
coordinate_6_D1_1= coordinate_scelta1_D1[which(coordinate_scelta1_D1[3]==6),c(1,2)]
dim(coordinate_6_D1_1) #34
# Rappresentazione geografica
ggmap(get_map(c(left = -73.5, bottom = 41.2, right = -69.9, top = 43),source = "stamen"))+
  geom_point(data = coordinate_1_D1_1, aes(x = lon, y = lat), color ="blue" , size = 1,na.rm =F)+
  geom_point(data = coordinate_2_D1_1, aes(x = lon, y = lat), color ="red" , size = 1,na.rm =F)+
  geom_point(data = coordinate_3_D1_1, aes(x = lon, y = lat), color ="yellow" , size = 1,na.rm =F)+
  geom_point(data = coordinate_4_D1_1, aes(x = lon, y = lat), color ="green" , size = 1,na.rm =F)+
  geom_point(data = coordinate_5_D1_1, aes(x = lon, y = lat), color ="purple" , size = 1,na.rm =F)+
  geom_point(data = coordinate_6_D1_1, aes(x = lon, y = lat), color ="orange" , size = 1,na.rm =F)

#TENTATIVO CON CD STANDARDIZZATA
v= seq(from=0.1,to=10, by=0.1)
U= matrix(data=0,ncol=6,nrow=6)
U= as.matrix(U)
V= matrix(data = c(1,-1,-1,-1,-1,-1,
                   -1, 1,-1,-1,-1,-1,
                   -1,-1, 1,-1,-1,-1,
                   -1,-1,-1, 1,-1,-1,
                   -1,-1,-1,-1, 1,-1,
                   -1,-1,-1,-1,-1, 1),ncol=6,nrow=6,byrow = T) 
V= as.matrix(V)

# Minimizzo Imbalance Index
imb_opt_z=Inf
sol_opt_z=NA
for(valore in v){
  params = list(U,valore,V)
  matrice_D1=additive1Distance(params, distanceMatrix = as.matrix(distanze), chargeMatrix=as.matrix(dataset[,4:9]) )
  funzione_D1=fairDistanceCalc(params,type = "additive1",distances=as.matrix(distanze),
                               chargeMatrix = as.matrix(dataset[,4:9]),Ks=6,totalPropMatrix = p_tot,
                               chargeDistance = as.matrix(CD_z) )
  clus_prop_D1= funzione_D1$averageHclust[[1]]$clusterProportions
  addendi=c()
  for(riga in 1:6){
    addendi=c(addendi,sqrt(sum((clus_prop_D1[riga,]-p_tot)^2)))
  }
  imb= sum(addendi)*(1/6)
  if (imb<imb_opt_z){
    imb_opt_z=imb
    sol_opt_z=valore
  }
  print(valore)
}
imb_opt_z 
sol_opt_z 
# Lavorare con CD standardizzata o meno porta alla stessa soluzione ottima.
# Scelgo di lavorare con CD non standardizzata da ora in poi


distanze=as.matrix(distanze)
dim(distanze)
# Alcune scuole sono molto vicine e vengono identificate con le stesse coordinate. Di conseguenza
# la distanza risulta nulla anche se le scuole sono distinte.
# Praticamente, si hanno degli 0 anche al di fuori della diagonale della matrice delle distanze.
# Questo crea difficoltà nell'applicazione delle dissimilarità moltiplicative.
# Soluzione pratica: al posto di tutti gli 0 sostituisco 1. In questo modo mantengo l'informazione
# di estrema vicinanza tra le scuole senza dare problemi alle dissimilarità.

distanze_bis=distanze
distanze_bis=as.matrix(distanze_bis)
for (riga in 1:1807){
  for (colonna in 1:1807){
    if (distanze_bis[riga,colonna]==0){
      if (riga!=colonna){
        distanze_bis[riga,colonna]=1
      }
    }
  }
}

# Delta 2
# per delta2 servono i parametri u e v
u= c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
v= c(0.1, 0.3, 0.5, 0.7, 0.9, 1, 3.25, 5.5, 7.75, 10)

# Minimizzo Imbalance Index
imb_opt=Inf
sol_opt=NA
for(valore in v){
  for (valore2 in u){
    params = c(valore2,valore)
    funzione_D2=fairDistanceCalc(params,type = "multiplicative",distances=as.matrix(distanze_bis),
                                 chargeMatrix = as.matrix(dataset[,4:9]),Ks=6,totalPropMatrix = p_tot,
                                 chargeDistance = as.matrix(CD)) 
    clus_prop_D2= funzione_D2$averageHclust[[1]]$clusterProportions
    addendi=c()
    for(riga in 1:6){
      addendi=c(addendi,sqrt(sum((clus_prop_D2[riga,]-p_tot)^2)))
    }
    imb= sum(addendi)*(1/6)
    if (imb<imb_opt){
      imb_opt=imb
      sol_opt=c(valore2,valore)
    }
  }
}

# Con Delta2 si è ottenuta la seguente soluzione ottima
SOL_D2_1=sol_opt #u=0.1  v=0.1
IMB_D2_1=imb_opt #0.1390062

funzione_D2_1=fairDistanceCalc(c(0.1,0.1),type = "multiplicative",distances=as.matrix(distanze_bis),
                             chargeMatrix = as.matrix(dataset[,4:9]),Ks=6,totalPropMatrix = p_tot,
                             chargeDistance = as.matrix(CD)) 
scelta1_D2=funzione_D2_1$clusterAverageHclust
SI_D2_1= summary(silhouette(scelta1_D2[[1]],dist=distanze_prova))$avg.width #0.44
coordinate_scelta1_D2= cbind(coordinate_scuole,scelta1_D2)

coordinate_1_D2_1= coordinate_scelta1_D2[which(coordinate_scelta1_D2[3]==1),c(1,2)]
dim(coordinate_1_D2_1) #1284 gruppo 1
coordinate_2_D2_1= coordinate_scelta1_D2[which(coordinate_scelta1_D2[3]==2),c(1,2)]
dim(coordinate_2_D2_1) #38
coordinate_3_D2_1= coordinate_scelta1_D2[which(coordinate_scelta1_D2[3]==3),c(1,2)]
dim(coordinate_3_D2_1) #164
coordinate_4_D2_1= coordinate_scelta1_D2[which(coordinate_scelta1_D2[3]==4),c(1,2)]
dim(coordinate_4_D2_1) #165
coordinate_5_D2_1= coordinate_scelta1_D2[which(coordinate_scelta1_D2[3]==5),c(1,2)]
dim(coordinate_5_D2_1) #122
coordinate_6_D2_1= coordinate_scelta1_D2[which(coordinate_scelta1_D2[3]==6),c(1,2)]
dim(coordinate_6_D2_1) #34

# Rappresentazione geografica
ggmap(get_map(c(left = -73.5, bottom = 41.2, right = -69.9, top = 43),source = "stamen"))+
  geom_point(data = coordinate_1_D2_1, aes(x = lon, y = lat), color ="blue" , size = 1,na.rm =F)+
  geom_point(data = coordinate_2_D2_1, aes(x = lon, y = lat), color ="red" , size = 1,na.rm =F)+
  geom_point(data = coordinate_3_D2_1, aes(x = lon, y = lat), color ="yellow" , size = 1,na.rm =F)+
  geom_point(data = coordinate_4_D2_1, aes(x = lon, y = lat), color ="green" , size = 1,na.rm =F)+
  geom_point(data = coordinate_5_D2_1, aes(x = lon, y = lat), color ="purple" , size = 1,na.rm =F)+
  geom_point(data = coordinate_6_D2_1, aes(x = lon, y = lat), color ="orange" , size = 1,na.rm =F)


# Delta 4
# Definisco i parametri
u= seq(0.1,1,by=0.2)
v= c(0.1, 1, 10)
w= c(0.1, 0.5, 0.9, 1, 5.5, 10)
V= matrix(data = c(1,-1,-1,-1,-1,-1,
                   -1, 1,-1,-1,-1,-1,
                   -1,-1, 1,-1,-1,-1,
                   -1,-1,-1, 1,-1,-1,
                   -1,-1,-1,-1, 1,-1,
                   -1,-1,-1,-1,-1, 1),ncol=6,nrow=6,byrow = T) 
V= as.matrix(V)

imb_opt=Inf
sol_opt=NA
for(valore in v){
  for (valore2 in u){
    for (valore3 in w){
      params = list(valore2,valore,valore3,V)
      funzione_D4=fairDistanceCalc(params,type = "local",distances=as.matrix(distanze_bis),
                                   chargeMatrix = as.matrix(dataset[,4:9]),Ks=6,totalPropMatrix = p_tot,
                                   chargeDistance = as.matrix(CD)) 
      clus_prop_D4= funzione_D4$averageHclust[[1]]$clusterProportions
      addendi=c()
      for(riga in 1:6){
        addendi=c(addendi,sqrt(sum((clus_prop_D4[riga,]-p_tot)^2)))
      }
      imb= sum(addendi)*(1/6)
      if (imb<imb_opt){
        imb_opt=imb
        sol_opt=c(valore2,valore,valore3)
      }
    }
  }
}  
# Con Delta4 si è ottenuta la seguente soluzione ottima
SOL_D4_1=sol_opt # u=0.1, v=0.1, w=0.1
IMB_D4_1=imb_opt # 0.12

funzione_D4_1=fairDistanceCalc(list(0.1,0.1,0.1,V),type = "local",distances=as.matrix(distanze_prova),
                             chargeMatrix = as.matrix(dataset[,4:9]),Ks=6,totalPropMatrix = p_tot,
                             chargeDistance = as.matrix(CD)) 
scelta1_D4=funzione_D4_1$clusterAverageHclust
SI_D4_1= summary(silhouette(scelta1_D4[[1]],dist=distanze_prova))$avg.width #0.41
coordinate_scelta1_D4= cbind(coordinate_scuole,scelta1_D4)

coordinate_1_D4_1= coordinate_scelta1_D4[which(coordinate_scelta1_D4[3]==1),c(1,2)]
dim(coordinate_1_D4_1) #1311 gruppo 1
coordinate_2_D4_1= coordinate_scelta1_D4[which(coordinate_scelta1_D4[3]==2),c(1,2)]
dim(coordinate_2_D4_1) #38
coordinate_3_D4_1= coordinate_scelta1_D4[which(coordinate_scelta1_D4[3]==3),c(1,2)]
dim(coordinate_3_D4_1) #143
coordinate_4_D4_1= coordinate_scelta1_D4[which(coordinate_scelta1_D4[3]==4),c(1,2)]
dim(coordinate_4_D4_1) #186
coordinate_5_D4_1= coordinate_scelta1_D4[which(coordinate_scelta1_D4[3]==5),c(1,2)]
dim(coordinate_5_D4_1) #95
coordinate_6_D4_1= coordinate_scelta1_D4[which(coordinate_scelta1_D4[3]==6),c(1,2)]
dim(coordinate_6_D4_1) #34


#il totale torna a 1807
ggmap(get_map(c(left = -73.5, bottom = 41.2, right = -69.9, top = 43),source = "stamen"))+
  geom_point(data = coordinate_1_D4_1, aes(x = lon, y = lat), color ="blue" , size = 1,na.rm =F)+
  geom_point(data = coordinate_2_D4_1, aes(x = lon, y = lat), color ="red" , size = 1,na.rm =F)+
  geom_point(data = coordinate_3_D4_1, aes(x = lon, y = lat), color ="yellow" , size = 1,na.rm =F)+
  geom_point(data = coordinate_4_D4_1, aes(x = lon, y = lat), color ="green" , size = 1,na.rm =F)+
  geom_point(data = coordinate_5_D4_1, aes(x = lon, y = lat), color ="purple" , size = 1,na.rm =F)+
  geom_point(data = coordinate_6_D4_1, aes(x = lon, y = lat), color ="orange" , size = 1,na.rm =F)

# Confronto tra le partizioni
adjustedRandIndex(cutree(medio,k=6),scelta1_D1[[1]]) # 0.99
adjustedRandIndex(cutree(medio,k=6),scelta1_D2[[1]]) # 0.95
adjustedRandIndex(cutree(medio,k=6),scelta1_D4[[1]]) # 1

# Ripeto per la scelta2: completo con 10 gruppi
# Delta1
v= seq(from=0.1,to=10, by=0.1)
U= matrix(data=0,ncol=6,nrow=6)
U= as.matrix(U)
V= matrix(data = c(1,-1,-1,-1,-1,-1,
                   -1, 1,-1,-1,-1,-1,
                   -1,-1, 1,-1,-1,-1,
                   -1,-1,-1, 1,-1,-1,
                   -1,-1,-1,-1, 1,-1,
                   -1,-1,-1,-1,-1, 1),ncol=6,nrow=6,byrow = T) 
#matrice ovvia in cui cerco di mixare tutte le etnie
V= as.matrix(V)

imb_opt=Inf
sol_opt=NA
for(valore in v){
  params = list(U,valore,V)
  funzione_D1=fairDistanceCalc(params,type = "additive1",distances=as.matrix(distanze),
                               chargeMatrix = as.matrix(dataset[,4:9]),Ks=10,totalPropMatrix = p_tot,
                               chargeDistance = as.matrix(CD)) 
  clus_prop_D1= funzione_D1$completeHclust[[1]]$clusterProportions
  addendi=c()
  for(riga in 1:10){
    addendi=c(addendi,sqrt(sum((clus_prop_D1[riga,]-p_tot)^2)))
  }
  imb= sum(addendi)*(1/10)
  if (imb<imb_opt){
    imb_opt=imb
    sol_opt=valore
  }
}
# Con Delta1 si è ottenuta la seguente soluzione ottima
SOL_D1_2=sol_opt #v=8.3
IMB_D1_2=imb_opt #0.14


funzione_D1_2=fairDistanceCalc(list(U,SOL_D1_2,V),type = "additive1",distances=as.matrix(distanze),
                                 chargeMatrix = as.matrix(dataset[,4:9]),Ks=10,totalPropMatrix = p_tot,
                                 chargeDistance = as.matrix(CD)) 
scelta2_D1=funzione_D1_2$clusterCompleteHclust
SI_D1_2= summary(silhouette(scelta2_D1[[1]],dist=distanze))$avg.width #0.38
coordinate_scelta2_D1= cbind(coordinate_scuole,scelta2_D1)

coordinate_1_D1_2= coordinate_scelta2_D1[which(coordinate_scelta2_D1[3]==1),c(1,2)]
dim(coordinate_1_D1_2) #519 gruppo 1
coordinate_2_D1_2= coordinate_scelta2_D1[which(coordinate_scelta2_D1[3]==2),c(1,2)]
dim(coordinate_2_D1_2) #239
coordinate_3_D1_2= coordinate_scelta2_D1[which(coordinate_scelta2_D1[3]==3),c(1,2)]
dim(coordinate_3_D1_2) #589
coordinate_4_D1_2= coordinate_scelta2_D1[which(coordinate_scelta2_D1[3]==4),c(1,2)]
dim(coordinate_4_D1_2) #55
coordinate_5_D1_2= coordinate_scelta2_D1[which(coordinate_scelta2_D1[3]==5),c(1,2)]
dim(coordinate_5_D1_2) #128
coordinate_6_D1_2= coordinate_scelta2_D1[which(coordinate_scelta2_D1[3]==6),c(1,2)]
dim(coordinate_6_D1_2) #137
coordinate_7_D1_2= coordinate_scelta2_D1[which(coordinate_scelta2_D1[3]==7),c(1,2)]
dim(coordinate_7_D1_2) #85
coordinate_8_D1_2= coordinate_scelta2_D1[which(coordinate_scelta2_D1[3]==8),c(1,2)]
dim(coordinate_8_D1_2) #24
coordinate_9_D1_2= coordinate_scelta2_D1[which(coordinate_scelta2_D1[3]==9),c(1,2)]
dim(coordinate_9_D1_2) #19
coordinate_10_D1_2= coordinate_scelta2_D1[which(coordinate_scelta2_D1[3]==10),c(1,2)]
dim(coordinate_10_D1_2) #12

# Rappresentazione geografica
ggmap(get_map(c(left = -73.5, bottom = 41.2, right = -69.9, top = 43),source = "stamen"))+
  geom_point(data = coordinate_1_D1_2, aes(x = lon, y = lat), color ="blue" , size = 1,na.rm =F)+
  geom_point(data = coordinate_2_D1_2, aes(x = lon, y = lat), color ="green" , size = 1,na.rm =F)+
  geom_point(data = coordinate_3_D1_2, aes(x = lon, y = lat), color ="purple" , size = 1,na.rm =F)+
  geom_point(data = coordinate_4_D1_2, aes(x = lon, y = lat), color ="orange" , size = 1,na.rm =F)+
  geom_point(data = coordinate_5_D1_2, aes(x = lon, y = lat), color ="yellow" , size = 1,na.rm =F)+
  geom_point(data = coordinate_6_D1_2, aes(x = lon, y = lat), color ="cyan" , size = 1,na.rm =F)+
  geom_point(data = coordinate_7_D1_2, aes(x = lon, y = lat), color ="red" , size = 1,na.rm =F)+
  geom_point(data = coordinate_8_D1_2, aes(x = lon, y = lat), color ="violet" , size = 1,na.rm =F)+
  geom_point(data = coordinate_9_D1_2, aes(x = lon, y = lat), color ="darkgreen" , size = 1,na.rm =F)+
  geom_point(data = coordinate_10_D1_2, aes(x = lon, y = lat), color ="black" , size = 1,na.rm =F)

# Delta2
u= c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
v= c(0.1, 0.3, 0.5, 0.7, 0.9, 1, 3.25, 5.5, 7.75, 10)

imb_opt=Inf
sol_opt=NA
for(valore in v){
  for (valore2 in u){
    params = c(valore2,valore)
    funzione_D2=fairDistanceCalc(params,type = "multiplicative",distances=as.matrix(distanze_prova),
                                 chargeMatrix = as.matrix(dataset[,4:9]),Ks=10,totalPropMatrix = p_tot,
                                 chargeDistance = as.matrix(CD)) 
    clus_prop_D2= funzione_D2$completeHclust[[1]]$clusterProportions
    addendi=c()
    for(riga in 1:10){
      addendi=c(addendi,sqrt(sum((clus_prop_D2[riga,]-p_tot)^2)))
    }
    imb= sum(addendi)*(1/10)
    if (imb<imb_opt){
      imb_opt=imb
      sol_opt=c(valore2,valore)
    }
  }
}


# Con Delta2 si è ottenuta la seguente soluzione ottima
SOL_D2_2=sol_opt # u=1, v=0.1
IMB_D2_2=imb_opt # 0.15

funzione_D2_2=fairDistanceCalc(SOL_D2_2,type = "multiplicative",distances=as.matrix(distanze_prova),
                                 chargeMatrix = as.matrix(dataset[,4:9]),Ks=10,totalPropMatrix = p_tot,
                                 chargeDistance = as.matrix(CD)) 
scelta2_D2=funzione_D2_2$clusterCompleteHclust
SI_D2_2= summary(silhouette(scelta2_D2[[1]],dist=distanze_prova))$avg.width #0.40
coordinate_scelta2_D2= cbind(coordinate_scuole,scelta2_D2)

coordinate_1_D2_2= coordinate_scelta2_D2[which(coordinate_scelta2_D2[3]==1),c(1,2)]
dim(coordinate_1_D2_2) #647
coordinate_2_D2_2= coordinate_scelta2_D2[which(coordinate_scelta2_D2[3]==2),c(1,2)]
dim(coordinate_2_D2_2) #179
coordinate_3_D2_2= coordinate_scelta2_D2[which(coordinate_scelta2_D2[3]==3),c(1,2)]
dim(coordinate_3_D2_2) #487
coordinate_4_D2_2= coordinate_scelta2_D2[which(coordinate_scelta2_D2[3]==4),c(1,2)]
dim(coordinate_4_D2_2) #55
coordinate_5_D2_2= coordinate_scelta2_D2[which(coordinate_scelta2_D2[3]==5),c(1,2)]
dim(coordinate_5_D2_2) #172
coordinate_6_D2_2= coordinate_scelta2_D2[which(coordinate_scelta2_D2[3]==6),c(1,2)]
dim(coordinate_6_D2_2) #137
coordinate_7_D2_2= coordinate_scelta2_D2[which(coordinate_scelta2_D2[3]==7),c(1,2)]
dim(coordinate_7_D2_2) #72
coordinate_8_D2_2= coordinate_scelta2_D2[which(coordinate_scelta2_D2[3]==8),c(1,2)]
dim(coordinate_8_D2_2) #30
coordinate_9_D2_2= coordinate_scelta2_D2[which(coordinate_scelta2_D2[3]==9),c(1,2)]
dim(coordinate_9_D2_2) #24
coordinate_10_D2_2= coordinate_scelta2_D2[which(coordinate_scelta2_D2[3]==10),c(1,2)]
dim(coordinate_10_D2_2) #4

# Rappresentazione geografica
ggmap(get_map(c(left = -73.5, bottom = 41.2, right = -69.9, top = 43),source = "stamen"))+
  geom_point(data = coordinate_1_D2_2, aes(x = lon, y = lat), color ="blue" , size = 1,na.rm =F)+
  geom_point(data = coordinate_2_D2_2, aes(x = lon, y = lat), color ="green" , size = 1,na.rm =F)+
  geom_point(data = coordinate_3_D2_2, aes(x = lon, y = lat), color ="purple" , size = 1,na.rm =F)+
  geom_point(data = coordinate_4_D2_2, aes(x = lon, y = lat), color ="orange" , size = 1,na.rm =F)+
  geom_point(data = coordinate_5_D2_2, aes(x = lon, y = lat), color ="yellow" , size = 1,na.rm =F)+
  geom_point(data = coordinate_6_D2_2, aes(x = lon, y = lat), color ="cyan" , size = 1,na.rm =F)+
  geom_point(data = coordinate_7_D2_2, aes(x = lon, y = lat), color ="red" , size = 1,na.rm =F)+
  geom_point(data = coordinate_8_D2_2, aes(x = lon, y = lat), color ="black" , size = 1,na.rm =F)+
  geom_point(data = coordinate_9_D2_2, aes(x = lon, y = lat), color ="violet" , size = 1,na.rm =F)+
  geom_point(data = coordinate_10_D2_2, aes(x = lon, y = lat), color ="darkgreen" , size = 1,na.rm =F)

# Delta 4
u= seq(0.1,1,by=0.2)
v= c(0.1, 1, 10)
w= c(0.1, 0.5, 0.9, 1, 5.5, 10)
str(v)

V= matrix(data = c(1,-1,-1,-1,-1,-1,
                   -1, 1,-1,-1,-1,-1,
                   -1,-1, 1,-1,-1,-1,
                   -1,-1,-1, 1,-1,-1,
                   -1,-1,-1,-1, 1,-1,
                   -1,-1,-1,-1,-1, 1),ncol=6,nrow=6,byrow = T) #matrice ovvia in cui cerco di mix tutte le etnie
V= as.matrix(V)
dim(V)


imb_opt=Inf
sol_opt=NA
for(valore in v){
  for (valore2 in u){
    for (valore3 in w){
      params = list(valore2,valore,valore3,V)
      funzione_D4=fairDistanceCalc(params,type = "local",distances=as.matrix(distanze_prova),
                                   chargeMatrix = as.matrix(dataset[,4:9]),Ks=10,totalPropMatrix = p_tot,
                                   chargeDistance = as.matrix(CD)) 
      clus_prop_D4= funzione_D4$completeHclust[[1]]$clusterProportions
      addendi=c()
      for(riga in 1:10){
        addendi=c(addendi,sqrt(sum((clus_prop_D4[riga,]-p_tot)^2)))
      }
      imb= sum(addendi)*(1/10)
      print(imb)
      if (imb<imb_opt){
        imb_opt=imb
        sol_opt=c(valore2,valore,valore3)
      }
      print(c(valore2,valore,valore3))
    }
  }
}  

# Con Delta4 si è ottenuta la seguente soluzione ottima
SOL_D4_2=sol_opt #u=0.1, v=0.1, w=0.1
IMB_D4_2=imb_opt #0.16

funzione_D4_2=fairDistanceCalc(list(...,V),type = "local",distances=as.matrix(distanze_prova),
                                 chargeMatrix = as.matrix(dataset[,4:9]),Ks=10,totalPropMatrix = p_tot,
                                 chargeDistance = as.matrix(CD)) 
scelta2_D4=funzione_D4_2$clusterCompleteHclust
SI_D4_2= summary(silhouette(scelta2_D4[[1]],dist=distanze_prova))$avg.width #0.43
coordinate_scelta2_D4= cbind(coordinate_scuole,scelta2_D4)

coordinate_1_D4_2= coordinate_scelta2_D4[which(coordinate_scelta2_D4[3]==1),c(1,2)]
dim(coordinate_1_D4_2) #519
coordinate_2_D4_2= coordinate_scelta2_D4[which(coordinate_scelta2_D4[3]==2),c(1,2)]
dim(coordinate_2_D4_2) #221
coordinate_3_D4_2= coordinate_scelta2_D4[which(coordinate_scelta2_D4[3]==3),c(1,2)]
dim(coordinate_3_D4_2) #408
coordinate_4_D4_2= coordinate_scelta2_D4[which(coordinate_scelta2_D4[3]==4),c(1,2)]
dim(coordinate_4_D4_2) #191
coordinate_5_D4_2= coordinate_scelta2_D4[which(coordinate_scelta2_D4[3]==5),c(1,2)]
dim(coordinate_5_D4_2) #146
coordinate_6_D4_2= coordinate_scelta2_D4[which(coordinate_scelta2_D4[3]==6),c(1,2)]
dim(coordinate_6_D4_2) #55
coordinate_7_D4_2= coordinate_scelta2_D4[which(coordinate_scelta2_D4[3]==7),c(1,2)]
dim(coordinate_7_D2_2) #137
coordinate_8_D4_2= coordinate_scelta2_D4[which(coordinate_scelta2_D4[3]==8),c(1,2)]
dim(coordinate_8_D2_2) #72
coordinate_9_D4_2= coordinate_scelta2_D4[which(coordinate_scelta2_D4[3]==9),c(1,2)]
dim(coordinate_9_D4_2) #34
coordinate_10_D4_2= coordinate_scelta2_D4[which(coordinate_scelta2_D4[3]==10),c(1,2)]
dim(coordinate_10_D4_2) #24

# Rappresentazione geografica
ggmap(get_map(c(left = -73.5, bottom = 41.2, right = -69.9, top = 43),source = "stamen"))+
  geom_point(data = coordinate_1_D4_2, aes(x = lon, y = lat), color ="blue" , size = 1,na.rm =F)+
  geom_point(data = coordinate_2_D4_2, aes(x = lon, y = lat), color ="green" , size = 1,na.rm =F)+
  geom_point(data = coordinate_3_D4_2, aes(x = lon, y = lat), color ="purple" , size = 1,na.rm =F)+
  geom_point(data = coordinate_4_D4_2, aes(x = lon, y = lat), color ="yellow" , size = 1,na.rm =F)+
  geom_point(data = coordinate_5_D4_2, aes(x = lon, y = lat), color ="red" , size = 1,na.rm =F)+
  geom_point(data = coordinate_6_D4_2, aes(x = lon, y = lat), color ="orange" , size = 1,na.rm =F)+
  geom_point(data = coordinate_7_D4_2, aes(x = lon, y = lat), color ="cyan" , size = 1,na.rm =F)+
  geom_point(data = coordinate_8_D4_2, aes(x = lon, y = lat), color ="black" , size = 1,na.rm =F)+
  geom_point(data = coordinate_9_D4_2, aes(x = lon, y = lat), color ="darkgreen" , size = 1,na.rm =F)+
  geom_point(data = coordinate_10_D4_2, aes(x = lon, y = lat), color ="violet" , size = 1,na.rm =F)

# Confronto tra le partizioni
adjustedRandIndex(cutree(completo,k=10),scelta2_D1[[1]]) # 0.79
adjustedRandIndex(cutree(completo,k=10),scelta2_D2[[1]]) # 0.61
adjustedRandIndex(cutree(completo,k=10),scelta2_D4[[1]]) # 1
