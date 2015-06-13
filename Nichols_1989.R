


m<-0.125  #dispersal per deme
s<-0.02   #selection against hybrids

##gradients
steep1<-c(20,12,12,6,6,6,2,6,6,6,12,12,20)
steep2<-c(50,30,30,15,15,15,5,15,15,15,30,30,50)
steep3<-c(100, 60, 60, 30, 30, 30, 10, 30, 30, 30, 60, 60, 100)
steep3.1<-c(100, 60, 60, 30, 30, 30, 30, 30, 10, 30, 30, 30, 100)
steep4<-c(500,300,300,150,150,150,50,150,150,150,300,300,500)
shallow1<-c(10,4,4,3,3,3,2,3,3,3,4,4,10)
shallow2<-c(250,100,100,75,75,75,50,75,75,75,100,100,250)

initialise<-function(gradient=steep3){
        
        ##abundance (density) matrix
        abu<-matrix(c(0,gradient,0), 15,20)
        abu[,1]<-0
        abu[,20]<-0
#        abu


        ##number of neighbours
        neigh<-matrix(4,13,18)
        neigh[1,]<-3
        neigh[13,]<-3
        neigh[,1]<-3
        neigh[,18]<-3
        neigh[1,1]<-2
        neigh[1,18]<-2
        neigh[13,1]<-2
        neigh[13,18]<-2
#        neigh
        neigh2<-matrix(0,15,20)
        neigh2[2:14,2:19]<-neigh2[2:14,2:19]+neigh
        neigh2
        ##initial allele ditribution
        
        pop<-abu*2*matrix(c(0,1,1,1,1,1,1,0.5,0,0,0,0,0,0,0), 15, 20)
#        pop

        image(t(pop/abu), asp=15/18, axes=F)
#        pop/abu
        return(list(abu,pop))
}




tgb<-function(generations=1000){        ##time goes by
        
        for(i in 1:generations){
                #dispersal
                disp<-matrix(0,15,20)
                #dispn<-matrix(0,15,20)
                
                #dispn[2:14,2:19]<-((1-m)*pop[2:14,2:19] + (pop[1:13,1:18] + pop[3:15,1:18] + pop[1:13,3:20] + pop[3:15,3:20])*/neigh2[2:14,2:19]*m)/
                #  ((1-m)*abu[2:14,2:19] + (abu[1:13,1:18] + abu[3:15,1:18] + abu[1:13,3:20] + abu[3:15,3:20])/neigh2[2:14,2:19]*m)
                
                disp[2:14,2:19]<-((1-m)*abupop[[2]][2:14,2:19] + (abupop[[2]][1:13,2:19] + abupop[[2]][3:15,2:19] + abupop[[2]][2:14,1:18] + abupop[[2]][2:14,3:20])*m/4)/
                        ((1-m)*abupop[[1]][2:14,2:19]*2 + (abupop[[1]][1:13,2:19] + abupop[[1]][3:15,2:19] + abupop[[1]][2:14,1:18] + abupop[[1]][2:14,3:20])*m*2/4)
                
                
                
                #round(disp, digits=2)
                #round(dispn, digits=2)
                #selection
                
                sel<-matrix(0,15,20)
                
                sel[2:14,2:19]<-(disp[2:14,2:19]**2 + disp[2:14,2:19] * (1 - disp[2:14,2:19]) * (1 - s))/
                        (disp[2:14,2:19]**2 + 2 * disp[2:14,2:19] * (1 - disp[2:14,2:19]) * (1 - s) + (1 - disp[2:14,2:19])**2)
                
                
                
                ##selection dummy
                #sel<-disp
                #round(sel, digits=2)
                ##drift
                dri<-matrix(0,15,20)
                for( i in 2:14){
                        for(j in 2: 19){
                                dri[i,j]<-rbinom(1,abupop[[1]][i,j]*2,sel[i,j])
                        }
                }
                ##drift dummy (no drift)
                #dri<-round(sel*abu*2)
                
                
                
                abupop[[2]]<-dri
        }
        #pop
        image(t(abupop[[2]]/abupop[[1]]), asp=15/18, axes=F)
        return(abupop)
}
 

abupop<-initialise(shallow1)


##run a specifies number of generations
abupop<-tgb(100)


#or loop and watch 

for(i in 1:1000){
        print(i)
        abupop<-tgb(100)
}

  

