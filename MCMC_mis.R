simulation_MCMC_mis<-function(loop=100,step=10^4,
                                   burnin=10^3,
                                   N=3*10^3,
                                   ProbA=c(0.1,0.7),
                                   ProbB=c(0.7,0.4),
                                   ProbC=c(0.6,0.2),
                                   ProbD=c(0.3,0.9),
                                   ProbX=0.6){
  N_list<-vector(mode = "numeric",length = loop)
  for (t in 1:loop) {
    sim_data<-simulation(N=N,ProbA=ProbA,ProbB=ProbB,ProbC=ProbC,ProbD=ProbD,ProbX=ProbX)
    n_obs<-sum(sim_data)
    N_list0<-vector(mode = "numeric",length = step)
    
    Pyx<-matrix(runif(8,0,1),2,4)
    Px<-runif(1,0,1)
    Pxy<-matrix(0,2,16)
    for (jA in 0:1) {
      for (jB in 0:1) {
        for (jC in 0:1) {
          for (jD in 0:1) {
            ind<-1+jA+2*jB+4*jC+8*jD
            Pxy[1,ind]<-Px*ifelse(jA==1,Pyx[1,1],1-Pyx[1,1])*ifelse(jB==1,Pyx[1,2],1-Pyx[1,2])*ifelse(jC==1,Pyx[1,3],1-Pyx[1,3])*ifelse(jD==1,Pyx[1,4],1-Pyx[1,4])
            Pxy[2,ind]<-(1-Px)*ifelse(jA==1,Pyx[2,1],1-Pyx[2,1])*ifelse(jB==1,Pyx[2,2],1-Pyx[2,2])*ifelse(jC==1,Pyx[2,3],1-Pyx[2,3])*ifelse(jD==1,Pyx[2,4],1-Pyx[2,4])
          }
        }
      }
    }
    Ps<-sum(sim_data[1,])/n_obs
    Pz1<-c(Px,1-Px)%*%((1-Pyx[,1])*(1-Pyx[,2])*(1-Pyx[,3])*(1-Pyx[,4]))*Ps
    Pz2<-c(Px,1-Px)%*%((1-Pyx[,1])*(1-Pyx[,2])*(1-Pyx[,3]))*(1-Ps)
    
    Nxy1<-Nxy2<-matrix(0,2,16)
    for (i in 1:burnin) {
      for (j in 2:16) {
        Nxy1[,j]<-rmultinom(1,sim_data[1,j],Pxy[,j])
      }
      for (j in 2:8) {
        P_tmp<-c(Pxy[1,j],Pxy[2,j],Pxy[1,j+8],Pxy[2,j+8])
        Nxy2[,c(j,j+8)]<-rmultinom(1,sim_data[2,j],P_tmp)
      }
      
      N_tmp<-n_obs+rnbinom(1,n_obs,1-Pz1-Pz2)
      Nz<-rmultinom(1,N_tmp-n_obs,c(Pz1,Pz2))
      Nxy1[,1]<-rmultinom(1,Nz[1],Pxy[,1])
      P_tmp<-c(Pxy[1,1],Pxy[2,1],Pxy[1,9],Pxy[2,9])
      Nxy2[,c(1,9)]<-rmultinom(1,Nz[2],P_tmp)
      
      Px<-rbeta(1,1+sum(Nxy1[1,])+sum(Nxy2[1,]),1+sum(Nxy1[2,])+sum(Nxy2[2,]))
      nyx<-matrix(0,2,4)
      for (jA in 0:1) {
        for (jB in 0:1) {
          for (jC in 0:1) {
            for (jD in 0:1) {
              ind<-1+jA+2*jB+4*jC+8*jD
              nyx[,1]<-nyx[,1]+jA*(Nxy1[,ind]+Nxy2[,ind])
              nyx[,2]<-nyx[,2]+jB*(Nxy1[,ind]+Nxy2[,ind])
              nyx[,3]<-nyx[,3]+jC*(Nxy1[,ind]+Nxy2[,ind])
              nyx[,4]<-nyx[,4]+jD*(Nxy1[,ind]+Nxy2[,ind])
            }
          }
        }
      }
      for (x in 1:2) {
        for(j in 1:4){
          Pyx[x,j]<-rbeta(1,1+nyx[x,j],1+sum(Nxy1[x,])+sum(Nxy2[x,])-nyx[x,j])
        }
      }
      for (jA in 0:1) {
        for (jB in 0:1) {
          for (jC in 0:1) {
            for (jD in 0:1) {
              ind<-1+jA+2*jB+4*jC+8*jD
              Pxy[1,ind]<-Px*ifelse(jA==1,Pyx[1,1],1-Pyx[1,1])*ifelse(jB==1,Pyx[1,2],1-Pyx[1,2])*ifelse(jC==1,Pyx[1,3],1-Pyx[1,3])*ifelse(jD==1,Pyx[1,4],1-Pyx[1,4])
              Pxy[2,ind]<-(1-Px)*ifelse(jA==1,Pyx[2,1],1-Pyx[2,1])*ifelse(jB==1,Pyx[2,2],1-Pyx[2,2])*ifelse(jC==1,Pyx[2,3],1-Pyx[2,3])*ifelse(jD==1,Pyx[2,4],1-Pyx[2,4])
            }
          }
        }
      }
      
      Ps<-sum(Nxy1)/N_tmp
      Pz1<-c(Px,1-Px)%*%((1-Pyx[,1])*(1-Pyx[,2])*(1-Pyx[,3])*(1-Pyx[,4]))*Ps
      Pz2<-c(Px,1-Px)%*%((1-Pyx[,1])*(1-Pyx[,2])*(1-Pyx[,3]))*(1-Ps)
    }
    
    for (i in 1:step) {
      for (j in 2:16) {
        Nxy1[,j]<-rmultinom(1,sim_data[1,j],Pxy[,j])
      }
      for (j in 2:8) {
        P_tmp<-c(Pxy[1,j],Pxy[2,j],Pxy[1,j+8],Pxy[2,j+8])
        Nxy2[,c(j,j+8)]<-rmultinom(1,sim_data[2,j],P_tmp)
      }
      
      N_tmp<-n_obs+rnbinom(1,n_obs,1-Pz1-Pz2)
      Nz<-rmultinom(1,N_tmp-n_obs,c(Pz1,Pz2))
      Nxy1[,1]<-rmultinom(1,Nz[1],Pxy[,1])
      P_tmp<-c(Pxy[1,1],Pxy[2,1],Pxy[1,9],Pxy[2,9])
      Nxy2[,c(1,9)]<-rmultinom(1,Nz[2],P_tmp)
      
      Px<-rbeta(1,1+sum(Nxy1[1,])+sum(Nxy2[1,]),1+sum(Nxy1[2,])+sum(Nxy2[2,]))
      nyx<-matrix(0,2,4)
      for (jA in 0:1) {
        for (jB in 0:1) {
          for (jC in 0:1) {
            for (jD in 0:1) {
              ind<-1+jA+2*jB+4*jC+8*jD
              nyx[,1]<-nyx[,1]+jA*(Nxy1[,ind]+Nxy2[,ind])
              nyx[,2]<-nyx[,2]+jB*(Nxy1[,ind]+Nxy2[,ind])
              nyx[,3]<-nyx[,3]+jC*(Nxy1[,ind]+Nxy2[,ind])
              nyx[,4]<-nyx[,4]+jD*(Nxy1[,ind]+Nxy2[,ind])
            }
          }
        }
      }
      for (x in 1:2) {
        for(j in 1:4){
          Pyx[x,j]<-rbeta(1,1+nyx[x,j],1+sum(Nxy1[x,])+sum(Nxy2[x,])-nyx[x,j])
        }
      }
      for (jA in 0:1) {
        for (jB in 0:1) {
          for (jC in 0:1) {
            for (jD in 0:1) {
              ind<-1+jA+2*jB+4*jC+8*jD
              Pxy[1,ind]<-Px*ifelse(jA==1,Pyx[1,1],1-Pyx[1,1])*ifelse(jB==1,Pyx[1,2],1-Pyx[1,2])*ifelse(jC==1,Pyx[1,3],1-Pyx[1,3])*ifelse(jD==1,Pyx[1,4],1-Pyx[1,4])
              Pxy[2,ind]<-(1-Px)*ifelse(jA==1,Pyx[2,1],1-Pyx[2,1])*ifelse(jB==1,Pyx[2,2],1-Pyx[2,2])*ifelse(jC==1,Pyx[2,3],1-Pyx[2,3])*ifelse(jD==1,Pyx[2,4],1-Pyx[2,4])
            }
          }
        }
      }
      
      Ps<-sum(Nxy1)/N_tmp
      Pz1<-c(Px,1-Px)%*%((1-Pyx[,1])*(1-Pyx[,2])*(1-Pyx[,3])*(1-Pyx[,4]))*Ps
      Pz2<-c(Px,1-Px)%*%((1-Pyx[,1])*(1-Pyx[,2])*(1-Pyx[,3]))*(1-Ps)
      
      N_list0[i]<-N_tmp
    }
    N_list[t]<-round(mean(N_list0))
  }
  return(N_list)
}