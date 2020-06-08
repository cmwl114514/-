simulation_MCMC_complete<-function(loop=100,step=10^4,
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
    n_obs<-apply(sim_data, 1, sum)
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
    Nxy<-matrix(0,2,16)
    for (i in 1:burnin) {
      for (j in 2:16) {
        Nxy[,j]<-rmultinom(1,sim_data[1,j],Pxy[,j])
      }
      P0<-Pxy[1,1]+Pxy[2,1]
      N1<-n_obs[1]+rnbinom(1,n_obs[1],1-P0)
      Nxy[,1]<-rmultinom(1,N1-n_obs[1],Pxy[,1])
      
      Px<-rbeta(1,1+sum(Nxy[1,]),1+sum(Nxy[2,]))
      nyx<-matrix(0,2,4)
      for (jA in 0:1) {
        for (jB in 0:1) {
          for (jC in 0:1) {
            for (jD in 0:1) {
              ind<-1+jA+2*jB+4*jC+8*jD
              nyx[,1]<-nyx[,1]+jA*Nxy[,ind]
              nyx[,2]<-nyx[,2]+jB*Nxy[,ind]
              nyx[,3]<-nyx[,3]+jC*Nxy[,ind]
              nyx[,4]<-nyx[,4]+jD*Nxy[,ind]
            }
          }
        }
      }
      
      for (x in 1:2) {
        for(j in 1:4){
          Pyx[x,j]<-rbeta(1,1+nyx[x,j],1+sum(Nxy[x,])-nyx[x,j])
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
    }
    
    for (i in 1:step) {
      for (j in 2:16) {
        Nxy[,j]<-rmultinom(1,sim_data[1,j],Pxy[,j])
      }
      P0<-Pxy[1,1]+Pxy[2,1]
      N1<-n_obs[1]+rnbinom(1,n_obs[1],1-P0)
      Nxy[,1]<-rmultinom(1,N1-n_obs[1],Pxy[,1])
      
      Px<-rbeta(1,1+sum(Nxy[1,]),1+sum(Nxy[2,]))
      nyx<-matrix(0,2,4)
      for (jA in 0:1) {
        for (jB in 0:1) {
          for (jC in 0:1) {
            for (jD in 0:1) {
              ind<-1+jA+2*jB+4*jC+8*jD
              nyx[,1]<-nyx[,1]+jA*Nxy[,ind]
              nyx[,2]<-nyx[,2]+jB*Nxy[,ind]
              nyx[,3]<-nyx[,3]+jC*Nxy[,ind]
              nyx[,4]<-nyx[,4]+jD*Nxy[,ind]
            }
          }
        }
      }
      
      for (x in 1:2) {
        for(j in 1:4){
          Pyx[x,j]<-rbeta(1,1+nyx[x,j],1+sum(Nxy[x,])-nyx[x,j])
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
      
      N_list0[i]<-N1
    }
    
    N_list[t]<-N_list[t]+round(mean(N_list0))
    
    N_list0<-vector(mode = "numeric",length = step)
    
    Pyx<-matrix(runif(6,0,1),2,3)
    Px<-runif(1,0,1)
    Pxy<-matrix(0,2,8)
    for (jA in 0:1) {
      for (jB in 0:1) {
        for (jC in 0:1) {
          ind<-1+jA+2*jB+4*jC
          Pxy[1,ind]<-Px*ifelse(jA==1,Pyx[1,1],1-Pyx[1,1])*ifelse(jB==1,Pyx[1,2],1-Pyx[1,2])*ifelse(jC==1,Pyx[1,3],1-Pyx[1,3])
          Pxy[2,ind]<-(1-Px)*ifelse(jA==1,Pyx[2,1],1-Pyx[2,1])*ifelse(jB==1,Pyx[2,2],1-Pyx[2,2])*ifelse(jC==1,Pyx[2,3],1-Pyx[2,3])
        }
      }
    }
    Nxy<-matrix(0,2,8)
    for (i in 1:burnin) {
      for (j in 2:8) {
        Nxy[,j]<-rmultinom(1,sim_data[2,j],Pxy[,j])
      }
      P0<-Pxy[1,1]+Pxy[2,1]
      N1<-n_obs[2]+rnbinom(1,n_obs[2],1-P0)
      Nxy[,1]<-rmultinom(1,N1-n_obs[2],Pxy[,1])
      
      Px<-rbeta(1,1+sum(Nxy[1,]),1+sum(Nxy[2,]))
      nyx<-matrix(0,2,3)
      for (jA in 0:1) {
        for (jB in 0:1) {
          for (jC in 0:1) {
              ind<-1+jA+2*jB+4*jC
              nyx[,1]<-nyx[,1]+jA*Nxy[,ind]
              nyx[,2]<-nyx[,2]+jB*Nxy[,ind]
              nyx[,3]<-nyx[,3]+jC*Nxy[,ind]
          }
        }
      }
      
      for (x in 1:2) {
        for(j in 1:3){
          Pyx[x,j]<-rbeta(1,1+nyx[x,j],1+sum(Nxy[x,])-nyx[x,j])
        }
      }
      
      for (jA in 0:1) {
        for (jB in 0:1) {
          for (jC in 0:1) {
              ind<-1+jA+2*jB+4*jC
              Pxy[1,ind]<-Px*ifelse(jA==1,Pyx[1,1],1-Pyx[1,1])*ifelse(jB==1,Pyx[1,2],1-Pyx[1,2])*ifelse(jC==1,Pyx[1,3],1-Pyx[1,3])
              Pxy[2,ind]<-(1-Px)*ifelse(jA==1,Pyx[2,1],1-Pyx[2,1])*ifelse(jB==1,Pyx[2,2],1-Pyx[2,2])*ifelse(jC==1,Pyx[2,3],1-Pyx[2,3])
          }
        }
      }
    }
    
    for (i in 1:step) {
      for (j in 2:8) {
        Nxy[,j]<-rmultinom(1,sim_data[2,j],Pxy[,j])
      }
      P0<-Pxy[1,1]+Pxy[2,1]
      N1<-n_obs[2]+rnbinom(1,n_obs[2],1-P0)
      Nxy[,1]<-rmultinom(1,N1-n_obs[2],Pxy[,1])
      
      Px<-rbeta(1,1+sum(Nxy[1,]),1+sum(Nxy[2,]))
      nyx<-matrix(0,2,3)
      for (jA in 0:1) {
        for (jB in 0:1) {
          for (jC in 0:1) {
            ind<-1+jA+2*jB+4*jC
            nyx[,1]<-nyx[,1]+jA*Nxy[,ind]
            nyx[,2]<-nyx[,2]+jB*Nxy[,ind]
            nyx[,3]<-nyx[,3]+jC*Nxy[,ind]
          }
        }
      }
      
      for (x in 1:2) {
        for(j in 1:3){
          Pyx[x,j]<-rbeta(1,1+nyx[x,j],1+sum(Nxy[x,])-nyx[x,j])
        }
      }
      
      for (jA in 0:1) {
        for (jB in 0:1) {
          for (jC in 0:1) {
            ind<-1+jA+2*jB+4*jC
            Pxy[1,ind]<-Px*ifelse(jA==1,Pyx[1,1],1-Pyx[1,1])*ifelse(jB==1,Pyx[1,2],1-Pyx[1,2])*ifelse(jC==1,Pyx[1,3],1-Pyx[1,3])
            Pxy[2,ind]<-(1-Px)*ifelse(jA==1,Pyx[2,1],1-Pyx[2,1])*ifelse(jB==1,Pyx[2,2],1-Pyx[2,2])*ifelse(jC==1,Pyx[2,3],1-Pyx[2,3])
          }
        }
      }
      
      N_list0[i]<-N1
    }
    
    N_list[t]<-N_list[t]+round(mean(N_list0))
  }
  return(N_list)
}