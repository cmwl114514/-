simulation_EM<-function(loop=200,step=60,
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
    if(sum(sim_data)==0) next()
    piA<-piB<-piC<-piD<-0.5
    n_com<-sim_data
    if(loop==1){
      N_list0<-vector(mode = "numeric",length = step)
    }
    for (s in 1:step) {
      tmp<-(1-piA)*(1-piB)*(1-piC)*(1-piD)
      n_com[1,1]<-tmp*n_obs[1]/(1-tmp)
      
      tmp<-(1-piA)*(1-piB)*(1-piC)
      n_com[2,1]<-tmp*n_obs[2]/(1-tmp)*(1-piD)
      n_com[2,9]<-tmp*n_obs[2]/(1-tmp)*piD
      for (i in 2:8) {
        n_com[2,i]<-piD*sim_data[2,i]
        n_com[2,i+8]<-sim_data[2,i]-n_com[2,i]
      }
      
      N_tmp<-sum(n_com)
      if(loop==1){
        N_list0[s]<-N_tmp
      }
      
      sumA<-sumB<-sumC<-sumD<-0
      for (jA in 0:1) {
        for (jB in 0:1) {
          for (jC in 0:1) {
            for (jD in 0:1) {
              ind<-1+jA+2*jB+4*jC+8*jD
              tmp<-sum(n_com[,ind])
              sumA=sumA+jA*tmp
              sumB=sumB+jB*tmp
              sumC=sumC+jC*tmp
              sumD=sumD+jD*tmp
            }
          }
        }
      }
      
      piA<-sumA/N_tmp
      piB<-sumB/N_tmp
      piC<-sumC/N_tmp
      piD<-sumD/N_tmp
    }
    if(loop==1){
      plot(N_list0,type = "l",xlab = "steps",ylab = "N")
    }
    N_list[t]<-round(N_tmp)
  }
  return(N_list)
}