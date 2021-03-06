BiMOL.density=function(Xs,Ys,s,t,xlims,ylims,N ,delt ,mu1,mu2,sig11,sig12,sig21,sig22,final.only=FALSE,show.lattice =FALSE)
{
   if((missing(mu1)&&missing(mu2)&&missing(sig11)&&missing(sig12)&&missing(sig21)&&missing(sig22)))
  {

    namess=c('mu1','mu2','sig11','sig12','sig21','sig22')
    namess2=c('muu1','muu2','sigg11','sigg12','sigg21','sigg22')
    txt =rep('+matrix(0,dim(X)[1],dim(Y)[2])',6)
    func.list=rep(0,length(namess))
    obs=objects(pos=1)
    muu1   =function(X,Y,t){}
    muu2   =function(X,Y,t){}
    sigg11 =function(X,Y,t){}
    sigg12 =function(X,Y,t){}
    sigg21 =function(X,Y,t){}
    sigg22 =function(X,Y,t){}
    for(i in 1:length(namess))
    {
      if(sum(obs==namess[i]))
      {
        txt[i] = paste0(body(namess[i])[2],txt[i])
      }
    }
    body(muu1)    =  parse(text=txt[1])
    body(muu2)    =  parse(text=txt[2])
    body(sigg11)  =  parse(text=txt[3])
    body(sigg12)  =  parse(text=txt[4])
    body(sigg21)  =  parse(text=txt[5])
    body(sigg22)  =  parse(text=txt[6])
   }else
   {
    muu1   =function(X,Y,t){}
    muu2   =function(X,Y,t){}
    sigg11 =function(X,Y,t){}
    sigg12 =function(X,Y,t){}
    sigg21 =function(X,Y,t){}
    sigg22 =function(X,Y,t){}
     body(muu1) =  parse(text=paste0(mu1,'+matrix(0,dim(X)[1],dim(Y)[2])'))
     body(muu1) =  parse(text=paste0(mu1,'+matrix(0,dim(X)[1],dim(Y)[2])'))
     body(sigg11) =  parse(text=paste0(sig11,'+matrix(0,dim(X)[1],dim(Y)[2])'))
     if(!missing(sigg12)){body(sigg12) =  parse(text=paste0(sig12,'+matrix(0,dim(X)[1],dim(Y)[2])'))}
     if(!missing(sigg21)){body(sigg21) =  parse(text=paste0(sig21,'+matrix(0,dim(X)[1],dim(Y)[2])'))}
     body(sigg22) =  parse(text=paste0(sig22,'+matrix(0,dim(X)[1],dim(Y)[2])'))
   }

  f=function(U,tme)
  {
    MU1= muu1(XX1,XX2,tme)*U
    MU2= muu2(XX1,XX2,tme)*U
    SU1= sigg11(XX1,XX2,tme)^2*U
    SU2= sigg22(XX1,XX2,tme)^2*U
    D1 = (MU1[-c(1,2),-c(1,N)]-MU1[-c(N-1,N),-c(1,N)])/dx1/2
    D2 = (MU2[-c(1,N),-c(1,2)]-MU2[-c(1,N),-c(N-1,N)])/dx2/2
    D3 =  1/2*((SU1[-c(1,N),-c(1,2)]-2*SU1[-c(1,N),-c(1,N)]+SU1[-c(1,N),-c(N-1,N)]))/dxx1
    D4 =  1/2*((SU2[-c(1,2),-c(1,N)]-2*SU2[-c(1,N),-c(1,N)]+SU2[-c(N-1,N),-c(1,N)]))/dxx2
    MMM1=-D1-D2+D3+D4
    return(cbind(0,rbind(0,MMM1,0),0))
  }
  
  
  xx1=seq(xlims[1],xlims[2],length=N)
  xx2=seq(ylims[1],ylims[2],length=N)
  
  wwwhhh1=min(which(abs(xx1-Xs)==min(abs(xx1-Xs))))
  wwwhhh2=min(which(abs(xx2-Ys)==min(abs(xx2-Ys))))
  
  xx1=xx1+(Xs-xx1[wwwhhh1])
  xx2=xx2+(Ys-xx2[wwwhhh2])
  
  dx1=diff(xx1)[1]
  dxx1=dx1^2
  dx2=diff(xx2)[1]
  dxx2=dx2^2
  XX1=outer(xx1,rep(1,N))
  XX2=t(outer(xx2,rep(1,N)))
  

  
  if(show.lattice)
  {
    plot(1,1,type='n',xlab='',ylab='',xlim=xlims,ylim=ylims,main='Lattice')
    #abline(h=xx2,lty='dotted')
    #abline(v=xx1,lty='dotted')
    xlims
    
    xlims=range(xx1)
    ylims=range(xx2)
    segments(xlims[1],ylims[1],xlims[1],ylims[2],lty='solid')
    segments(xlims[2],ylims[1],xlims[2],ylims[2],lty='solid')
    segments(xlims[1],ylims[2],xlims[2],ylims[2],lty='solid')
    segments(xlims[1],ylims[1],xlims[2],ylims[1],lty='solid')
    for(i in 1:length(xx2))
    {
      points(rep(xx2[i],length(xx1))~xx1,pch=19,cex=0.5,col='grey50')
    }
    points(xx2[wwwhhh2]~xx1[wwwhhh1],col='red',pch=19)
    axis(2,at=xx2,tcl=-0.25,labels=NA)
    axis(1,at=xx1,tcl=-0.25,labels=NA)
  }
  M=matrix(0,N,N)
  M[wwwhhh1,wwwhhh2]=1
  #print(c(xx1[wwwhhh1],xx2[wwwhhh2]))
  
  N.mesh=round((t-s)/delt)+1
  if(!final.only)
  {
    MM = array(0,dim=c(N,N,N.mesh))
  }
  ttt=s
  
  t.alpha=
    c(0.000000000000000000000000000000000000000000000000000000000000,
      0.100000000000000000000000000000000000000000000000000000000000,
      0.539357840802981787532485197881302436857273449701009015505500,
      0.809036761204472681298727796821953655285910174551513523258250,
      0.309036761204472681298727796821953655285910174551513523258250,
      0.981074190219795268254879548310562080489056746118724882027805,
      0.833333333333333333333333333333333333333333333333333333333333,
      0.354017365856802376329264185948796742115824053807373968324184,
      0.882527661964732346425501486979669075182867844268052119663791,
      0.642615758240322548157075497020439535959501736363212695909875,
      0.357384241759677451842924502979560464040498263636787304090125,
      0.117472338035267653574498513020330924817132155731947880336209,
      0.833333333333333333333333333333333333333333333333333333333333,
      0.309036761204472681298727796821953655285910174551513523258250,
      0.539357840802981787532485197881302436857273449701009015505500,
      0.100000000000000000000000000000000000000000000000000000000000,
      1.00000000000000000000000000000000000000000000000000000000000)
  
  
  for(i in 2:N.mesh)
  {
    
    x0=M
    fx0=f(x0,ttt)
    x1=x0+delt*(0.1*fx0)
    fx1=f(x1,ttt+t.alpha[2]*delt)
    x2=x0+delt*(-0.915176561375291*fx0  +1.45453440217827*fx1)
    fx2=f(x2,ttt+t.alpha[3]*delt)
    x3=x0+delt*( 0.202259190301118*fx0  +0.606777570903354*fx2)
    fx3=f(x3,ttt+t.alpha[4]*delt)
    x4=x0+delt*( 0.184024714708644*fx0  +0.197966831227192*fx2-0.0729547847313633*fx3)
    fx4=f(x4,ttt+t.alpha[5]*delt)
    x5=x0+delt*( 0.0879007340206681*fx0 +0.410459702520261*fx3+0.482713753678866*fx4)
    fx5=f(x5,ttt+t.alpha[6]*delt)
    x6=x0+delt*(0.085970050490246*fx0   +0.330885963040722*fx3+0.48966295730945*fx4-0.0731856375070851*fx5)
    fx6=f(x6,ttt+t.alpha[7]*delt)
    x7=x0+delt*(0.120930449125334*fx0   +0.260124675758296*fx4+0.0325402621549091*fx5-0.0595780211817361*fx6)
    fx7=f(x7,ttt+t.alpha[8]*delt)
    x8=x0+delt*(0.110854379580391*fx0   -0.0605761488255006*fx5+0.321763705601778*fx6+0.510485725608063*fx7)
    fx8=f(x8,ttt+t.alpha[9]*delt)
    x9=x0+delt*(0.112054414752879*fx0   -0.144942775902866*fx5-0.333269719096257*fx6+0.49926922955688*fx7+0.509504608929686*fx8)
    fx9=f(x9,ttt+t.alpha[10]*delt)
    x10=x0+delt*(0.113976783964186*fx0  -0.0768813364203357*fx5+0.239527360324391*fx6+0.397774662368095*fx7+0.0107558956873607*fx8-0.327769124164019*fx9)
    fx10=f(x10,ttt+t.alpha[11]*delt)
    x11=x0+delt*(0.0798314528280196*fx0 -0.0520329686800603*fx5-0.0576954146168549*fx6+0.194781915712104*fx7+0.145384923188325*fx8-0.0782942710351671*fx9-0.114503299361099*fx10)
    fx11=f(x11,ttt+t.alpha[12]*delt)
    x12=x0+delt*(0.985115610164857*fx0  +0.330885963040722*fx3+0.48966295730945*fx4-1.37896486574844*fx5-0.861164195027636*fx6+5.78428813637537*fx7+3.28807761985104*fx8-2.38633905093136*fx9-3.25479342483644*fx10-2.16343541686423*fx11)
    fx12=f(x12,ttt+t.alpha[13]*delt)
    x13=x0+delt*(0.895080295771633*fx0  +0.197966831227192*fx2-0.0729547847313633*fx3-0.851236239662008*fx5+0.398320112318533*fx6+3.63937263181036*fx7+1.5482287703983*fx8-2.12221714704054*fx9-1.58350398545326*fx10-1.71561608285936*fx11-0.0244036405750127*fx12)
    fx13=f(x13,ttt+t.alpha[14]*delt)
    x14=x0+delt*(-0.915176561375291*fx0+1.45453440217827*fx1+0*fx2+0*fx3-0.777333643644968*fx4+0*fx5-0.0910895662155176*fx6+0.0910895662155176*fx12+0.777333643644968*fx13)
    fx14=f(x14,ttt+t.alpha[15]*delt)
    x15=x0+delt*(0.1*fx0-0.157178665799771*fx2+0.157178665799771*fx14)
    fx15=f(x15,ttt+t.alpha[16]*delt)
    x16=x0+delt*(0.181781300700095*fx0+0.675*fx1+0.34275815984719*fx2+0*fx3+0.259111214548323*fx4-0.358278966717952*fx5-1.04594895940883*fx6+0.930327845415627*fx7+1.77950959431708*fx8+0.1*fx9-0.282547569539044*fx10-0.159327350119973*fx11-0.145515894647002*fx12-0.259111214548323*fx13-0.34275815984719*fx14-0.675*fx15)
    fx16=f(x16,ttt+t.alpha[17]*delt)
    
    M=M+(0.0333333333333333333333333333333333333333333333333333333333333*fx0
         +0.0250000000000000000000000000000000000000000000000000000000000*fx1
         +0.0333333333333333333333333333333333333333333333333333333333333*fx2
         +0.000000000000000000000000000000000000000000000000000000000000*fx3
         +0.0500000000000000000000000000000000000000000000000000000000000*fx4
         +0.000000000000000000000000000000000000000000000000000000000000*fx5
         +0.0400000000000000000000000000000000000000000000000000000000000*fx6
         +0.000000000000000000000000000000000000000000000000000000000000*fx7
         +0.189237478148923490158306404106012326238162346948625830327194*fx8
         +0.277429188517743176508360262560654340428504319718040836339472*fx9
         +0.277429188517743176508360262560654340428504319718040836339472*fx10
         +0.189237478148923490158306404106012326238162346948625830327194*fx11
         -0.0400000000000000000000000000000000000000000000000000000000000*fx12
         -0.0500000000000000000000000000000000000000000000000000000000000*fx13
         -0.0333333333333333333333333333333333333333333333333333333333333*fx14
         -0.0250000000000000000000000000000000000000000000000000000000000*fx15
         +0.0333333333333333333333333333333333333333333333333333333333333*fx16)*delt
    ttt=ttt+delt
    if(!final.only)
    {
      MM[,,i]=M
    }
    if(any(is.na(M))){break; warning(paste("NAs produced at time",ttt,". Perhaps increase NN and decrease delt. Alternatively, widen xlims."))}
  }
  if(final.only)
  {
     MM=M
  }
  ret =(list(density=MM/dx1/dx2,Xt=xx1,Yt=xx2,time=seq(s,t,delt)))
  class(ret) = 'BiMOL.density'
  return(ret)
}