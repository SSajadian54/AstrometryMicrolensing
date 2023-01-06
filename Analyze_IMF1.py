import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
rcParams["font.size"] = 13
rcParams["font.family"] = "sans-serif"
rcParams["font.sans-serif"] = ["Computer Modern Sans"]
rcParams["text.usetex"] = True
rcParams["text.latex.preamble"] = r"\usepackage{cmbright}"
###############################################################
yr=float(365.2421875)
kp=3.08568025*pow(10.,19);
au=1.4960*pow(10.0,11.0);
pp=float(100.0)
###############################################################

f1=open("./files/MONT1/IBH1.txt","r")
f2=open("./files/MONT1/IBH2.txt","r")
f3=open("./files/MONT1/IBH3.txt","r")
nm= sum(1 for line in f1) 

par=np.zeros((nm,49)) 
dat=np.zeros((nm,36)) 
mat=np.zeros((nm,10)) 

par= np.loadtxt("./files/MONT1/IBH1.txt")
dat= np.loadtxt("./files/MONT1/IBH2.txt")
mat= np.loadtxt("./files/MONT1/IBH3.txt")

fij=open("./Histo_BHs1/results1.txt","w")
fij.close(); 
fij=open("./Histo_BHs1/results2.txt","w")
fij.close(); 

arry2=np.zeros((151746, 16))
arry2=np.loadtxt("/home/sajadian/0_ARCHIVED/10_disk/Histo_BHs1/without_Sparse/results1.txt")

##############################################################
nd=7;  na=3;  nq=16; nv=11
F = np.zeros([nd , nd]);  
G = np.zeros([na , na])
Era=np.zeros((nd)); 
Erb=np.zeros((na));
resu= np.zeros((nq))
numc= np.zeros(nq)
numa= np.zeros(nq)
numb= np.zeros(nq)
hisd0= np.zeros((nm,10))
hisd1= np.zeros((nm,10))
epsi= np.zeros(3)
ngood=0.0;
arry=np.zeros((nm, nq+4))
l=0;  nc=0; 
xar=np.zeros((nv,4))
#xar[:,0]=np.array([2.0, 6.8 , 11.6, 16.4, 21.2, 26.0, 30.8, 35.6, 40.4, 42.2, 49.0])
xar[:,0]=np.array([2.0, 4.8 , 8.0, 14.0, 20.0, 28.0, 33.0, 38.0, 42.5, 46.5, 50.0])
xar[:,1]=np.array([ 0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])
xar[:,2]=np.array([0.1,2.0,4.5,6.0,7.0, 7.7, 8.5, 9.2, 10.0, 11.0, 12.5])
xar[:,3]=np.array([15.5,16.4,17.3,18.2,19.1,20.0,20.9,21.8, 22.7, 23.6, 24.5])
##############################################################
def func(para,cc): 
    i1=-1
    for j in range(nv-1):  
        if( float((para-xar[j,cc])*(para-xar[j+1,cc]))<0.0  or  para==xar[j,cc]):
            i1=j;  
            break
    return(i1)        
##############################################################
## plotting the histograms of all parameters
for i in range(49):
    plt.clf()
    fig= plt.figure(figsize=(8,6))
    ax= plt.gca()              
    plt.hist(par[:,i],30,histtype='bar',ec='darkgreen',facecolor='green',alpha=0.7, rwidth=1.5)
    y_vals = ax.get_yticks()
    ax.set_yticklabels(['{:.2f}'.format(1.0*x*(1.0/nm)) for x in y_vals]) 
    y_vals = ax.get_yticks()
    plt.ylim([np.min(y_vals), np.max(y_vals)])
    ax.set_ylabel(r"$\rm{Normalized}~\rm{Distribution}$",fontsize=19,labelpad=0.1)
    plt.xticks(fontsize=17, rotation=0)
    plt.yticks(fontsize=17, rotation=0)
    plt.grid("True")
    plt.grid(linestyle='dashed')
    fig=plt.gcf()
    fig.savefig("./Histo_BHs1/hist_{0:d}.jpg".format(i),dpi=200)
print ("****  All histos are plotted *****************************" )   
##############################################################
'''
plt.clf()
fig= plt.figure(figsize=(8,6))
ax= plt.gca()              
plt.plot(par[:,10],par[:,29], "ro", ms=0.2)
plt.xticks(fontsize=17, rotation=0)
plt.yticks(fontsize=17, rotation=0)
plt.grid("True")
plt.grid(linestyle='dashed')
fig=plt.gcf()
fig.savefig("./Histo_BHs1/Vtxls.jpg",dpi=200)
'''
#############################################################
## reading all raws  in the files 
mapp=np.zeros((nv,nv,2)); 
for i in range(nm): 
    coun1, lat, lon= par[i,0], par[i,1], par[i,2]
    strucl, Ml, Dl, vl= par[i,3], par[i,4], par[i,5], par[i,6]
    strucs, cl, mass, Ds,Tstar, Rstar, logl= par[i,7], par[i,8], par[i,9], par[i,10], par[i,11], par[i,12], par[i,13]
    types, col, vs, MI, MW149, mI, mW149=    par[i,14], par[i,15],par[i,16], par[i,17], par[i,18], par[i,19], par[i,20]
    magbI, mbs, blendI, fb, NbI, Nbw, ExI, ExW= par[i,21], par[i,22],par[i,23], par[i,24], par[i,25], par[i,26], par[i,27], par[i,28]
    tE, RE, t0, mul, Vt, u0, opd,ros,tetE=par[i,29],par[i,30],par[i,31], par[i,32], par[i,33], par[i,34], par[i,35],par[i,36], par[i,37]
    flagf,flagD, dchi, ndw,li, mus1, mus2=par[i,38], par[i,39],par[i,40], par[i,41], par[i,42], par[i,43], par[i,44]
    xi, mul1, mul2, piE=  par[i,45], par[i,46],par[i,47], par[i,48]
    mus=np.sqrt(mus1*mus1+mus2*mus2)
    mul=np.sqrt(mul1*mul1+mul2*mul2)
    coun2,t0b, u0b, tEb, xib, fbb, mbsb, piEb= dat[i,0], dat[i,1], dat[i,2], dat[i,3], dat[i,4], dat[i,5], dat[i,6], dat[i,7];
    F[0,0]=dat[i,8]
    F[1,0], F[1,1]=dat[i,9], dat[i,10];
    F[2,0], F[2,1], F[2,2]= dat[i,11], dat[i,12], dat[i,13];
    F[3,0], F[3,1], F[3,2], F[3,3]=dat[i,14], dat[i,15], dat[i,16], dat[i,17];
    F[4,0], F[4,1], F[4,2], F[4,3], F[4,4]= dat[i,18], dat[i,19], dat[i,20], dat[i,21], dat[i,22];
    F[5,0], F[5,1], F[5,2], F[5,3], F[5,4], F[5,5]=dat[i,23], dat[i,24], dat[i,25], dat[i,26], dat[i,27], dat[i,28];
    F[6,0], F[6,1], F[6,2], F[6,3], F[6,4], F[6,5], F[6,6]= dat[i,29], dat[i,30], dat[i,31], dat[i,32], dat[i,33], dat[i,34], dat[i,35]   
    coun3, tetEb, mus1b, mus2b = mat[i,0], mat[i,1], mat[i,2], mat[i,3]
    G[0][0]= mat[i,4]
    G[1][0], G[1][1]= mat[i,5], mat[i,6]
    G[2][0], G[2][1], G[2][2]= mat[i,7], mat[i,8], mat[i,9]
    
    if(coun1!=coun2 or coun2!=coun3 or coun3!=coun1 or abs(t0b-t0)>0.001 or abs(u0b-u0)>0.001 or abs(tEb-tE)>0.001 or abs(xib-xi)>0.001 
    or abs(fb-fbb)>0.001 or abs(mbs-mbsb)>0.001 or abs(piEb-piE)>0.001 or abs(tetEb-tetE)>0.001 or abs(mus1-mus1b)>0.001 or 
    abs(mus2-mus2b)>0.001 or Dl>Ds or u0>1.0 or t0>(5.0*yr) or tE<0.0 or Ml<2.0  or Ml>50.0 or flagD<1 
    or dchi<800.0  or tetE<0.0  or vl<0.0  or Vt<0.0  or mass<0.0 or Nbw<1.0 ): 
        print ("Error:  ",  coun1, coun2, coun3, t0, t0b,  u0, u0b,  tE, tEb,  xi, xib, fb, fbb, mbs, mbsb)
        print ("Error ",  piE, piEb, tetE, tetEb, mus1, mus1b, mus2, mus2b,  Dl, Ds, Ml )
        print ("Error ",  flagD, dchi, vl, vs, Vt, mass, Nbw )
        input("Enter a number !!!")
    
    
    for k in range(nd): 
        F[k,(k+1):]= F[(k+1):,k]
    for k in range(na): 
        G[k,(k+1):]= G[(k+1):,k]
       
    
    corr1=0.0;
    inverF= np.mat(F).I
    inverG= np.mat(G).I
    for k in range(nd):  
        Era[k]=np.sqrt(abs(inverF[k,k]))  
    corr1= inverF[2,3]/(Era[2]*Era[3])## correlation between tE, xi         
    for k in range(na):  
        Erb[k]=np.sqrt(abs(inverG[k,k]))     
    if(all(Era)>0.0 and all(Erb)>0.0):    
        resu[0]=abs(Era[0]*1.0/t0)  
        resu[1]=abs(Era[1]*1.0/u0)
        resu[2]=abs(Era[2]*1.0/tE)  
        resu[3]=abs(Era[3]*1.0/xi)
        resu[4]=abs(Era[4]*1.0/fb)   
        resu[5]=abs(Era[5]*1.0/mbs)  
        resu[6]=abs(Era[6]*1.0/piE)
        resu[7]=abs(Erb[0]*1.0/tetE)
        resu[8]=abs(Erb[1]*1.0/mus1)
        resu[9]=abs(Erb[2]*1.0/mus2)
        resu[10]=np.sqrt((resu[7])**2.0 + (resu[6])**2.0 )### relative error in mass 
        resu[11]=abs(resu[10]*(Ds-Dl)/Ds)## relative error in Dl 
        f1= abs(resu[7]**2.0 + resu[2]**2.0 + (Era[3]*np.tan(xi))**2.0 -2.0*resu[2]*Era[3]*np.tan(xi)*corr1 )
        f2= abs(resu[7]**2.0 + resu[2]**2.0 + (Era[3]/np.tan(xi))**2.0 -2.0*resu[2]*Era[3]/np.tan(xi)*corr1 )
        sig_mul1=np.sqrt( Erb[1]**2.0 + mul**2.0 *(np.cos(xi))**2.0 *f1)
        sig_mul2=np.sqrt( Erb[2]**2.0 + mul**2.0 *(np.sin(xi))**2.0 *f2)  
        resu[12]= abs(sig_mul1)/abs(mul1)
        resu[13]= abs(sig_mul2)/abs(mul2)
        #resu[14]=np.sqrt(Erb[1]**2.0 + Erb[2]**2.0+mul**2.0*(resu[7]**2.0+resu[2]**2.0+Era[3]**2.0))/mul
        resu[14]=np.sqrt(sig_mul1*sig_mul1 + sig_mul2*sig_mul2)/mul#relative error in Lens proper motion
        resu[15]=np.sqrt(Erb[1]**2.0 + Erb[2]**2.0)/mus#            relative error in source proper motion
        resu[:]=resu[:]*pp### in persent 
        '''
        print ("******************** counter:  ", i )
        print ("Err(t0)*100.0/t0:  ",     resu[0] )
        print ("Err(u0)*100.0/u0:  ",     resu[1] ) 
        print ("Err(tE)*100.0/tE:  ",     resu[2]  )
        print ("Err(ksi)*100.0/ksi:  ",   resu[3] )
        print ("Err(fb)*100.0/fb:  ",     resu[4]   )          
        print ("Err(mbs)*100.0/mbs:  ",   resu[5] )
        print ("Err(piE)*100.0/piE:  ",   resu[6]   )
        print ("Err(tetE)*100.0/tetE:  ", resu[7] )
        print ("Err(mus1)*100.0/mus1:  ", resu[8])
        print ("Err(mus2)*100.0/mus2:  ", resu[9])
        print ("Err(Mass)*100.0/Mass: " , resu[10]) 
        print ("Err(Dl)*100.0/Dl: " ,     resu[11] )
        print ("Err(mul1)*100.0/mul1: " , resu[12] )
        print ("Err(mul2)*100.0/mul2: " , resu[13] ) 
        print ("Err(mul)*100.0/mul : " ,  resu[14] ) 
        print ("Err(mus)*100.0/mus : " ,  resu[15] ) 
        '''
        fij=open("./Histo_BHs1/results1.txt","a")
        np.savetxt(fij,resu.reshape((1,nq)),fmt="%.7e %.7e  %.7e  %.7e  %.7e  %.7e  %.7e  %.7e %.7e  %.7e  %.7e   %.7e   %.7e  %.7e  %.7e  %.7e") 
        fij.close()   
      
         
        if(np.isfinite(resu).all and all(resu)>0.0):     
            for j in range(nq): 
                if(resu[j]<10.0):  numc[j]+= 1.0
                if(resu[j]< 5.0):  numb[j]+= 1.0
                if(resu[j]< 1.0):  numa[j]+= 1.0
                if(resu[j]==0.0):  resu[j]= +1.0e-10
            
            arry[l,:nq]=resu
            arry[l, nq]=   func(Ml,0)
            arry[l, nq+1]= func(float(Dl/Ds),1)
            arry[l, nq+2]= func(Ds,2)
            arry[l, nq+3]= func(mbs,3)
            if(resu[10]<5.0 and resu[11]<5.0): 
                mapp[int(arry[l,nq+2]) , int(arry[l,nq]) , 0 ] += Dl/Ds #abs(resu[10]) + abs(resu[11])
                mapp[int(arry[l,nq+2]) , int(arry[l,nq]) , 1 ] += 1.0
                
            hisd0[l,:]= np.array([np.log10(tE), mbs, np.log10(piE), np.log10(tetE), t0/yr,  np.log10(mus), np.log10(mul), u0, Dl, Ds ])
            l+=1
            if(resu[10]<1.0 and resu[11]<1.0 and resu[14]<1.0):     epsi[0]+=1.0
            if(resu[10]<10.0 and resu[11]<10.0 and resu[14]<10.0):  epsi[2]+=1.0
            if(resu[10]<5.0 and resu[11]<5.0 and resu[14]<5.0):     
                epsi[1]+=1.0 
                hisd1[nc,:]= np.array([np.log10(tE), mbs, np.log10(piE), np.log10(tetE), t0/yr,  np.log10(mus), np.log10(mul), u0, Dl, Ds ])
                nc+=1
print ("*****************************************************" )
print ("tot_number:  ", nm, " accepted number:  ",  l, "detected number:  ",  nc) 
print ("epsilon:  ",   epsi[0],       epsi[1],    epsi[2]) 
epsi= epsi*(100.0/l)       
print ("epsilon[\%]:  ",  epsi[0],    epsi[1],    epsi[2] )
print ("*****************************************************" )
######################################################################
xdet=[r"$\log_{10}[t_{\rm E}(days)]$",  r"$m_{\rm{base}}(\rm{mag})$", r"$\log_{10}[\pi_{\rm{E}}]$",  r"$\log_{10}[\theta_{\rm{E}}(\rm{mas})]$",  r"$t_{0}(years)$", r"$\log_{10}[\mu_{s}(\rm{mas/days})]$", r"$\log_{10}[\mu_{\rm{l}}(\rm{mas/days})]$",  r"$u_{0}$",  r"$D_{\rm{l}}(\rm{kpc})$",  r"$D_{\rm{s}}(\rm{kpc})$"]

#xdet=[r"$\left<t_{\rm E}(days)\right>=$",  r"$\left<m_{\rm{base}}(\rm{mag})$", r"$\log_{10}[\pi_{\rm{E}}]$",  r"$\log_{10}[\theta_{\rm{E}}(\rm{mas})]$",  r"$t_{0}(years)$", r"$\log_{10}[\mu_{s}(\rm{mas/days})]$", r"$\log_{10}[\mu_{\rm{l}}(\rm{mas/days})]$",  r"$u_{0}$",  r"$D_{\rm{l}}(\rm{kpc})$",  r"$D_{\rm{s}}(\rm{kpc})$"]

### histo of events with measured errors and detected events
xd1=[1.1,15.5, -2.8,-0.7, 0.0, -3.0,-3.5, 0.0,  0.0,   2.0]
xd2=[3.8,24.0, -0.6, 1.5, 5.0, -1.35,-1.0, 1.0, 10.0, 13.0]

lab1=[np.mean(np.power(10,hisd0[:l,0])), np.mean(hisd0[:l,1]), np.mean(np.power(10,hisd0[:l,2])), 
      np.mean(np.power(10,hisd0[:l,3])), np.mean(hisd0[:l,4]), np.mean(np.power(10,hisd0[:l,5])), 
      np.mean(np.power(10,hisd0[:l,6])), np.mean(hisd0[:l,7]), np.mean(hisd0[:l,8]), np.mean(hisd0[:l,9])]

lab2=[np.mean(np.power(10,hisd1[:nc,0])), np.mean(hisd1[:nc,1]), np.mean(np.power(10,hisd1[:nc,2])), 
      np.mean(np.power(10,hisd1[:nc,3])), np.mean(hisd1[:nc,4]), np.mean(np.power(10,hisd1[:nc,5])), 
      np.mean(np.power(10,hisd1[:nc,6])), np.mean(hisd1[:nc,7]), np.mean(hisd1[:nc,8]), np.mean(hisd1[:nc,9])]

for i in range(10):
    plt.clf()
    plt.clf()
    fig= plt.figure(figsize=(8,6))
    ax1= plt.gca()
    plt.hist(hisd0[:l,i], 26,histtype='bar',ec='darkgreen',facecolor='green',alpha=0.55,rwidth=1.5,label=str(round(lab1[i],2))  )
    plt.hist(hisd1[:nc,i],26,histtype='step',color='k',alpha=0.9, lw=2.2,label= str(round(lab2[i],2))   )
    y_vals=ax1.get_yticks()
    ax1.set_yticklabels(['{:.2f}'.format(1.0*x*(1.0/len(hisd0[:l,i]))) for x in y_vals]) 
    y_vals = ax1.get_yticks()
    ax1.set_ylim([np.min(y_vals)+0.01, np.max(y_vals)-0.01])
    plt.yticks(fontsize=18, rotation=0)
    plt.xticks(fontsize=18, rotation=0)
    ax1.set_ylabel(r"$\rm{Normalized}~\rm{Distribution}$",fontsize=19,labelpad=0.1)
    ax1.set_xlabel(str(xdet[i]) , fontsize=19,  labelpad=0.1)
    plt.xlim([ xd1[i],  xd2[i] ])
    #plt.axvline(x=np.mean(hisd0[:l ,i]) , color='darkgreen', linestyle='--', lw=2.3)
    #plt.axvline(x=np.mean(hisd1[:nc,i]) , color='k',        linestyle='--', lw=2.3)
    plt.grid("True")
    plt.grid(linestyle='dashed')
    plt.legend()
    plt.legend(loc='best',fancybox=True, shadow=True)
    plt.legend(prop={"size":18})
    fig= plt.gcf()
    fig.tight_layout()
    fig.savefig("./Histo_BHs1/histdet_{0:d}.jpg".format(i),dpi=200)
print ("****  All histos are plotted *****************************" )   

######################################################################        
#pi1=np.zeros((nv,3))        
m1=np.zeros((nv,4,4))
for i in range(nv): 
    for k in range(4):
        arr2=np.zeros((l,4));
        k1=0; k2=0; k3=0; k4=0; ##k5=0; k6=0; 
        for j in range(l): 
            if(arry[j,nq+k]==i):
                if(arry[j,10]<5.0):   arr2[k1,0]=arry[j,10]; k1+=1 ##mass
                if(arry[j,11]<5.0):   arr2[k2,1]=arry[j,11]; k2+=1 ##Dl
                if(arry[j,14]<5.0):   arr2[k3,2]=arry[j,14]; k3+=1 ##mul
                if(arry[j,6]< 5.0):   arr2[k4,3]=arry[j, 6]; k4+=1 ##piE
                #if(k==0 and arry[j,6]<5.0 ):  arr2[k4,3]=arry[j,6]; k4+=1
                #if(k==0 and arry[j,2]<5.0 ):  arr2[k5,4]=arry[j,2]; k5+=1
                #if(k==0 and arry[j,3]<5.0 ):  arr2[k6,5]=arry[j,3]; k6+=1
        m1[i,0,k]=np.mean(arr2[:k1,0]); 
        m1[i,1,k]=np.mean(arr2[:k2,1]);  
        m1[i,2,k]=np.mean(arr2[:k3,2]); 
        m1[i,3,k]=np.mean(arr2[:k4,3]); 
        #if(k==0): 
        #    pi1[i,0]= np.mean(arr2[:k4,3]);
        #    pi1[i,1]= np.mean(arr2[:k5,4]);
        #    pi1[i,2]= np.mean(arr2[:k6,5]);
        print ("<Err_mass>,  <Err_Dl>,  <Err_mul> :   ",  m1[i,0,k],   m1[i,1,k],   m1[i,2,k],  k,  i)

###############################################################      
for j in range(nv-1): 
    print("m_Dls:  ",   np.sum(mapp[:,j,0])/np.sum(mapp[:,j,1]) )

mapp[:,:,0]= mapp[:,:,0]/mapp[:,:,1]  
cmap=plt.get_cmap('viridis')
plt.clf()
plt.clf()
plt.figure()
fig= plt.figure(figsize=(8,8))
plot=plt.imshow(mapp[:,:,0] ,interpolation='nearest',extent=(2.0,45.0, 0.0, 12.5), origin='lower', cmap=cmap, aspect='auto')
ax=plt.gca()
ax.set_xticks(xar[:(nv-1),0] )
ax.set_yticks(xar[:(nv-1),2] )
ax.tick_params(direction='out',pad=5,top=False, right=False)
plt.tight_layout()
cb=plt.colorbar(plot)
cb.set_label(r"$\rm{Relative}~\rm{Error}$", fontsize=17, labelpad=0.0)
plt.xticks(fontsize=17, rotation=0)
plt.yticks(fontsize=17, rotation=0)
plt.xlabel(r"$M_{\rm l}(M_{\odot})$", fontsize=17)
plt.ylabel(r"$x_{\rm ls}$",  fontsize=17)
fig=plt.gcf()
fig.savefig("mapp2D.jpg", dpi=200) 
################################################################### 
nam=[r"$\log_{10}[\sigma_{t_{0}}\big/t_{0}(\%) ]$", 
     r"$\log_{10}[\sigma_{u_{0}}\big/u_{0}(\%)]$", 
     r"$\log_{10}[\sigma_{t_{\rm E}}\big/t_{\rm E}(\%)]$", 
     r"$\log_{10}[\sigma_{\xi}\big/ \xi(\%)]$",     
     r"$\log_{10}[\sigma_{f_{b}}\big/f_{b}(\%)]$", 
     r"$\log_{10}[\sigma_{m_{\rm base}}\big/m_{\rm base}(\%)]$", 
     r"$\log_{10}[\sigma_{\pi_{\rm E}}\big/ \pi_{\rm E}(\%)]$",  
     r"$\log_{10}[\sigma_{\theta_{\rm E}}\big/ \theta_{\rm E}(\%)]$",  
     r"$\log_{10}[\sigma_{\mu_{\rm s, n1}}\big/\mu_{\rm s, n1}(\%)]$", 
     r"$\log_{10}[\sigma_{\mu_{\rm s, n2}}\big/ \mu_{\rm s, n2}(\%)]$", 
     r"$\log_{10}[\sigma_{M_{\rm l}}\big/ M_{\rm l}(\%)]$", 
     r"$\log_{10}[\sigma_{D_{\rm l}}\big/ D_{\rm l}(\%)]$", 
     r"$\log_{10}[\sigma_{\mu_{\rm l, n1}}\big/ \mu_{\rm l, n1}(\%)]$", 
     r"$\log_{10}[\sigma_{\mu_{\rm l, n2}}\big/ \mu_{\rm l, n2}(\%)]$",
     r"$\log_{10}[\sigma_{\mu_{\rm l}}\big/ \mu_{\rm l}(\%)]$", 
     r"$\log_{10}[\sigma_{\mu_{\rm s}}\big/ \mu_{\rm s}(\%)]$"]        
################################################################### 

x1=[-5.0,-2.0,-2.0, -3.0,-2.0, -5.0, -2.0,-2.5, -3.0, -3.0, -1.5, -2.5, -3.0 ,  -2.0, -1.5,-4.5]
x2=[ 3.0, 6.0, 3.5,  5.0, 4.5,  1.3,  5.0, 1.3,  1.3,  1.3,  5.0 , 3.5 , 6.0,    5.9,  5.0, 1.3]    
for i in range(nq): 
    plt.clf()
    fig= plt.figure(figsize=(8,6))
    ax= plt.gca()  
    ll=np.min((len(arry2[:,0]),l))         
    #plt.hist(np.log10(arry2[:,i]),33,histtype='bar',ec='darkred',facecolor='red',   alpha=0.35, rwidth=1.5)
    
    plt.hist(np.log10(arry[:ll,i]),35, histtype='bar',ec='darkgreen',facecolor='green',alpha=0.55, rwidth=1.5)
    #y_vals = ax.get_yticks()
    #ax.set_yticklabels(['{:.2f}'.format(1.0*x*(1.0/l)) for x in y_vals]) 
    
    
    plt.hist(np.log10(arry2[:ll,i]),35, histtype='step',color='k',alpha=0.9, lw=2.2  )
    y_vals = ax.get_yticks()
    ax.set_yticklabels(['{:.2f}'.format(1.0*x*(1.0/ll)) for x in y_vals]) 
    
    y_vals = ax.get_yticks()
    plt.ylim([np.min(y_vals), np.max(y_vals)])
    
    #plt.axvline(x= np.log10(np.mean(arry[:ll,i])) , color='darkgreen', linestyle='--', lw=2.3)
    #plt.axvline(x=np.log10(np.mean(arry2[:ll,i])) , color='k',        linestyle='--', lw=2.3)
    
    plt.xlim([x1[i], x2[i]])
    plt.axvline(x=1.000 ,  color='k', linestyle='-', lw=1.6)
    plt.axvline(x=0.6989 , color='k', linestyle='--', lw=1.6)
    plt.axvline(x=0.0 ,    color='k', linestyle=':', lw=1.6)
    plt.xlabel(str(nam[i]), fontsize=19,labelpad=0.15)
    ax.set_ylabel(r"$\rm{Normalized}~\rm{Distribution}$",fontsize=19,labelpad=0.1)
    plt.xticks(fontsize=17, rotation=0)
    plt.yticks(fontsize=17, rotation=0)
    #plt.grid("True")
    #plt.grid(linestyle='dashed')
    fig=plt.gcf()
    fig.savefig("./Histo_BHs1/Histo_{0:d}.jpg".format(i),dpi=200)
    
    numa[i]=float(numa[i]*100.0/l)
    numb[i]=float(numb[i]*100.0/l)
    numc[i]=float(numc[i]*100.0/l)
    #print "No. E(\sigma< 10\%):  ", str(nam[i]),   numa[i], numb[i], numc[i]
    print ("&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&")    

###################################################################      
        
fii=open("./Histo_BHs1/results2.txt","a")
param=np.array([numa[2], numa[6], numa[7], numa[10], numa[11], numa[15], numa[14], epsi[0] ]) 
np.savetxt(fii,param.reshape((1,8)),fmt="$%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$\n") 

param=np.array([numb[2], numb[6], numb[7], numb[10], numb[11], numb[15], numb[14], epsi[1] ]) 
np.savetxt(fii,param.reshape((1,8)),fmt="$%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$\n") 

param=np.array([numc[2], numc[6], numc[7], numc[10], numc[11], numc[15], numc[14], epsi[2] ]) 
np.savetxt(fii,param.reshape((1,8)),fmt="$%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$\n") 
fii.close()         
        
###################################################################

for i in range(4):      
    plt.clf()
    fig= plt.figure(figsize=(8,6))
    ax1= plt.gca()
    plt.step(xar[:,i],m1[:,0,i],where='mid',lw=2.1,linestyle='-',alpha=0.85,color="g",label=r"$\left<\sigma_{\rm M_{\rm l}}/M_{\rm l}[\%] \right>$")
    plt.step(xar[:,i],m1[:,1,i],where='mid',lw=2.1,linestyle='--',alpha=0.85,color="b",label=r"$\left<\sigma_{\rm D_{\rm l}}/D_{\rm l}[\%]\right>$")
    plt.step(xar[:,i],m1[:,2,i],where='mid',lw=2.1,linestyle='-.',alpha=0.85,color="m",label=r"$\left<\sigma_{\mu_{\rm l}}/\mu_{\rm l}[\%]\right>$")
    plt.step(xar[:,i],m1[:,3,i],where='mid',lw=2.1,linestyle=':',alpha=0.85,color="r",label=r"$\left<\sigma_{\pi_{\rm E}}/\pi_{\rm E}[\%]\right>$")
    
    plt.scatter(xar[:,i],m1[:,0,i],marker= "o",facecolors='darkgreen',   edgecolors='g', s= 34.0)
    plt.scatter(xar[:,i],m1[:,1,i],marker= "o",facecolors='darkblue',    edgecolors='b', s= 34.0)
    plt.scatter(xar[:,i],m1[:,2,i],marker= "o",facecolors='darkmagenta', edgecolors='m', s= 34.0)
    plt.scatter(xar[:,i],m1[:,3,i],marker= "o",facecolors='darkred',     edgecolors='r', s= 34.0)
    
    plt.xticks(fontsize=19, rotation=0)
    plt.yticks(fontsize=19, rotation=0)
    
    if(i==0): plt.xlim([2.0,49.8])
    if(i==1): plt.xlim([0.0, 0.9])
    if(i==2): plt.xlim([0.0, 12.0])
    if(i==3): plt.xlim([15.5, 24.5])
    
    if(i==0): plt.xlabel(r"$M_{\rm l}[M_{\odot}]$", fontsize=18.5)
    if(i==1): plt.xlabel(r"$x_{\rm ls}$", fontsize=18.5)
    if(i==2): plt.xlabel(r"$D_{\rm s}(\rm{kpc})$", fontsize=18.5)
    if(i==3): plt.xlabel(r"$m_{\rm{base}}(\rm{mag})$", fontsize=18.5)
    
    plt.ylabel(r"$\rm{Relative}~\rm{Error}$", fontsize=18.5)
    
    plt.grid("True")
    plt.grid(linestyle='dashed')
    ax1.legend(prop={"size":18.5})
    if(i==0 or i==1):  ax1.legend().set_visible(False)
    fig3= plt.gcf()
    fig3.savefig("./Histo_BHs1/Err_{0:d}.jpg".format(i),dpi=200)
    
###################################################################    



'''
plt.clf()
fig= plt.figure(figsize=(8,6))
ax1= plt.gca()
plt.step(xar[:,0],pi1[:,0],where='mid',lw=2.1,linestyle='-',alpha=0.85,color="g",label=r"$\left<\sigma_{\rm M_{\rm l}}/M_{\rm l}[\%] \right>$")
plt.step(xar[:,0],pi1[:,1],where='mid',lw=2.1,linestyle='--',alpha=0.85,color="b",label=r"$\left<\sigma_{\rm D_{\rm l}}/D_{\rm l}[\%]\right>$")
plt.step(xar[:,0],pi1[:,2],where='mid',lw=2.1,linestyle='-.',alpha=0.85,color="m",label=r"$\left<\sigma_{\mu_{\rm l}}/\mu_{\rm l}[\%]\right>$")
    
plt.scatter(xar[:,0],pi1[:,0],marker= "o",facecolors='darkgreen',   edgecolors='g', s= 34.0)
plt.scatter(xar[:,0],pi1[:,1],marker= "o",facecolors='darkblue',    edgecolors='b', s= 34.0)
plt.scatter(xar[:,0],pi1[:,2],marker= "o",facecolors='darkmagenta', edgecolors='m', s= 34.0)
    
plt.xticks(fontsize=18, rotation=0)
plt.yticks(fontsize=18, rotation=0)
    
plt.xlim([1.5,19.0])
plt.xlabel(r"$M_{\rm l}[M_{\odot}]$", fontsize=18.5)
    
plt.grid("True")
plt.grid(linestyle='dashed')
ax1.legend(prop={"size":18.5})
fig3= plt.gcf()
fig3.savefig("./Histo_BHs1/piEdep.jpg",dpi=200)
'''   


  

