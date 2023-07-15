int EnergyFound()

{

int locals, StopFlagg, nmatchout,nmatchin, Nlimit;   

StopFlagg = FALSE;

for(locals=0;locals<ntot+5;locals++) rr[locals] =  locals *H;

for(locals=1;locals<30;locals++)     /*   Starting energy search   */
{

      Rmax=  (ntot-2)*H;
      Npoints=0;
          Rold=Rstart;
          Fold=Fstart;
          Gold=Gstart;
          Fnew=0.0;
          Furt= Func(Rold);
          p[0]=0.0;
          rr[0]=Rold;
          Npoints=0;

          for(;;)
          {
             Ruku(Ix);
             Rold=Rnew;
             Fold=Fnew;
             Gold=Gnew;

             Npoints++;

             p[Npoints]=Fnew;
             rr[Npoints]=Rnew;  

            if ( Rnew > Rmatch )  break;

         }
     
                
        
      nmatchout =  Npoints;
      DerOut=Gnew/Fnew;   /*  logarithmic derivative  */
      for(kg=0;kg<nmatchout+1;kg++  )    p[kg]=p[kg]/Fold;
                        
/*
     inward integration
*/ 
      Hsave = H;
      H=(-H);
     
     
        Rold =  rr[ntot+1] ;
        Npoints=ntot+1;
        Nlimit=Npoints;
        Fold=0.0001;
        Fnew=0.0001;
        Gold = (-Fnew) ;

        
          for(;;)
          {
             Ruku(Ix);
             Rold=Rnew; 
             Fold=Fnew;
             Gold=Gnew;

             Npoints--; 
      
             p[Npoints]=Fnew;
             rr[Npoints]=Rnew;
            if ( Rnew < (Rmatch+Hsave) )  break;
         }

       nmatchin=Npoints;
                          /*   End of RE CALCULATING  */ 
                          
     for(kg=Npoints;kg<Nlimit;kg++  )    p[kg]=p[kg]/Fold;
      Rin=Rold;
      DerIn=Gnew/Fnew;  /*  logarithmic derivative  */

      H=Hsave;
      for(kst=0;kst<Nlimit;kst++  )   pp2[kst]=p[kst]*p[kst];                   
      a=sumd(pp2,0,Nlimit,H);  
      a=sqrt(a);
      acon=1.0/a;

         Xp=Xscal*rr[0];; Ix=Xp;
         Yp=WaveY;        Iy=Yp;
         MoveTo( Ix,Iy) ; 
         
      COLOR(blueColor);
      for(kst=0;kst<ntot;kst++  )
       { p[kst]=p[kst]*acon;
        
         Xp=Xscal*rr[kst];                     Ix=Xp;
         Yp=WaveY-(p[kst]*Scal*UnitY)*13.0;     Iy=Yp;
         LineTo( Ix,Iy) ; 
       }
       
/*
       estimate of new energy after Hartree's book
*/             
       dele=   p[nmatchout-1]*p[nmatchout-1]* (DerOut-DerIn)  ;
        En=En+dele;
      
        
        
       COLOR(blueColor);
       Iy=BaseY-En*Vsc*Scal*UnitY;     
       MoveTo(0,Iy);
       LineTo(Scale1,Iy);             /*   Energy Line  */
  
      
     if(  (  (dele/En)<0.0001  ) && ((dele/En)> (-0.0001))  )  break;
     if(  StopFlagg )    break; 

 
 
  }   /*   Ending energy search   */
  
  
        Xp=Xscal*Rstart; Ix=Xp;
         Yp=WaveY;      Iy=Yp;
         MoveTo( Ix,Iy) ; 

       COLOR(redColor);
      for(kst=0;kst<ntot;kst++  )
       { 
         Xp=Xscal*kst*H;                        Ix=Xp;
         Yp=WaveY+(p[kst]*Scal*UnitY)*16.0;     Iy=Yp;  /* inversed */
         LineTo( Ix,Iy) ; 
       }        
}
      
double sumd(vec,m1,m2,xh) double *vec, xh; int m1,m2;
       /*
            simpson rule,  from  m1-th  to  m2-th  point
        */
 {     double summl; int m3,m4, kst; 
      
      summl=vec[m1]+vec[m2];
      m3=m2+1;
      m4=m1+1;
      for(kst=m4;kst<m2;   kst+=2  )
             summl=summl+4.0*vec[kst];
   if( (m2-m1) > 2 ) 
    {  
      m3=m1+2;
      m4=m2-1;   
      for(kst=m3;kst<m4;   kst+=2  )
          summl=summl+2.0*vec[kst];
    }
   summl=summl*xh/3;
   return (summl);
 }

