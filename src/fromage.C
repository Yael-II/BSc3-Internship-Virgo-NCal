/*-----------FROMAGE-----------*/
/*---Finite element analysis of ROtating MAsses for Gravitational Effects---*/

/*-------------------------*/
/*    ___                  */
/*   |   |       O         */
/*   |   |        \        */
/*   |   |         \       */
/*   |   |          \      */
/*   |___|           O     */
/*                         */   
/*-------------------------*/

/*---Developer: Dimitri Estevez, dimitri.estevez@iphc.cnrs.fr---*/
/*---Version v1r2 from May 2022---*/

// This software is a simulation tool to compute gravitational effects of rotating masses
// on an object.
// It has been primarily developed to be used for calibration purpose of gravitational-wave
// interferometric detectors (currently LIGO-Virgo-KAGRA collaboration) to simulate 
// NCal (Newtonian Calibrator) or GCal (Gravity Calibrator) effects on the end mirrors
// (or end test masses) of the interferometers.

// The definition of the objects is done using a configuration file (.cfg file) using
// specific keywords defined in ~/FROMAGE/v1r2/doc/FroDoc.pdf


/*-----------------------------*/

#include <iostream>
#include <iomanip>
#include <string.h>
#include <math.h>
#include <vector>

using namespace std;

struct FroRoot {
  vector<int> nelements; 		   // list of number of elements per object 
  vector<vector<int> > grid; 		   // linked list of parameters in the grid 
  vector<vector<vector<double> > > mirror; // linked list of mirror objects 
  vector<vector<vector<double> > > rotor;  // linked list of rotor objects 
  int rescale_0=1; 			   // number to rescale the first parameter of the grids
  int rescale_1=1; 			   // number to rescale the second parameter of the grids
  int rescale_2=1; 			   // number to rescale the third parameter of the grids
  int nmir_obj; 			   // total number of mirror objects defined
  int nrot_obj; 			   // total number of rotor objects defined
  int nstep; 			   	   // number of steps for the rotor motion 
  int rotorindex = 0; 			   // index defining the rotor 
  int dir; 			 	   // direction of the rotor axis 
  int signal = 0; 			   // number N to compute the amplitude residuals with respect to a Fourier series at N theta 
  int debug = 1; 			   // level of debug
  double theta; 			   // angle to spin the rotor sectors 
  double G = 6.67430e-11; 		   // gravitational constant in SI unit 
  double arm_length; 			   // arm length of the interferometer
  double d; 			      	   // distance from the center of the mirror to the center of the rotor (m)
  double phi; 				   // angle of the rotor with respect to the ITF beam axis 
  double offset; 			   // offset of the rotor with respect to the plane of the ITF 
  double force[3]; 			   // list of forces projected on (x,y,z) 
  double torque[2]; 			   // list of torques along the orthogonal axes (y,z) to the interferometer beam (x)  
  double Iy_cyl;			   // moment of inertia of the mirror around the y axis
  double Iz_cyl;			   // moment of inertia of the mirror around the z axis
  double beam_y = 0.; 			   // offset of the ITF beam on the y-axis 
  double beam_z = 0.; 			   // offset of the ITF beam on the z-axis 
  double mirror_mass = 0.; 		   // total mass of the mirror 
  char userfile0[256]; 			   // file name where to save the displacement values 
  char userfile1[256]; 			   // file name where to save the mirror elements positions 
  char userfile2[256]; 			   // file name where to save the rotor elements positions
  };

struct FroRoot Fromage;

/*-----------------------------*/

void DefineGrid(char* line, int cpt_grid){
  int a,b,c;
  if(sscanf(line,"GRID_SIZE %d %d %d", &a, &b, &c)!=3){
    cout<<"ERROR: Wrong number of parameters for GRID_SIZE!"<<endl;
    exit(-4);
  }
  if(a<=0 || b<=0 || c<=0){
    cout<<"ERROR: Grid parameters must be strictly positive!"<<endl;
    exit(-5);
  }
  Fromage.grid.push_back(vector<int>(3));
  Fromage.grid[cpt_grid][0] = a*Fromage.rescale_0;
  Fromage.grid[cpt_grid][1] = b*Fromage.rescale_1;
  Fromage.grid[cpt_grid][2] = c*Fromage.rescale_2;
  return;
}

/*-----------------------------*/

void DefineCylinder(char* line, int cpt_grid, int cpt_cyl){

  Fromage.nelements.push_back(Fromage.grid[cpt_grid-1][0]*Fromage.grid[cpt_grid-1][1]*Fromage.grid[cpt_grid-1][2]);

  double density, rmin, rmax, thick, alpha, distObj_r, distObj_x, angle;
  double x_pos, distG, mass_partcyl, R2, r2, R3, r3;
  int cpt = 0;
  int nx = Fromage.grid[cpt_grid-1][0], nt = Fromage.grid[cpt_grid-1][1], nr = Fromage.grid[cpt_grid-1][2];

  if(sscanf(line,"CYLINDER %lf %lf %lf %lf %lf %lf %lf %lf", &density, &rmin, &rmax, &thick, &alpha, &distObj_r, &distObj_x, &angle)!=8){
    cout<<"ERROR: Wrong number of parameters for CYLINDER!"<<endl;
    exit(-6);
  }
  if(rmin<0 || rmax<=rmin || thick<0 || alpha<=0){
    cout<<"ERROR: The cylinder parameters are not well defined! Please refer to the documentation."<<endl;
    exit(-7);
  }
  alpha = alpha*M_PI/180.;
  angle = angle*M_PI/180.;
  if(cpt_cyl==0 && Fromage.rotorindex==0){
    Fromage.Iy_cyl = 0.25*density*rmax*rmax*thick*M_PI*(rmax*rmax+thick*thick/3.);
    Fromage.Iz_cyl = 0.25*density*rmax*rmax*thick*M_PI*(rmax*rmax+thick*thick/3.);
  }
  if(Fromage.rotorindex==0) Fromage.mirror.push_back(vector<vector<double> >());
  else Fromage.rotor.push_back(vector<vector<double> >());

  for(int a=0;a<nx;a++){
    R2 = rmax*rmax;
    R3 = R2*rmax;

    if(nx%2==0){
      if(a<nx/2) x_pos = thick/nx/2.+a*thick/nx;
      else x_pos = -thick/nx/2.-(a-nx/2.)*thick/nx;
    }
    else{
      if(a<nx/2.) x_pos = a*thick/nx;
      else x_pos = -(a-nx/2)*thick/nx;
    }

    for(int b=0; b<nr; b++){

      r2 = 1./(double)nr*((nr-(b+1))*rmax*rmax+(b+1)*rmin*rmin);
      r3 = r2*sqrt(1./(double)nr*((nr-(b+1))*rmax*rmax+(b+1)*rmin*rmin));
      distG = 4./3./(alpha/nt)*sin(alpha/2./nt)*(R3-r3)/(R2-r2);

      mass_partcyl = alpha/2./nt*thick/nx*(R2-r2)*density;
      R2 = r2;
      R3 = r3;

      for(int c=0; c<nt; c++){
        if(Fromage.rotorindex==0){
  	  Fromage.mirror[cpt_cyl].push_back(vector<double>(4));
 	  Fromage.mirror[cpt_cyl][cpt][0] = x_pos+distObj_x;
          Fromage.mirror[cpt_cyl][cpt][1] = distG*cos(alpha/2./nt+c*alpha/nt-alpha/2.+angle)+distObj_r*cos(angle);
          Fromage.mirror[cpt_cyl][cpt][2] = distG*sin(alpha/2./nt+c*alpha/nt-alpha/2.+angle)+distObj_r*sin(angle);
 	  Fromage.mirror[cpt_cyl][cpt][3] = mass_partcyl;
          cpt++;
        }
	else{
  	  Fromage.rotor[cpt_cyl].push_back(vector<double>(4));
 	  Fromage.rotor[cpt_cyl][cpt][0] = x_pos+distObj_x;
          Fromage.rotor[cpt_cyl][cpt][1] = distG*cos(alpha/2./nt+c*alpha/nt-alpha/2.+angle)+distObj_r*cos(angle);
          Fromage.rotor[cpt_cyl][cpt][2] = distG*sin(alpha/2./nt+c*alpha/nt-alpha/2.+angle)+distObj_r*sin(angle);
 	  Fromage.rotor[cpt_cyl][cpt][3] = mass_partcyl;
          cpt++;
        }
        
      }
    }
  }
  return;
}

/*-----------------------------*/

void DefineCylinderCut(char* line, int cpt_grid, int cpt_cyl){

  Fromage.nelements.push_back(Fromage.grid[cpt_grid-1][0]);

  double x_pos, density, radius, thick, cut, distObj_x, angle;
  double mass_cut, cutG;
  int cpt = 0, nx = Fromage.grid[cpt_grid-1][0];

  if(sscanf(line,"CUT_CYL %lf %lf %lf %lf %lf %lf", &density, &radius, &thick, &cut, &distObj_x, &angle)!=6){
    cout<<"ERROR: Wrong number of parameters for CUT_CYL!"<<endl;
    exit(-8);
  }
  if(radius<=0 || thick<0 || cut<=0){
    cout<<"ERROR: The flats parameters on the cylinder are not well defined! Please refer to the documentation."<<endl;
    exit(-9);
  }
  angle = angle*M_PI/180.;
  if(Fromage.rotorindex==0) Fromage.mirror.push_back(vector<vector<double> >());
  else Fromage.rotor.push_back(vector<vector<double> >());

  mass_cut = -density*radius*radius/2.*thick*(2.*asin(cut/2./radius)-sin(2.*asin(cut/2./radius)))/nx;
  cutG = 4./3.*radius*cut/2./radius*cut/2./radius*cut/2./radius/(2.*asin(cut/2./radius)-sin(2.*asin(cut/2./radius)));
 
  for(int a=0;a<nx;a++){
    if(nx%2==0){
       if(a<nx/2) x_pos = thick/nx/2.+a*thick/nx;
       else x_pos = -thick/nx/2.-(a-nx/2.)*thick/nx;
     }
     else{
       if(a<nx/2.) x_pos = a*thick/nx;
       else x_pos = -(a-nx/2)*thick/nx;
     }
    if(Fromage.rotorindex==0){
      Fromage.mirror[cpt_cyl].push_back(vector<double>(4));
      Fromage.mirror[cpt_cyl][cpt][0] = x_pos+distObj_x;
      Fromage.mirror[cpt_cyl][cpt][1] = cutG*cos(angle);
      Fromage.mirror[cpt_cyl][cpt][2] = cutG*sin(angle);
      Fromage.mirror[cpt_cyl][cpt][3] = mass_cut;
      cpt++;
    }
    else{
      Fromage.rotor[cpt_cyl].push_back(vector<double>(4));
      Fromage.rotor[cpt_cyl][cpt][0] = x_pos+distObj_x;
      Fromage.rotor[cpt_cyl][cpt][1] = cutG*cos(angle);
      Fromage.rotor[cpt_cyl][cpt][2] = cutG*sin(angle);
      Fromage.rotor[cpt_cyl][cpt][3] = mass_cut;
      cpt++;
    }
  }

  return;
}

/*-----------------------------*/

void DefineCuboid(char* line, int cpt_grid, int cpt_cub){

  Fromage.nelements.push_back(Fromage.grid[cpt_grid-1][0]*Fromage.grid[cpt_grid-1][1]*Fromage.grid[cpt_grid-1][2]);

  double density, length, thick, height;
  double x_pos, y_pos, z_pos, mass_partcub, distObj_x, distObj_y, distObj_z;
  int cpt = 0;
  int nx = Fromage.grid[cpt_grid-1][0], nt = Fromage.grid[cpt_grid-1][1], nr = Fromage.grid[cpt_grid-1][2];

  if(sscanf(line,"CUBOID %lf %lf %lf %lf %lf %lf %lf", &density, &length, &thick, &height, &distObj_x, &distObj_y, &distObj_z)!=7){
    cout<<"ERROR: Wrong number of parameters for CUBOID!"<<endl;
    exit(-10);
  }
  if(length<0 || height<0 || thick<0){
    cout<<"ERROR: The cuboid parameters are not well defined! Please refer to the documentation."<<endl;
    exit(-11);
  }
  if(Fromage.rotorindex==0) Fromage.mirror.push_back(vector<vector<double> >());
  else Fromage.rotor.push_back(vector<vector<double> >());


  for(int a=0;a<nx;a++){

    if(nx%2==0){
      if(a<nx/2) y_pos = thick/nx/2.+a*thick/nx;
      else y_pos = -thick/nx/2.-(a-nx/2.)*thick/nx;
    }
    else{
      if(a<nx/2.) y_pos = a*thick/nx;
      else y_pos = -(a-nx/2)*thick/nx;
    }

    for(int b=0; b<nr; b++){

      if(nr%2==0){
        if(b<nr/2) z_pos = height/nr/2.+b*height/nr;
        else z_pos = -height/nr/2.-(b-nr/2.)*height/nr;
      }
      else{
        if(b<nx/2.) z_pos = b*height/nr;
        else z_pos = -(b-nr/2)*height/nr;
      }

      mass_partcub = density*length*height*thick/(nx*nr*nt);

      for(int c=0; c<nt; c++){
        if(nt%2==0){
          if(c<nt/2) x_pos = length/nt/2.+c*length/nt;
          else x_pos = -length/nt/2.-(c-nt/2.)*length/nt;
        }
        else{
          if(c<nt/2.) x_pos = c*length/nt;
          else x_pos = -(c-nt/2)*length/nt;
        }
        if(Fromage.rotorindex==0){
  	Fromage.mirror[cpt_cub].push_back(vector<double>(4));
 	Fromage.mirror[cpt_cub][cpt][0] = x_pos+distObj_x;
        Fromage.mirror[cpt_cub][cpt][1] = y_pos+distObj_y;
        Fromage.mirror[cpt_cub][cpt][2] = z_pos+distObj_z;
 	Fromage.mirror[cpt_cub][cpt][3] = mass_partcub;
        cpt++;  
        }  
        else{
  	Fromage.rotor[cpt_cub].push_back(vector<double>(4));
 	Fromage.rotor[cpt_cub][cpt][0] = x_pos+distObj_x;
        Fromage.rotor[cpt_cub][cpt][1] = y_pos+distObj_y;
        Fromage.rotor[cpt_cub][cpt][2] = z_pos+distObj_z;
 	Fromage.rotor[cpt_cub][cpt][3] = mass_partcub;
        cpt++;  
        }  
      }
    }
  }
  return;
}

/*-----------------------------*/

void DefineOuterFillet(char* line, int cpt_grid, int cpt_cyl){

  Fromage.nelements.push_back(Fromage.grid[cpt_grid-1][0]);

  double x_pos, density, radius, thick, distObj_x, distObj_r, alpha_demi, angle;
  double mass_fil, filG_y, filG_z, area, beta, gamma, g_beta, g_gamma, area_beta, area_gamma;
  int cpt = 0, nx = Fromage.grid[cpt_grid-1][0];

  if(sscanf(line,"OUTER_FILLET %lf %lf %lf %lf %lf %lf %lf", &density, &distObj_r, &thick, &distObj_x, &radius, &alpha_demi, &angle)!=7){
    cout<<"ERROR: Wrong number of parameters for OUTER_FILLET!"<<endl;
    exit(-12);
  }
  if(distObj_r<0 || thick<0 || alpha_demi==0 || radius<0){
    cout<<"ERROR: The cylinder outer fillet parameters are not well defined! Please refer to the documentation."<<endl;
    exit(-13);
  }

  alpha_demi = alpha_demi*M_PI/180.;
  angle = angle*M_PI/180.;

  beta = asin(radius/(distObj_r+radius));
  gamma = M_PI/2.-beta;

  area = 0.5*(radius*distObj_r*sqrt(1.+2.*radius/distObj_r)-distObj_r*distObj_r*beta-radius*radius*gamma);

  mass_fil = area*thick*density/nx;

  g_beta = 2.*distObj_r*sin(beta/2.)/(3.*beta/2.);
  g_gamma = 2.*radius*sin(gamma/2.)/(3.*gamma/2.);

  area_beta = beta*distObj_r*distObj_r/2.;
  area_gamma = gamma*radius*radius/2.;

  double y1 = 0., z1 = 0.;
  double y2 = distObj_r*sqrt(1.+2.*radius/distObj_r), z2 = radius;
  double y3 = distObj_r*sqrt(1.+2.*radius/distObj_r), z3 = 0.;

  double area_tr = 0.5*z2*y3;

  filG_y = 1./area*(1./3.*(y1+y2+y3)*area_tr-area_beta*g_beta*cos(beta/2.)-area_gamma*(y2-g_gamma*sin(gamma/2.)));
  filG_z = 1./area*(1./3.*(z1+z2+z3)*area_tr-area_beta*g_beta*sin(beta/2.)-area_gamma*(radius-g_gamma*cos(gamma/2.)));

  if(Fromage.rotorindex==0) Fromage.mirror.push_back(vector<vector<double> >());
  else Fromage.rotor.push_back(vector<vector<double> >());

  if(alpha_demi<0) filG_z=-filG_z;

  for(int a=0;a<nx;a++){
    if(nx%2==0){
      if(a<nx/2) x_pos = thick/nx/2.+a*thick/nx;
      else x_pos = -thick/nx/2.-(a-nx/2.)*thick/nx;
    }
    else{
      if(a<nx/2.) x_pos = a*thick/nx;
      else x_pos = -(a-nx/2)*thick/nx;
    }
    if(Fromage.rotorindex==0){
      Fromage.mirror[cpt_cyl].push_back(vector<double>(4));
      Fromage.mirror[cpt_cyl][cpt][0] = x_pos+distObj_x;
      Fromage.mirror[cpt_cyl][cpt][1] = filG_y*cos(alpha_demi+angle)-filG_z*sin(alpha_demi+angle);
      Fromage.mirror[cpt_cyl][cpt][2] = filG_y*sin(alpha_demi+angle)+filG_z*cos(alpha_demi+angle);
      Fromage.mirror[cpt_cyl][cpt][3] = mass_fil;
      cpt++;
    }
    else{
      Fromage.rotor[cpt_cyl].push_back(vector<double>(4));
      Fromage.rotor[cpt_cyl][cpt][0] = x_pos+distObj_x;
      Fromage.rotor[cpt_cyl][cpt][1] = filG_y*cos(alpha_demi+angle)-filG_z*sin(alpha_demi+angle);
      Fromage.rotor[cpt_cyl][cpt][2] = filG_y*sin(alpha_demi+angle)+filG_z*cos(alpha_demi+angle);
      Fromage.rotor[cpt_cyl][cpt][3] = mass_fil;
      cpt++;
    }
  }

  return;
}

/*-----------------------------*/

void DefineInnerFillet(char* line, int cpt_grid, int cpt_cyl){

  Fromage.nelements.push_back(Fromage.grid[cpt_grid-1][0]);

  double x_pos, density, radius, thick, distObj_x, distObj_r, alpha_demi, angle;
  double mass_fil, filG_y, filG_z, area, beta, gamma, g_beta, g_gamma, area_beta, area_gamma;
  int cpt = 0, nx = Fromage.grid[cpt_grid-1][0];

  if(sscanf(line,"INNER_FILLET %lf %lf %lf %lf %lf %lf %lf", &density, &distObj_r, &thick, &distObj_x, &radius, &alpha_demi, &angle)!=7){
    cout<<"ERROR: Wrong number of parameters for INNER_FILLET!"<<endl;
    exit(-14);
  }
  if(distObj_r<0 || thick<0 || alpha_demi==0 || radius<0){
    cout<<"ERROR: The cylinder inner fillet parameters are not well defined! Please refer to the documentation."<<endl;
    exit(-15);
  }

  alpha_demi = alpha_demi*M_PI/180.;
  angle = angle*M_PI/180.;

  beta = asin(radius/(distObj_r-radius));
  gamma = M_PI/2.+beta;

  area = 0.5*(distObj_r*distObj_r*beta-radius*distObj_r*sqrt(1.-2.*radius/distObj_r)-radius*radius*gamma);

  mass_fil = area*thick*density/nx;

  g_beta = 2.*distObj_r*sin(beta/2.)/(3.*beta/2.);
  g_gamma = 2.*radius*sin(gamma/2.)/(3.*gamma/2.);

  area_beta = beta*distObj_r*distObj_r/2.;
  area_gamma = gamma*radius*radius/2.;

  double y1 = 0., z1 = 0.;
  double y2 = distObj_r*sqrt(1.-2.*radius/distObj_r), z2 = radius;
  double y3 = distObj_r*sqrt(1.-2.*radius/distObj_r), z3 = 0.;

  double area_tr = 0.5*z2*y3;

  filG_y = 1./area*(area_beta*g_beta*cos(beta/2.)-1./3.*(y1+y2+y3)*area_tr-area_gamma*(y2+g_gamma*sin(gamma/2.)));
  filG_z = 1./area*(area_beta*g_beta*sin(beta/2.)-1./3.*(z1+z2+z3)*area_tr-area_gamma*(radius-g_gamma*cos(gamma/2.)));

  if(Fromage.rotorindex==0) Fromage.mirror.push_back(vector<vector<double> >());
  else Fromage.rotor.push_back(vector<vector<double> >());

  if(alpha_demi<0) filG_z=-filG_z;

  for(int a=0;a<nx;a++){
    if(nx%2==0){
      if(a<nx/2) x_pos = thick/nx/2.+a*thick/nx;
      else x_pos = -thick/nx/2.-(a-nx/2.)*thick/nx;
    }
    else{
      if(a<nx/2.) x_pos = a*thick/nx;
      else x_pos = -(a-nx/2)*thick/nx;
    }
    if(Fromage.rotorindex==0){
      Fromage.mirror[cpt_cyl].push_back(vector<double>(4));
      Fromage.mirror[cpt_cyl][cpt][0] = x_pos+distObj_x;
      Fromage.mirror[cpt_cyl][cpt][1] = filG_y*cos(alpha_demi+angle)-filG_z*sin(alpha_demi+angle);
      Fromage.mirror[cpt_cyl][cpt][2] = filG_y*sin(alpha_demi+angle)+filG_z*cos(alpha_demi+angle);
      Fromage.mirror[cpt_cyl][cpt][3] = mass_fil;
      cpt++;
    }
    else{
      Fromage.rotor[cpt_cyl].push_back(vector<double>(4));
      Fromage.rotor[cpt_cyl][cpt][0] = x_pos+distObj_x;
      Fromage.rotor[cpt_cyl][cpt][1] = filG_y*cos(alpha_demi+angle)-filG_z*sin(alpha_demi+angle);
      Fromage.rotor[cpt_cyl][cpt][2] = filG_y*sin(alpha_demi+angle)+filG_z*cos(alpha_demi+angle);
      Fromage.rotor[cpt_cyl][cpt][3] = mass_fil;
      cpt++;
    }
  }

  return;
}

/*-----------------------------*/

void ComputeFourierForce(vector<vector<double> > buffer, double force_tot){
  double coeffRe, coeffIm, res, res_max = 0;
  vector<double> CoeffFour(Fromage.signal+1);
  vector<double> ArgFour(Fromage.signal+1);

  for(int i=0;i<=Fromage.signal;i++){
    coeffRe=0;
    coeffIm=0;
    for(int j=0;j<Fromage.nstep;j++){
      coeffRe += buffer[1][j]*cos(i*buffer[0][j]);
      coeffIm += -buffer[1][j]*sin(i*buffer[0][j]);
    }
    CoeffFour[i]=sqrt(coeffRe*coeffRe+coeffIm*coeffIm);
    if(coeffRe<0 && coeffIm>=0){
      ArgFour[i] = atan(coeffIm/coeffRe)*180./M_PI+180.;
    }
    else if(coeffRe<0 && coeffIm<0){
      ArgFour[i] = atan(coeffIm/coeffRe)*180./M_PI-180.;
    }
    else if(coeffRe==0 && coeffIm>0){
      ArgFour[i] = 90.;
    }
    else if(coeffRe==0 && coeffIm<0){
      ArgFour[i] = -90.;
    }
    else if(coeffRe==0 && coeffIm==0){
      ArgFour[i] = 0.;
    }
    else{
      ArgFour[i] = atan(coeffIm/coeffRe)*180./M_PI;
    }

    if(Fromage.debug>=1){
      if(i==0) cout<<setprecision(40)<<"Force offset is "<<1./Fromage.nstep*CoeffFour[i]<<" N"<<", Phase shift is "<<ArgFour[i]<<" deg"<<endl;
      else cout<<setprecision(40)<<"Force at "<<i<<"f is "<<2./Fromage.nstep*CoeffFour[i]<<" N_pk (strain h = "<<2./Fromage.nstep*CoeffFour[i]/(Fromage.mirror_mass*4*M_PI*M_PI*Fromage.arm_length)<<"/("<<i<<"f_rot)^2), Phase shift is "<<ArgFour[i]<<" deg"<<endl;
    }
  }
  for(int i=0;i<Fromage.nstep;i++){
    res = buffer[1][i];
    for(int j=0;j<=Fromage.signal;j++){
      if(j==0) res -= CoeffFour[j]*cos(j*buffer[0][i]+ArgFour[j]*M_PI/180.)/Fromage.nstep;
      else res -= 2*CoeffFour[j]*cos(j*buffer[0][i]+ArgFour[j]*M_PI/180.)/Fromage.nstep;
    }
    if(i==0) res_max = fabs(res);
    if(fabs(res)>res_max) res_max = fabs(res);
  }
  cout<<endl;
  cout<<"Residual between numerical force and Fourier series of order "<<Fromage.signal<<" is "<<res_max/force_tot*100.<<"%"<<endl;
  cout<<endl;

  return;
}

/*-----------------------------*/

void ComputeFourierTorque(vector<vector<double> > buffer, double force_tot, double tory_tot, double torz_tot){
  double coeffRe, coeffIm, res, res_max = 0;
  vector<double> CoeffFour(Fromage.signal+1);
  vector<double> ArgFour(Fromage.signal+1);
  double coeffRey, coeffImy, resy, resy_max = 0;
  vector<double> CoeffFoury(Fromage.signal+1);
  vector<double> ArgFoury(Fromage.signal+1);
  double coeffRez, coeffImz, resz, resz_max = 0;
  vector<double> CoeffFourz(Fromage.signal+1);
  vector<double> ArgFourz(Fromage.signal+1);
  vector<vector<double> > amp_tot(Fromage.nstep,vector<double>(Fromage.signal+1));
  vector<double> max_tot(Fromage.signal+1);
  vector<double> min_tot(Fromage.signal+1);

  for(int i=0;i<=Fromage.signal;i++){
    coeffRe=0;
    coeffIm=0;
    coeffRey=0;
    coeffImy=0;
    coeffRez=0;
    coeffImz=0;
    for(int j=0;j<Fromage.nstep;j++){
      coeffRe += buffer[1][j]*cos(i*buffer[0][j]);
      coeffIm += -buffer[1][j]*sin(i*buffer[0][j]);
      coeffRey += buffer[2][j]*cos(i*buffer[0][j]);
      coeffImy += -buffer[2][j]*sin(i*buffer[0][j]);
      coeffRez += buffer[3][j]*cos(i*buffer[0][j]);
      coeffImz += -buffer[3][j]*sin(i*buffer[0][j]);
    }
    CoeffFour[i]=sqrt(coeffRe*coeffRe+coeffIm*coeffIm);
    CoeffFoury[i]=sqrt(coeffRey*coeffRey+coeffImy*coeffImy);
    CoeffFourz[i]=sqrt(coeffRez*coeffRez+coeffImz*coeffImz);

    if(coeffRe<0 && coeffIm>=0){
      ArgFour[i] = atan(coeffIm/coeffRe)*180./M_PI+180.;
    }
    else if(coeffRe<0 && coeffIm<0){
      ArgFour[i] = atan(coeffIm/coeffRe)*180./M_PI-180.;
    }
    else if(coeffRe==0 && coeffIm>0){
      ArgFour[i] = 90.;
    }
    else if(coeffRe==0 && coeffIm<0){
      ArgFour[i] = -90.;
    }
    else if(coeffRe==0 && coeffIm==0){
      ArgFour[i] = 0.;
    }
    else{
      ArgFour[i] = atan(coeffIm/coeffRe)*180./M_PI;
    }

    if(coeffRey<0 && coeffImy>=0){
      ArgFoury[i] = atan(coeffImy/coeffRey)*180./M_PI+180.;
    }
    else if(coeffRey<0 && coeffImy<0){
      ArgFoury[i] = atan(coeffImy/coeffRey)*180./M_PI-180.;
    }
    else if(coeffRey==0 && coeffImy>0){
      ArgFoury[i] = 90.;
    }
    else if(coeffRey==0 && coeffImy<0){
      ArgFoury[i] = -90.;
    }
    else if(coeffRey==0 && coeffImy==0){
      ArgFoury[i] = 0.;
    }
    else{
      ArgFoury[i] = atan(coeffImy/coeffRey)*180./M_PI;
    }

    if(coeffRez<0 && coeffImz>=0){
      ArgFourz[i] = atan(coeffImz/coeffRez)*180./M_PI+180.;
    }
    else if(coeffRez<0 && coeffImz<0){
      ArgFourz[i] = atan(coeffImz/coeffRez)*180./M_PI-180.;
    }
    else if(coeffRez==0 && coeffImz>0){
      ArgFourz[i] = 90.;
    }
    else if(coeffRez==0 && coeffImz<0){
      ArgFourz[i] = -90.;
    }
    else if(coeffRez==0 && coeffImz==0){
      ArgFourz[i] = 0.;
    }
    else{
      ArgFourz[i] = atan(coeffImz/coeffRez)*180./M_PI;
    }

    if(Fromage.debug>=1){
      if(i==0){
        cout<<"Force offset is "<<1./Fromage.nstep*CoeffFour[i]<<" N"<<", Phase shift is "<<ArgFour[i]<<" deg"<<endl;
        cout<<"Torque offset around y is "<<1./Fromage.nstep*CoeffFoury[i]<<" N.m"<<", Phase shift is "<<ArgFour[i]<<" deg"<<endl;
        cout<<"Torque offset around z is "<<1./Fromage.nstep*CoeffFourz[i]<<" N.m"<<", Phase shift is "<<ArgFour[i]<<" deg"<<endl;
      }
      else{
        cout<<endl;
	cout<<"Force at "<<i<<"f is "<<2./Fromage.nstep*CoeffFour[i]<<" N_pk (strain h = "<<2./Fromage.nstep*CoeffFour[i]/(Fromage.mirror_mass*4*M_PI*M_PI*Fromage.arm_length)<<"/("<<i<<"f_rot)^2, Phase shift is "<<ArgFour[i]<<" deg"<<endl;
        cout<<"Torque around y at "<<i<<"f is "<<2./Fromage.nstep*CoeffFoury[i]<<" N_pk.m (strain h = "<<2./Fromage.nstep*CoeffFoury[i]/(Fromage.Iy_cyl*4*M_PI*M_PI*Fromage.arm_length)*Fromage.beam_z<<"/("<<i<<"f_rot)^2), Phase shift is "<<ArgFoury[i]<<" deg"<<endl;
        cout<<"Torque around z at "<<i<<"f is "<<2./Fromage.nstep*CoeffFourz[i]<<" N_pk.m (strain h = "<<2./Fromage.nstep*CoeffFourz[i]/(Fromage.Iz_cyl*4*M_PI*M_PI*Fromage.arm_length)*Fromage.beam_y<<"/("<<i<<"f_rot)^2), Phase shift is "<<ArgFourz[i]<<" deg"<<endl;
      }
    }
  }
  for(int i=0;i<Fromage.nstep;i++){
    res = buffer[1][i];
    resy = buffer[2][i];
    resz = buffer[3][i];
    for(int j=0;j<=Fromage.signal;j++){
      if(j==0){
        res -= CoeffFour[j]*cos(j*buffer[0][i]+ArgFour[j]*M_PI/180.)/Fromage.nstep;
        resy -= CoeffFoury[j]*cos(j*buffer[0][i]+ArgFoury[j]*M_PI/180.)/Fromage.nstep;
        resz -= CoeffFourz[j]*cos(j*buffer[0][i]+ArgFourz[j]*M_PI/180.)/Fromage.nstep;
      }
      else{
        res -= 2*CoeffFour[j]*cos(j*buffer[0][i]+ArgFour[j]*M_PI/180.)/Fromage.nstep;
        resy -= 2*CoeffFoury[j]*cos(j*buffer[0][i]+ArgFoury[j]*M_PI/180.)/Fromage.nstep;
        resz -= 2*CoeffFourz[j]*cos(j*buffer[0][i]+ArgFourz[j]*M_PI/180.)/Fromage.nstep;
        amp_tot[i][j] = 2./(Fromage.nstep*4*M_PI*M_PI)*(CoeffFour[j]*cos(j*buffer[0][i]+ArgFour[j]*M_PI/180.)/Fromage.mirror_mass+CoeffFoury[j]*cos(j*buffer[0][i]+ArgFoury[j]*M_PI/180.)/Fromage.Iy_cyl*Fromage.beam_z+CoeffFourz[j]*cos(j*buffer[0][i]+ArgFourz[j]*M_PI/180.)/Fromage.Iz_cyl*Fromage.beam_y);
      }
      if(amp_tot[i][j]>=max_tot[j]) max_tot[j]=amp_tot[i][j];
      if(amp_tot[i][j]<=min_tot[j]) min_tot[j]=amp_tot[i][j];
    }
    if(i==0){
      res_max = fabs(res);
      resy_max = fabs(resy);
      resz_max = fabs(resz);
    }
    if(fabs(res)>res_max) res_max = fabs(res);
    if(fabs(resy)>resy_max) resy_max = fabs(resy);
    if(fabs(resz)>resz_max) resz_max = fabs(resz);
  }
  cout<<endl;
  for(int k=1;k<=Fromage.signal;k++) cout<<"Strain h with force and torques at "<<k<<"f is "<<(max_tot[k]-min_tot[k])/2./Fromage.arm_length<<"/("<<k<<"f_rot)^2"<<endl;

  cout<<endl;
  cout<<"Residual between numerical force and Fourier series of order "<<Fromage.signal<<" is "<<res_max/force_tot*100.<<"%"<<endl;
  cout<<"Residual between numerical torque around y (resp. z) and Fourier series of order "<<Fromage.signal<<" is "<<resy_max/tory_tot*100.<<"% (resp. "<<resz_max/torz_tot*100.<<"%)"<<endl;
  cout<<endl;

  return;
}

/*-----------------------------*/

void ReadCfg(char* filename){

  char* line = NULL;
  size_t len = 0;
  int cpt_grid = 0, cpt_obj = 0, cpt_mir = 0, cpt_rot = 0;
  double x, y;

  FILE *fp = fopen(filename,"r");
  if(fp == NULL) {
    printf("File does not exist!");
    exit(-1);
  }

  while(getline(&line,&len,fp) !=-1){
    if(strncmp(line,"#",strlen("#"))==0) continue;
    if(strncmp(line,"RESCALE_GRID",strlen("RESCALE_GRID"))==0){
      if(sscanf(line,"RESCALE_GRID %d %d %d", &Fromage.rescale_0, &Fromage.rescale_1, &Fromage.rescale_2)!=3){
        cout<<"ERROR: Wrong number of parameters for RESCALE_GRID!"<<endl;
        exit(-16);
      }
      if(Fromage.rescale_0<=0 || Fromage.rescale_1<=0 || Fromage.rescale_2<=0){
        cout<<"ERROR: Rescale grid parameters must be strictly positive!"<<endl;
        exit(-17);
      }
    }
    if(strncmp(line,"GRID_SIZE",strlen("GRID_SIZE"))==0){
      DefineGrid(line, cpt_grid);
      cpt_grid++;
    }
    if(strncmp(line,"CYLINDER",strlen("CYLINDER"))==0){
      DefineCylinder(line,cpt_grid,cpt_obj);
      if(Fromage.rotorindex==0){ cpt_mir++; cpt_obj=cpt_mir; }
      else{ cpt_rot++; cpt_obj=cpt_rot; }
    }
    if(strncmp(line,"CUT_CYL",strlen("CUT_CYL"))==0){
      DefineCylinderCut(line,cpt_grid,cpt_obj);
      if(Fromage.rotorindex==0){ cpt_mir++; cpt_obj=cpt_mir; }
      else{ cpt_rot++; cpt_obj=cpt_rot; }
    }
    if(strncmp(line,"CUBOID",strlen("CUBOID"))==0){
      DefineCuboid(line,cpt_grid,cpt_obj);
      if(Fromage.rotorindex==0){ cpt_mir++; cpt_obj=cpt_mir; }
      else{ cpt_rot++; cpt_obj=cpt_rot; }
    }
    if(strncmp(line,"OUTER_FILLET",strlen("OUTER_FILLET"))==0){
      DefineOuterFillet(line,cpt_grid,cpt_obj);
      if(Fromage.rotorindex==0){ cpt_mir++; cpt_obj=cpt_mir; }
      else{ cpt_rot++; cpt_obj=cpt_rot; }
    }
    if(strncmp(line,"INNER_FILLET",strlen("INNER_FILLET"))==0){
      DefineInnerFillet(line,cpt_grid,cpt_obj);
      if(Fromage.rotorindex==0){ cpt_mir++; cpt_obj=cpt_mir; }
      else{ cpt_rot++; cpt_obj=cpt_rot; }
    }
    if(strncmp(line,"ROTOR_CYLINDRICAL",strlen("ROTOR_CYLINDRICAL"))==0){
      if(sscanf(line,"ROTOR_CYLINDRICAL %lf %lf %lf %d", &Fromage.d, &Fromage.phi, &Fromage.offset, &Fromage.dir)!=4){
        cout<<"ERROR: Wrong number of parameters for ROTOR_CYLINDRICAL!"<<endl;
        exit(-18);
      }
      Fromage.phi = Fromage.phi*M_PI/180.;
      Fromage.rotorindex = 1;
      cpt_obj =0;
    }
    if(strncmp(line,"ROTOR_CARTESIAN",strlen("ROTOR_CARTESIAN"))==0){
      if(sscanf(line,"ROTOR_CARTESIAN %lf %lf %lf %d", &x, &y, &Fromage.offset, &Fromage.dir)!=4){
        cout<<"ERROR: Wrong number of parameters for ROTOR_CARTESIAN!"<<endl;
        exit(-19);
      }
      Fromage.d = sqrt(x*x+y*y);
      if(x<0){
        Fromage.phi = atan(y/x)+M_PI;
      }
      else if(x==0 && y>0){
        Fromage.phi = M_PI/2.;
      }
      else if(x==0 && y<0){
        Fromage.phi = -M_PI/2;
      }
      else if(x==0 && y==0){
        Fromage.phi = 0.;
      }
      else{
        Fromage.phi= atan(y/x);
      }
      Fromage.rotorindex = 1;
      cpt_obj =0;
    }
    if(strncmp(line,"G_GRAVITY",strlen("G_GRAVITY"))==0){
      if(sscanf(line,"G_GRAVITY %lf" , &Fromage.G)!=1){
        cout<<"ERROR: Wrong number of parameters for G_GRAVITY!"<<endl;
        exit(-20);
      } 
    }
    if(strncmp(line,"STEP",strlen("STEP"))==0){
      if(sscanf(line,"STEP %lf %d" , &Fromage.theta, &Fromage.nstep)!=2){
        cout<<"ERROR: Wrong number of parameters for STEP!"<<endl;
        exit(-21);
      }
      Fromage.theta = Fromage.theta*M_PI/180.;
      if(Fromage.nstep<=0){
        cout<<"ERROR: Number of rotor angle steps must be strictly positive!"<<endl;
        exit(-22);
      }
    }
    if(strncmp(line,"ARM_LENGTH",strlen("ARM_LENGTH"))==0){
      if(sscanf(line,"ARM_LENGTH %lf" , &Fromage.arm_length)!=1){
        cout<<"ERROR: Wrong number of parameters for ARM_LENGTH!"<<endl;
        exit(-23);
      }
      if(Fromage.arm_length<=0){
        cout<<"ERROR: The arm length must be strictly positive!"<<endl;
        exit(-24);
      }
    }
    if(strncmp(line,"BEAM_OFFSET",strlen("BEAM_OFFSET"))==0){
      if(sscanf(line,"BEAM_OFFSET %lf %lf" , &Fromage.beam_y, &Fromage.beam_z)!=2){
        cout<<"ERROR: Wrong number of parameters for BEAM_OFFSET!"<<endl;
        exit(-25);
      }
    }
    if(strncmp(line,"SIGNAL",strlen("SIGNAL"))==0){
      if(sscanf(line,"SIGNAL %d" , &Fromage.signal)!=1){
        cout<<"ERROR: Wrong number of parameters for SIGNAL!"<<endl;
        exit(-26);
      }
      if(Fromage.signal<=0){
        cout<<"ERROR: Fourier series order must be strictly positive!"<<endl;
        exit(-27);
      }
    }
    if(strncmp(line,"SAVE_DISPLACEMENT",strlen("SAVE_DISPLACEMENT"))==0){
      if(sscanf(line,"SAVE_DISPLACEMENT %s", Fromage.userfile0)!=1){
        cout<<"ERROR: Wrong number of parameters for SAVE_DISPLACEMENT!"<<endl;
        exit(-28);
      }
    }
    if(strncmp(line,"SAVE_MIRROR",strlen("SAVE_MIRROR"))==0){
      if(sscanf(line,"SAVE_MIRROR %s", Fromage.userfile1)!=1){
        cout<<"ERROR: Wrong number of parameters for SAVE_MIRROR!"<<endl;
        exit(-29);
      }
    }
    if(strncmp(line,"SAVE_ROTOR",strlen("SAVE_ROTOR"))==0){
      if(sscanf(line,"SAVE_ROTOR %s", Fromage.userfile2)!=1){
        cout<<"ERROR: Wrong number of parameters for SAVE_ROTOR!"<<endl;
        exit(-30);
      }
    }
    if(strncmp(line,"DEBUG",strlen("DEBUG"))==0){
      if(sscanf(line,"DEBUG %d", &Fromage.debug)!=1){
        cout<<"ERROR: Wrong number of parameters for DEBUG!"<<endl;
        exit(-31);
      }
      if(Fromage.debug<0 || Fromage.debug>4){
        cout<<"ERROR: The debug level can be either 0, 1, 2, 3 or 4."<<endl;
        exit(-32);
      }
    }
  }
  Fromage.nmir_obj = cpt_mir;
  Fromage.nrot_obj = cpt_rot;
  free(line);
  return;

}

/*-----------------------------*/

void GetRotorCoord(vector<vector<vector<double> > > rotor_init, double theta){

 double a=0, b=1, c=2, temp0, temp1, temp2, dtheta, sign0=1., sign2=1.;
  double phi = Fromage.phi;
  double tilt=0.;
  double off = 0.;
  if(Fromage.dir==0){
    a=1;
    b=0;
    c=2;
    sign0=-1.;
  }
  if(Fromage.dir==2){
    a=2;
    b=1;
    c=0;
    sign2=-1.;
  }
  for(int l=0; l<Fromage.nrot_obj; l++){
    for(int m=0; m<Fromage.nelements[l+Fromage.nmir_obj] ; m++){
      temp0 = rotor_init[l][m][0];
      temp1 = rotor_init[l][m][1];
      temp2 = rotor_init[l][m][2];
      Fromage.rotor[l][m][a]= sign0*temp0;
      Fromage.rotor[l][m][b]= temp1;
      Fromage.rotor[l][m][c]= sign2*temp2;

      if(Fromage.rotor[l][m][b]<0 && Fromage.rotor[l][m][c]>=0){
        dtheta = atan(Fromage.rotor[l][m][c]/Fromage.rotor[l][m][b])+M_PI;
      }
      else if(Fromage.rotor[l][m][b]<0 && Fromage.rotor[l][m][c]<0){
        dtheta = atan(Fromage.rotor[l][m][c]/Fromage.rotor[l][m][b])-M_PI;
      }
      else if(Fromage.rotor[l][m][b]==0 && Fromage.rotor[l][m][c]>0){
        dtheta = M_PI/2.;
      }
      else if(Fromage.rotor[l][m][b]==0 && Fromage.rotor[l][m][c]<0){
        dtheta = -M_PI/2.;
      }
      else if(Fromage.rotor[l][m][b]==0 && Fromage.rotor[l][m][c]==0){
        dtheta = 0.;
      }
      else{
        dtheta = atan(Fromage.rotor[l][m][c]/Fromage.rotor[l][m][b]);
      }

      Fromage.rotor[l][m][b] = sqrt(temp1*temp1+temp2*temp2)*cos(dtheta+theta);
      Fromage.rotor[l][m][c] = sqrt(temp1*temp1+temp2*temp2)*sin(dtheta+theta);

      temp0 = Fromage.rotor[l][m][0];
      temp1 = Fromage.rotor[l][m][1];
      temp2 = Fromage.rotor[l][m][2];

      Fromage.rotor[l][m][0] = cos(phi)*(cos(tilt)*temp0-sin(tilt)*temp1)-sin(phi)*(sin(tilt)*temp0+cos(tilt)*temp1)+Fromage.d*cos(phi);//cos(phi)*temp0-sin(phi)*temp1+Fromage.d*cos(phi);
      Fromage.rotor[l][m][1] = sin(phi)*(cos(tilt)*temp0-sin(tilt)*temp1)+cos(phi)*(sin(tilt)*temp0+cos(tilt)*temp1)+Fromage.d*sin(phi)+off;//sin(phi)*temp0+cos(phi)*temp1+Fromage.d*sin(phi);
      Fromage.rotor[l][m][2] = temp2+Fromage.offset;
    }
  }
  return;
}

/*-----------------------------*/

int main(int argc, char *argv[]){

  double element_dist, d3, dx, dy, dz, mass_mir, mass_rot, displacement, max_for=-1e10, min_for = 1e10, force_tot, tory_tot = 0., torz_tot = 0.;
  double max_tory=-1e10, min_tory = 1e10, max_torz=-1e10, min_torz = 1e10;
  vector<vector<vector<double> > > rotor_init;
  FILE *fsave_disp = NULL;
  FILE *fsave_mir = NULL;
  FILE *fsave_rot = NULL;
  char tempuser[256];
  cout<<endl;

  if(argc > 2){ 
    cout<<"ERROR: Number of arguments is invalid"<<endl;
    exit(-3);
  }

  if(Fromage.signal!=0){
    if(Fromage.nstep*Fromage.theta!=2*M_PI){
      cout<<"ERROR: Rotor period must be 360° to compute Fourier coefficients!"<<endl;
      exit(-2);
    }
  }

  cout<<"Reading cfg file..."<<endl;

  ReadCfg(argv[1]);
  vector<vector<double> > buffer(4,vector<double>(Fromage.nstep));

  if(strlen(Fromage.userfile0)!=0) {
    fsave_disp = fopen(Fromage.userfile0,"w+");
    fprintf(fsave_disp,"Theta rotor [deg] Mirror displacement [m.(n.f_rot)^2]\n\n");
  }
  if(strlen(Fromage.userfile1)!=0) {
    fsave_mir = fopen(Fromage.userfile1,"w+");
    fprintf(fsave_mir,"x [m] y [m] z [m] mass [kg] (in the mirror's frame)\n\n");
  }

  rotor_init = Fromage.rotor;

  for(int i=0;i<Fromage.nstep;i++){

    if(Fromage.debug>=1){
      if(i==0) cout<<"Computing NCal signal..."<<endl;
      cout<<"\033[K"<<i+1<<"/"<<Fromage.nstep<<endl;
      if(i<Fromage.nstep-1) cout<<"\033[F";
    }
    
    Fromage.force[0] = 0.;
    Fromage.force[1] = 0.;
    Fromage.force[2] = 0.;
    Fromage.torque[0] = 0.;
    Fromage.torque[1] = 0.;

    GetRotorCoord(rotor_init,i*Fromage.theta);

    if(strlen(Fromage.userfile2)!=0) {
      sprintf(tempuser,"%s_%d.txt",Fromage.userfile2,i);
      fsave_rot = fopen(tempuser,"w+");
      fprintf(fsave_rot,"x [m] y [m] z [m] mass [kg] (in the mirror's frame)\n\n");
    }

    for(int p=0; p<Fromage.nmir_obj;p++){ 
      for(int q=0; q<Fromage.nelements[p];q++){
        if(i==0) Fromage.mirror_mass += Fromage.mirror[p][q][3];

        for(int s=0; s<Fromage.nrot_obj;s++){
          for(int u=0; u<Fromage.nelements[s+Fromage.nmir_obj];u++){

	    dx = Fromage.rotor[s][u][0]-Fromage.mirror[p][q][0];
	    dy = Fromage.rotor[s][u][1]-Fromage.mirror[p][q][1];
	    dz = Fromage.rotor[s][u][2]-Fromage.mirror[p][q][2];
            mass_mir = Fromage.mirror[p][q][3];
            mass_rot = Fromage.rotor[s][u][3];

            element_dist = sqrt(dx*dx+dy*dy+dz*dz);
            d3 = element_dist*element_dist*element_dist;

  	    Fromage.force[0] += Fromage.G*mass_rot*mass_mir*dx/d3; 

	    if(Fromage.beam_y!=0 || Fromage.beam_z!=0){
  	      Fromage.force[1] += Fromage.G*mass_rot*mass_mir*dy/d3;
  	      Fromage.force[2] += Fromage.G*mass_rot*mass_mir*dz/d3;
              Fromage.torque[0] += Fromage.G*mass_rot*mass_mir*dx/d3*Fromage.mirror[p][q][2]-Fromage.G*mass_rot*mass_mir*dz/d3*Fromage.mirror[p][q][0];
              Fromage.torque[1] += Fromage.G*mass_rot*mass_mir*dy/d3*Fromage.mirror[p][q][0]-Fromage.G*mass_rot*mass_mir*dx/d3*Fromage.mirror[p][q][1];
	    }

            if(strlen(Fromage.userfile2)!=0){
	      if(p==0 && q==0) fprintf(fsave_rot,"%.5f %.5f %.5f %.5f \n", Fromage.rotor[s][u][0], Fromage.rotor[s][u][1], Fromage.rotor[s][u][2], Fromage.rotor[s][u][3]);
              if(p==0 && q==0 && s==Fromage.nrot_obj-1 && u==Fromage.nelements[Fromage.nrot_obj-1+Fromage.nmir_obj]-1){
                fclose(fsave_rot);
                cout<<"Rotor elements position successfully saved in "<<tempuser<<endl;
	      }
            }
	  } 
        }
        if(strlen(Fromage.userfile1)!=0){
          if(i==0) fprintf(fsave_mir,"%.5f %.5f %.5f %.5f \n", Fromage.mirror[p][q][0], Fromage.mirror[p][q][1], Fromage.mirror[p][q][2], Fromage.mirror[p][q][3]);
          if(i==0 && q==Fromage.nelements[Fromage.nmir_obj-1]-1 && p==Fromage.nmir_obj-1){
            fclose(fsave_mir);
            cout<<"Mirror elements position successfully saved in "<<Fromage.userfile1<<endl;
	  }
        }
      }
    }
    buffer[0][i] = i*Fromage.theta;
    buffer[1][i] = Fromage.force[0];
    buffer[2][i] = Fromage.torque[0];
    buffer[3][i] = Fromage.torque[1];
    displacement = Fromage.force[0]/(Fromage.mirror_mass*4*M_PI*M_PI);

    if(Fromage.debug>=2){
      cout<<"Theta rotor = "<<i*Fromage.theta*180./M_PI<<" deg --> Longitudinal force is "<<Fromage.force[0]<<" N_pk"<<endl;
      cout<<endl;
    }
    
    if(strlen(Fromage.userfile0)!=0){
      fprintf(fsave_disp,"%.3f    %e \n",i*Fromage.theta*180./M_PI, displacement);
      if(i==Fromage.nstep-1){
        fclose(fsave_disp);
        cout<<"Displacement successfully saved in "<<Fromage.userfile0<<endl;
      }
    }

    if(Fromage.beam_y!=0 || Fromage.beam_z!=0){
      if(Fromage.torque[0]>=max_tory) max_tory = Fromage.torque[0];
      if(Fromage.torque[0]<=min_tory) min_tory = Fromage.torque[0];
      if(Fromage.torque[1]>=max_torz) max_torz = Fromage.torque[1];
      if(Fromage.torque[1]<=min_torz) min_torz = Fromage.torque[1];
    }

    if(Fromage.force[0]>=max_for) max_for = Fromage.force[0];
    if(Fromage.force[0]<=min_for) min_for = Fromage.force[0];

  }

  force_tot = (max_for - min_for)/2.;

  if(Fromage.beam_y!=0 || Fromage.beam_z!=0){
    tory_tot = (max_tory - min_tory)/2.;
    torz_tot = (max_torz - min_torz)/2.;
  }

  cout<<endl;
  cout<<setprecision(40)<<"Amplitude of the longitudinal force is "<<force_tot<<" N_pk"<<endl;
  cout<<endl;

  if(Fromage.signal!=0){
    if(Fromage.beam_y!=0 || Fromage.beam_z!=0){ 
      cout<<"Amplitude of the torque around y is "<<tory_tot<<" N_pk.m"<<endl;
      cout<<"Amplitude of the torque around z is "<<torz_tot<<" N_pk.m"<<endl;
      cout<<endl;
      ComputeFourierTorque(buffer, force_tot, tory_tot, torz_tot);
      }
    else ComputeFourierForce(buffer, force_tot);
  }

  cout<<"Done!"<<endl;
  return 0;
}

