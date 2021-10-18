%***************************************************************************************************
%*   Calculate adiabatic flame temperature by presented code.
%*   I take no responsibilities for any errors in the code or damage thereby.
%*   Please notify me at zolfaghari1992iut@gmail.com if the code is used in any type of application.
%***************************************************************************************************
%*   Developer   : Ali Zolfaghari Sichani (18-10-2021)
%***************************************************************************************************
%*   References  : 
%*   https://web.mit.edu/16.unified/www/FALL/thermodynamics/notes/node111.html
%*   https://www.ohio.edu/mechanical/thermo/property_tables/combustion
%***************************************************************************************************
%*   Equation of complete combustion reaction   :   
%*   natural gas (methane + ethane + propane + butane) + air (O2 + N2) ----> CO2 + H2O + N2
%*   Inputs      :
%*   ExcessAir  (percent of excess air in combustion)  (% )
%*   X_METHANE  (percent of methane in natural gas  )  (% )
%*   X_ETHANE   (percent of ethane in natural gas   )  (% )
%*   X_PROPANE  (percent of propane in natural gas  )  (% )
%*   X_BUTANE   (percent of butane in natural gas   )  (% )
%*   Outputs      :
%*   T_res      (adiabatic flame temperature        )  (oC)
%***************************************************************************************************


clear,clc
close all
format compact
format long

load('DATA_LIB');

ExcessAir = 20.0;

X_METHANE = 93.90;
X_ETHANE = 3.60;
X_PROPANE = 1.20;
X_BUTANE = 1.30;




Hf_METHANE = -74850;
Hf_ETHANE = -84680;
Hf_PROPANE = -103850;
Hf_BUTANE = -126125;

Hf_CO2 = -393520;
Hf_H2O = -241820;

Hstd_CO2 = 0.0;
Hstd_H2O = 0.0;
Hstd_N2 = 0.0;
Hstd_O2 = 0.0;

M_H = 1.0;
M_O = 16.0;
M_C = 12.0;
M_N = 14.0;

X_CARBON = (X_METHANE+2.0*X_ETHANE+3.0*X_PROPANE+4.0*X_BUTANE)*0.01;
X_HYDROGEN = (4.0*X_METHANE+6.0*X_ETHANE+8.0*X_PROPANE+10.0*X_BUTANE)*0.01;

alpha = (1.0+0.01*ExcessAir)*(X_CARBON+0.25*X_HYDROGEN);

M_O2 = 2.0*M_O;
M_N2 = 2.0*M_N;
M_CO2 = M_C+2.0*M_O;
M_H2O = 2.0*M_H+M_O;

n_O2 = alpha-X_CARBON-0.25*X_HYDROGEN;
n_N2 = 3.76*alpha;
n_CO2 = X_CARBON;
n_H2O = 0.5*X_HYDROGEN;

HR = 0.01*(X_METHANE*Hf_METHANE+X_ETHANE*Hf_ETHANE+X_PROPANE*Hf_PROPANE+X_BUTANE*Hf_BUTANE);
HP = n_O2*(0.0-Hstd_O2)+n_H2O*(Hf_H2O-Hstd_H2O)+n_CO2*(Hf_CO2-Hstd_CO2)+n_N2*(0.0-Hstd_N2);
Q = HP-HR;

i = 1;
for T = 350.0:1.0:2900.0
    
    h_CO2 = pchip(TH_CO2(:,1),TH_CO2(:,2),T);
    h_H2O = pchip(TH_H2O(:,1),TH_H2O(:,2),T);
    h_O2 = pchip(TH_O2(:,1),TH_O2(:,2),T);
    h_N2 = pchip(TH_N2(:,1),TH_N2(:,2),T);
    
    H = n_O2*h_O2+n_H2O*h_H2O+n_CO2*h_CO2+n_N2*h_N2;
    
    T_ad(i) = T;
    Q_ad(i) = Q+H;

    i = i+1;
end

T_res = pchip(Q_ad,T_ad,0.0)-273.15;


fprintf(['  Adiabatic Temperature : %10.2f   (',char(0176),'C) \n'],T_res);




