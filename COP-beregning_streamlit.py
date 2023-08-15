import numpy as np
import streamlit as st
import plotly.express as px
import csv
import numpy as np
import matplotlib.pyplot as plt
from GHEtool import Borefield, FluidData, GroundData, PipeData
import pygfunction as gt


st.set_page_config(page_title="O store COP-beregning", page_icon="🔥")

with open("styles/main.css") as f:
    st.markdown("<style>{}</style>".format(f.read()), unsafe_allow_html=True)

st.title('O store COP-beregning 🤯')
st.markdown('Laget av Åsmund Fossum 👨🏼‍💻')

dekgrad = 0.90
virkgrad = 0.80  #Virkningsgrad til kompressor i varmepumpen
startCOP = st.number_input('Startgjett for COP',value = 3) #Startgjett for COP-verdi for hver time gjennom (hvert) år(et).

utetemp_for_maks_turtemp = -10
utetemp_for_min_turtemp = 0
maks_turtemp = 45
min_turtemp = 35

# COP-Verdier fra datablad
databladtemp35 = np.array([-5,-2,0,2,5,10,15])
COP_data35 = np.array([3.68, 4.03, 4.23, 4.41, 4.56, 5.04, 5.42])
databladtemp45 = np.array([-2,0,2,5,10,15])
COP_data45 = np.array([3.3, 3.47, 3.61, 3.77, 4.11, 4.4])

#Kjører lineær regresjon på COP-verdier fra datablad:
from lin_reg import *
lin_COP_data35 = lin_reg(databladtemp35,COP_data35)
lin_COP_data45 = lin_reg(databladtemp45,COP_data45)

# Plotter COP-verdier fra datablad

fig = plt.figure()
plt.plot(databladtemp35,COP_data35,'o')
plt.plot(databladtemp45,COP_data45,'o')
plt.plot(databladtemp35,lin_COP_data35)
plt.plot(databladtemp45,lin_COP_data45)
plt.legend(['Datablad 35 degC', 'Datablad 45 degC','Lineær 35 degC','Lineær 45 degC'])
plt.xlabel('Brønntemperatur (degC)')
plt.ylabel('COP')
plt.title('COP-verdier fra datablad, ved 100 % kapasitet', fontsize = 20)
st.pyplot(fig)

# Parametre i uttrykket for den lineære regresjonen av COP fra datablad
stigtall35 = (lin_COP_data35[-1]-lin_COP_data35[0])/(databladtemp35[-1]-databladtemp35[0])
konstledd35 = lin_COP_data35[-1]-stigtall35*databladtemp35[-1]

stigtall45 = (lin_COP_data45[-1]-lin_COP_data45[0])/(databladtemp45[-1]-databladtemp45[0])
konstledd45 = lin_COP_data45[-1]-stigtall45*databladtemp45[-1]

# Funksjon for Lineær interpolering:
def lin_interp(x,x1,x2,y1,y2):
    y = y1+(x-x1)*(y2-y1)/(x2-x1)
    return y

timer=np.linspace(1,8760,8760)

dato = []
varmelast = []
utetemp =[]
  
with open('timedata2.csv','r') as csvfile:
    lines = csv.reader(csvfile, delimiter='\t')
    for row in lines:
        dato.append(row[0])
        utetemp.append(float(row[4]))
        varmelast.append(float(row[7]))

grunnlast = 1*varmelast
spisslast = 1*varmelast

maks = np.max(varmelast)
mini = np.min(varmelast)

for kap in np.arange(0.8*maks,mini,-0.01*maks):  #Sjekker fra 80% av makslast og nedover med steglengde 1% av denne.
    for i in range(0, len(grunnlast)):
        if grunnlast[i]>=kap:
            grunnlast[i]=kap
        else:
            grunnlast[i]=grunnlast[i]
        
    for j in range(0,len(spisslast)):
            if grunnlast[j]>=kap:
                spisslast[j]=varmelast[j]-kap
            else:
                spisslast[j]=0

    if np.sum(grunnlast)/(np.sum(varmelast))<dekgrad:
        break
print(kap)

def GHE_tool_bronndybde(bronnlast,meter_init):
    #meter_init = 100

    data = GroundData(3.5, 7.5, 0.10, 2.4 * 10**6)

    borefield_gt = gt.boreholes.rectangle_field(N_1=1, N_2=1, B_1=6, B_2=6, H=meter_init, D = 10, r_b = 0.10)

    # create the borefield object

    borefield = Borefield(simulation_period=20)
    borefield.set_ground_parameters(data)
    borefield.set_borefield(borefield_gt)        
    borefield.set_hourly_heating_load(bronnlast)

    #st.write(borefield._check_hourly_load())

    borefield.set_max_ground_temperature(16)   # maximum temperature
    borefield.set_min_ground_temperature(0)    # minimum temperature
    meter = borefield.size(meter_init, L3_sizing=True)

    print('Nødvendig brønndybde (GHE-tool):',meter,'m.')
    return meter

grunnlast=np.array(grunnlast)

def kjor_pygf(bronnlast_data,bronndybde,antall_aar):
    # Beregning av brønntemperatur vha. Magnes pygfunction-kode:
    import pygfunction_Magne
    pygf = pygfunction_Magne.Simulation()                           # For at dette skal fungere, må det velges "Rektangulær" i pygfunction-filen
    pygf.select_borehole_field(1)               #Antall brønner 
    pygf.YEARS = antall_aar
    pygf.U_PIPE = "Single"  # Choose between "Single" and "Double"
    pygf.R_B = 0.114  # Radius (m)
    pygf.R_OUT = 0.020  # Pipe outer radius (m)
    pygf.R_IN = 0.0176  # Pipe inner radius (m)
    pygf.D_S = 0.067/2  # Shank spacing (m)
    pygf.EPSILON = 1.0e-6  # Pipe roughness (m)
    pygf.ALPHA = 1.39e-6  # Ground thermal diffusivity (m2/s)
    pygf.K_S = 3.5  # Ground thermal conductivity (W/m.K)            
    pygf.T_G = 7.5  # Undisturbed ground temperature (degrees)   
    pygf.K_G = 2  # Grout thermal conductivity (W/m.K)
    pygf.K_P = 0.42  # Pipe thermal conductivity (W/m.K)
    pygf.H = bronndybde  # Borehole depth (m)
    pygf.B = 15  # Distance between boreholes (m)
    pygf.D = 10  # Borehole buried depth
    pygf.FLOW_RATE = 0.5  # Flow rate (kg/s)
    pygf.FLUID_NAME = "MPG"  # The fluid is propylene-glycol 
    pygf.FLUID_DEGREES = 5  # at 20 degC
    pygf.BOUNDARY_CONDITION = 'MIFT'
    pygf.run_simulation(np.array(bronnlast_data)) #Grunnlast i enhet kW

    nybronntemp = pygf.tf_out
    return nybronntemp, pygf.tf_in


def finn_ny_COP(bronntemp_vektor):
    # Definerer COP-verdier for de faktiske brønntemperaturer og plotter disse.
    cop35 = stigtall35*bronntemp_vektor+konstledd35   #Her interpoleres det mellom maks. og min. verdier, i stedet for mellom de to nærmeste, som er en grov tilnærming.
    cop45 = stigtall45*bronntemp_vektor+konstledd45

    turtemp = np.zeros(len(utetemp))
    for i in range(0,len(utetemp)):
        if utetemp[i]<utetemp_for_maks_turtemp:
            turtemp[i] = maks_turtemp
        elif utetemp[i]>utetemp_for_min_turtemp:
            turtemp[i] = min_turtemp
        else:
            #Lineær interpolering:
            turtemp[i] = lin_interp(utetemp[i],utetemp_for_maks_turtemp,utetemp_for_min_turtemp,maks_turtemp,min_turtemp)

    # COP som funksjon av turtemp (basert på COP som funksjon av brønntemp)
    
    nyCOP = np.zeros(len(turtemp))
    for i in range(0,len(turtemp)):
        if turtemp[i] == maks_turtemp:
            nyCOP[i] = stigtall45*bronntemp_vektor[i]+konstledd45
            #np.append(nyCOP,stigtall45*bronntemp[i]+konstledd45)
        elif turtemp[i] == min_turtemp:
            nyCOP[i] = stigtall35*bronntemp_vektor[i]+konstledd35
            #np.append(nyCOP,stigtall35*bronntemp[i]+konstledd35)
        else:
            stigtall_interp = lin_interp(turtemp[i],35,45,stigtall35,stigtall45)
            konstledd_interp = lin_interp(turtemp[i],35,45,konstledd35,konstledd45)
            COP_interp = stigtall_interp*bronntemp_vektor[i]+konstledd_interp
            nyCOP[i] = COP_interp
            #np.append(nyCOP,COP_interp)
        nyCOP=np.array(nyCOP)
        #if np.sum(np.absolute(nyCOP-COP))/len(COP)<0.01:  #Stopper dersom gjennomsnittlig endring i COP-verdi fra forrige iterasjon er mindre enn 0.01
        #    break
        COP=nyCOP
    return cop35, cop45, COP

def finn_ny_COP_konst_turtemp(bronntemp_vektor):
    nyCOP = (stigtall35+stigtall45)/2 * bronntemp_vektor + (konstledd35+konstledd45)/2
    return nyCOP

ukjent_param = 'Brønndybde'  # Eller 'COP' eller 'Brønntemp'

antall_aar = 25
grunnlast = np.hstack(antall_aar*[grunnlast])
utetemp = np.hstack(antall_aar*[utetemp])

startCOP = np.array([startCOP]*8760*antall_aar)
COP = startCOP
bronndybde = 50

ellast = grunnlast/COP*virkgrad
bronnlast = grunnlast-ellast

if ukjent_param == 'Validering Brønndybde':
    for k in range(0,10):
        [nybronntemp, t_inn] = kjor_pygf(bronnlast,bronndybde,antall_aar)

        print(nybronntemp)

        if np.min(nybronntemp) <= 0:
            bronndybde = bronndybde+0.25*bronndybde
        elif np.min(nybronntemp) >= 2:
            bronndybde = bronndybde-0.5*bronndybde
        else:
            break

    print('Antall iterasjoner:',k+1)
    print('Nødvendig brønndybde:',bronndybde)
    print('Minste brønntemperatur:',np.min(nybronntemp))
    print(COP)

    bronntemp = nybronntemp
    dybde_GHE = GHE_tool_bronndybde(bronnlast,bronndybde)



elif ukjent_param == 'Brønndybde':
    dybde_GHE = GHE_tool_bronndybde(bronnlast,bronndybde)
    #dybde_GHE = 108
    
    for k in range(0,5):
        [T_ut_her,T_inn] = kjor_pygf(bronnlast,dybde_GHE,antall_aar)
        print('Brønntemp ut fra pygf:')
        print(T_ut_her)
        print('Minste Bronntemp:',np.min(T_ut_her))
        #[COP_35,COP_45,COP_faktisk] = finn_ny_COP(T_ut_her)
        COP_faktisk = finn_ny_COP_konst_turtemp(T_ut_her)
        print('Ny COP basert på brønntemp fra pygf:')
        print(COP_faktisk)
        ellast = grunnlast/COP_faktisk*virkgrad
        
        if np.mean((np.abs(bronnlast - (grunnlast-ellast)))) <= 0.01:
            break
        else:
            bronnlast = grunnlast-ellast

