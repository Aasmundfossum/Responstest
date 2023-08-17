import numpy as np
import streamlit as st
import plotly.express as px
import csv
import numpy as np
import matplotlib.pyplot as plt
from GHEtool import Borefield, FluidData, GroundData, PipeData
import pygfunction as gt
import math
import itertools
from lin_reg import *


st.set_page_config(page_title="O store COP-beregning", page_icon="🔥")

with open("styles/main.css") as f:
    st.markdown("<style>{}</style>".format(f.read()), unsafe_allow_html=True)


## ------------------------------------------------------------------------------------------------------##

@st.cache_data
def plot_datablad():
    # COP-Verdier fra datablad
    databladtemp35 = np.array([-5,-2,0,2,5,10,15])
    COP_data35 = np.array([3.68, 4.03, 4.23, 4.41, 4.56, 5.04, 5.42])
    databladtemp45 = np.array([-2,0,2,5,10,15])
    COP_data45 = np.array([3.3, 3.47, 3.61, 3.77, 4.11, 4.4])

    #Kjører lineær regresjon på COP-verdier fra datablad:
    
    lin_COP_data35 = lin_reg(databladtemp35,COP_data35)
    lin_COP_data45 = lin_reg(databladtemp45,COP_data45)

    # Plotter COP-verdier fra datablad

    fig = plt.figure()
    plt.plot(databladtemp35,COP_data35,'o')
    plt.plot(databladtemp45,COP_data45,'o')
    plt.plot(databladtemp35,lin_COP_data35)
    plt.plot(databladtemp45,lin_COP_data45)
    plt.legend(['Datablad 35 \u2103', 'Datablad 45 \u2103','Lineær 35 \u2103','Lineær 45 \u2103'])
    plt.xlabel('Brønntemperatur (\u2103)')
    plt.ylabel('COP')
    plt.title('COP-verdier fra datablad, ved 100 % kapasitet', fontsize = 20)
    st.pyplot(fig)

    # Parametre i uttrykket for den lineære regresjonen av COP fra datablad
    stigtall35 = (lin_COP_data35[-1]-lin_COP_data35[0])/(databladtemp35[-1]-databladtemp35[0])
    konstledd35 = lin_COP_data35[-1]-stigtall35*databladtemp35[-1]

    stigtall45 = (lin_COP_data45[-1]-lin_COP_data45[0])/(databladtemp45[-1]-databladtemp45[0])
    konstledd45 = lin_COP_data45[-1]-stigtall45*databladtemp45[-1]

    return stigtall35,konstledd35,stigtall45,konstledd45


# Funksjon for Lineær interpolering:
def lin_interp(x,x1,x2,y1,y2):
    y = y1+(x-x1)*(y2-y1)/(x2-x1)
    return y

def bronnlast_fra_COP(grunnlast,cop,virkgrad):
    ellast = grunnlast/cop*virkgrad
    bronnlast = grunnlast-ellast

    return bronnlast


def GHE_tool_bronndybde(bronnlast,min_bronntemp,dybde_startgjett,ledningsevne,uforst_temp,term_motstand,antall_aar,bronnfelt):
    data = GroundData(ledningsevne, uforst_temp, term_motstand, 2.518 * 10**6)    # Siste parameter: Volumetric heat capacity of ground
    
    #borefield_gt = gt.boreholes.rectangle_field(N_1=ant_bronner1, N_2=ant_bronner2, B_1=10, B_2=10, H=dybde_startgjett, D = 10, r_b = 0.114) # Siste to parametre: Boreholde buried depth og borehole radius (m)
    borefield_gt = bronnfelt

    borefield = Borefield(simulation_period=antall_aar)
    borefield.set_ground_parameters(data)
    borefield.set_borefield(borefield_gt)        
    #borefield.set_hourly_heating_load(bronnlast)

    borefield.hourly_heating_load = bronnlast[-8760:]
    borefield.hourly_cooling_load = np.zeros(8760)

    borefield.set_max_ground_temperature(16)   # maximum temperature   Utgjør ingen forskjell å endre på denne.
    borefield.set_min_ground_temperature(min_bronntemp)    # minimum temperature
    dybde = borefield.size(dybde_startgjett, L4_sizing=True)
    
    snitt_koll_vaeske_temp = borefield.results_peak_heating
    bronntemp_vegg = borefield.Tb
    #SETT INN DELTA T BEREGNING
    bronntemp_tur = snitt_koll_vaeske_temp + 1.5
    bronntemp_retur = snitt_koll_vaeske_temp - 1.5


    print('Nødvendig brønndybde (GHE-tool):',dybde,'m.')
    return dybde,bronntemp_vegg,bronntemp_tur,bronntemp_retur


def kjor_pygf(bronnlast_data,bronndybde,antall_aar,ledningsevne,uforst_temp):
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
    pygf.K_S = ledningsevne  # Ground thermal conductivity (W/m.K)            
    pygf.T_G = uforst_temp  # Undisturbed ground temperature (degrees)   
    pygf.K_G = 2  # Grout thermal conductivity (W/m.K)
    pygf.K_P = 0.42  # Pipe thermal conductivity (W/m.K)
    pygf.H = bronndybde  # Borehole depth (m)
    pygf.B = 15  # Distance between boreholes (m)
    pygf.D = 10  # Borehole buried depth
    pygf.FLOW_RATE = 0.5  # Flow rate (kg/s)
    pygf.FLUID_NAME = "MPG"  # The fluid is propylene-glycol 
    pygf.FLUID_DEGREES = 5  # at 20 \u2103
    pygf.BOUNDARY_CONDITION = 'MIFT'
    pygf.run_simulation(np.array(bronnlast_data)) #Grunnlast i enhet kW

    print(pygf.R_B)
    nybronntemp = pygf.tf_out
    returtemp = pygf.tf_in
    term_motstand = pygf.R_B
    return nybronntemp,returtemp,term_motstand


def bestem_turtemp(utetemp,utetemp_for_maks_turtemp,utetemp_for_min_turtemp,maks_turtemp,min_turtemp):
    turtemp = np.zeros(len(utetemp))
    for i in range(0,len(utetemp)):
        if utetemp[i]<utetemp_for_maks_turtemp:
            turtemp[i] = maks_turtemp
        elif utetemp[i]>utetemp_for_min_turtemp:
            turtemp[i] = min_turtemp
        else:
            #Lineær interpolering:
            turtemp[i] = lin_interp(utetemp[i],utetemp_for_maks_turtemp,utetemp_for_min_turtemp,maks_turtemp,min_turtemp)
    return turtemp


def finn_ny_COP(bronntemp_vektor,stigtall35,konstledd35,stigtall45,konstledd45,turtemp,maks_turtemp,min_turtemp):
    # COP som funksjon av turtemp (basert på COP som funksjon av brønntemp)
    nyCOP = np.zeros(len(turtemp))
    for i in range(0,len(turtemp)):
        if turtemp[i] == maks_turtemp:
            nyCOP[i] = stigtall45*bronntemp_vektor[i]+konstledd45
        elif turtemp[i] == min_turtemp:
            nyCOP[i] = stigtall35*bronntemp_vektor[i]+konstledd35
        else:
            stigtall_interp = lin_interp(turtemp[i],min_turtemp,maks_turtemp,stigtall35,stigtall45)
            konstledd_interp = lin_interp(turtemp[i],min_turtemp,maks_turtemp,konstledd35,konstledd45)
            COP_interp = stigtall_interp*bronntemp_vektor[i]+konstledd_interp
            nyCOP[i] = COP_interp
    nyCOP=np.array(nyCOP)

    COP=nyCOP
    return COP


## ------------------------------------------------------------------------------------------------------##

class O_store_COP_beregning:
    def __init__(self):
        pass
    
    def kjor_hele(self):
        self.streamlit_input()
        if self.kjor_knapp == True:
            self.last_inn_varmebehov()
            self.grunnlast_fra_varmelast()
            self.dybde_COP_sloyfe()


    def streamlit_input(self): 

        st.title('O store COP-beregning 🤯')
        st.markdown('Laget av Åsmund Fossum 👨🏼‍💻')
        st.markdown('---')

        st.subheader('Data om grunnen')
        c1, c2 = st.columns(2)
        with c1:
            self.LEDNINGSEVNE = st.number_input('Termisk ledningsevne til grunnen (W/mK)',value=3.5, max_value=float(10), min_value=0.1, step=0.1)
        with c2:
            self.UFORST_TEMP = st.number_input('Uforstyrret temperatur (\u2103)',value=7.5, step=0.5)
        st.markdown('---')


        st.subheader('Brønnkonfigurasjon')
        c1, c2 = st.columns(2)
        with c1:
            self.MIN_BRONNTEMP = st.number_input('Laveste tillatte gjennomsnittlige kollektorvæsketemperatur etter 25 år (\u2103)', value = float(0), step=0.1)
        with c2:
            self.MAKS_DYBDE = st.number_input('.Maksimal tillatte dybde per brønn (m)', value=300, step=10)

        self.onsket_konfig = st.selectbox(label='Foretrukket brønnkonfigurasjon', options=['Linje','Kvadrat (tilnærmet)','Rektangel med fastsatt side'], index=1)
        
        c1, c2 = st.columns(2)
        with c1:
            self.avstand = st.number_input('Avstand mellom nærliggende brønner (m)', value=15, min_value=1, step=1)
        with c2: 
            if self.onsket_konfig == 'Rektangel med fastsatt side':
                self.fastsatt_side = st.number_input('Antall brønner langs en side av rektangelet', value=2, min_value=2, step=1)
        st.markdown('---')

        st.subheader('Valg av varmepumpe')
        c1, c2, = st.columns(2)
        with c1:
            type_VP = st.selectbox(label='Velg varmepumpe',options=['Eksempel-VP', 'VP 2', 'VP 3', 'VP 4'])
        with c2:
            self.VIRKGRAD = st.number_input('Virkningsgrad til kompressor (%)',value=80, max_value=100, min_value=10, step=1)/100

        if type_VP == 'Eksempel-VP':
            [self.stigtall35,self.konstledd35,self.stigtall45,self.konstledd45] = plot_datablad()
            self.MAKS_TURTEMP = 45
            self.MIN_TURTEMP = 35
        elif type_VP == 'VP 2':
            pass

        
        c1, c2 = st.columns(2)
        with c1:
            st.metric("Maksimal turtemperatur", f"{self.MAKS_TURTEMP} \u2103")
        with c2:
            st.metric("Minimal turtemperatur", f"{self.MIN_TURTEMP} \u2103")
        
        c1,c2 = st.columns(2)
        with c1:
            self.UTETEMP_FOR_MAKS_TURTEMP = st.number_input('Høyeste utetemperatur med maksimal turtemperatur', value=-15)
        with c2:
            self.UTETEMP_FOR_MIN_TURTEMP = st.number_input('Laveste utetemperatur med minimal turtemperatur', value=15)
        st.markdown('---')

        st.subheader('Varmebehov')
        self.varmelast_fil = st.file_uploader('CSV-fil med timesoppløst varmeenergibehov for et år',type='csv')
        self.DEKGRAD = st.number_input('Ønsket dekningsgrad (%)',value=90, max_value=100, min_value=1, step=1)/100

        self.kjor_knapp = st.checkbox('Kjør :steam_locomotive:')
        
        st.markdown('---')


        self.ANTALL_AAR = 25
        self.COP = 3.5                       # Resultat uavhengig av denne
        self.DYBDE_STARTGJETT = 250          # Resultat uavhengig av denne
        self.TERM_MOTSTAND = 0.08

    def last_inn_varmebehov(self):
        @st.cache_data
        def last_inn_energibehov_csv(filnavn):
            dato = []
            varmelast = []
            utetemp =[]
            
            with open(filnavn,'r') as csvfile:
                lines = csv.reader(csvfile, delimiter='\t')
                for row in lines:
                    dato.append(row[0])
                    utetemp.append(float(row[4]))
                    varmelast.append(float(row[7]))

            return dato,utetemp,varmelast

        [self.dato,self.utetemp,self.varmelast] = last_inn_energibehov_csv(self.varmelast_fil.name)

    def grunnlast_fra_varmelast(self):
        @st.cache_data
        def grunnlast_fra_varmelast(varmelast,DEKGRAD,ANTALL_AAR):
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

                if np.sum(grunnlast)/(np.sum(varmelast))<DEKGRAD:
                    break

            grunnlast = np.array(grunnlast)*20
            GRUNNLAST = np.hstack(ANTALL_AAR*[grunnlast])
            return GRUNNLAST

        self.GRUNNLAST = grunnlast_fra_varmelast(self.varmelast,self.DEKGRAD,self.ANTALL_AAR)

    def dybde_COP_sloyfe(self):
        
        turtemp = bestem_turtemp(self.utetemp,self.UTETEMP_FOR_MAKS_TURTEMP,self.UTETEMP_FOR_MIN_TURTEMP,self.MAKS_TURTEMP,self.MIN_TURTEMP)
        TURTEMP = np.hstack(self.ANTALL_AAR*[turtemp])

        COP = np.array([self.COP]*8760*self.ANTALL_AAR)

        dybde_GHE = self.DYBDE_STARTGJETT
        ant_bronner1 = 1
        
        if self.onsket_konfig == 'Rektangel med fastsatt side':
            ant_bronner2 = self.fastsatt_side
        else:
            ant_bronner2 = 1

        for k in range(0,20):

            BRONNLAST = bronnlast_fra_COP(self.GRUNNLAST,COP,self.VIRKGRAD)

            bronnfelt = gt.boreholes.rectangle_field(N_1=ant_bronner1, N_2=ant_bronner2, B_1=self.avstand, B_2=self.avstand, H=dybde_GHE, D = 10, r_b = 0.114) # Siste to parametre: Boreholde buried depth og borehole radius (m)

            [dybde_GHE,bronntemp_vegg,bronntemp_tur,bronntemp_retur] = GHE_tool_bronndybde(BRONNLAST,self.MIN_BRONNTEMP,dybde_GHE,self.LEDNINGSEVNE,self.UFORST_TEMP,self.TERM_MOTSTAND,self.ANTALL_AAR,bronnfelt)

            if k==0 and dybde_GHE >= self.MAKS_DYBDE:
                ant_bronner1 = math.ceil(dybde_GHE/self.MAKS_DYBDE) #runder alltid opp
                print(ant_bronner1)

            print('DYBDE:',dybde_GHE)

            if k > 0 and dybde_GHE >= self.MAKS_DYBDE:
                ant_bronner1 = ant_bronner1+1
            
            print(ant_bronner1)
            print(ant_bronner2)
            
            #[bronntemp,returtemp,TERM_MOTSTAND] = kjor_pygf(BRONNLAST,dybde_GHE,ANTALL_AAR,LEDNINGSEVNE,UFORST_TEMP)

            nyCOP = finn_ny_COP(bronntemp_tur,self.stigtall35,self.konstledd35,self.stigtall45,self.konstledd45,TURTEMP,self.MAKS_TURTEMP,self.MIN_TURTEMP)
            
            if np.mean(np.abs(nyCOP-COP))<0.001 and dybde_GHE <= self.MAKS_DYBDE: 
                COP = nyCOP
                if self.onsket_konfig == 'Kvadrat (tilnærmet)' and ant_bronner2 == 1:
                    ant_per_side = math.ceil(np.sqrt(ant_bronner1))
                    ant_bronner1 = ant_per_side
                    ant_bronner2 = ant_per_side
                
                elif self.onsket_konfig == 'Rektangel med fastsatt side':
                    break
                elif self.onsket_konfig == 'Sirkel':
                    break
                else:  # Hvis ønkset konfig. er linje
                    break
            else:
                COP = nyCOP


        c1, c2, c3 = st.columns(3)
        with c1:
            st.metric("Brønnkonfigurasjon", f"{ant_bronner1} x {ant_bronner2}")
        with c2:
            st.metric("Totalt antall brønner", f"{ant_bronner1*ant_bronner2}")
        with c3:
            st.metric("Dybden til hver brønn", f"{round(dybde_GHE,1)} m")

        
        r = ant_bronner2
        c = ant_bronner1
        x = [i % c * self.avstand for i in range(r*c)]
        y = [i // c * self.avstand for i in range(r*c)]
        fig2, ax2 = plt.subplots()
        ax2.scatter(x,y,color = '#367A2F')

        x_ticks = [i * self.avstand for i in range(c)]
        y_ticks = [i * self.avstand for i in range(r)]
        ax2.set_xticks(x_ticks)
        ax2.set_yticks(y_ticks)
        ax2.grid(True)
        ax2.set_title('Brønnkonfigurasjon')
        ax2.set_xlabel('Avstand (m)')
        ax2.set_ylabel('Avstand (m)')
        ax2.axis('equal')
        st.pyplot(fig2)


        c1, c2, c3 = st.columns(3)
        with c1:
            st.metric("Laveste turtemp. (fra brønn)", f"{round(np.min(bronntemp_tur),2)} \u2103")
        with c2:
            st.metric("Laveste returtemp. (til brønn)", f"{round(np.min(bronntemp_retur),2)} \u2103")
        with c3:
            st.metric("Laveste veggtemp. i brønnen", f"{round(np.min(bronntemp_vegg),2)} \u2103")

        snittCOP = round(np.mean(COP),2)
        c1, c2, c3 = st.columns(3)
        with c1:
            st.metric("Maksimal COP", f"{round(np.max(COP),2)}")
        with c2:
            st.metric("Minimal COP", f"{round(np.min(COP),2)}")
        with c3:
            st.metric("Gjennomsnittlig COP", f"{snittCOP}")
            


        
        #print('')
        #print('---------------------------------------------')
        #print('Brønndybde: \t \t',round(dybde_GHE,3),'m')
        #print('Minste turtemp: \t',round(np.min(bronntemp_tur),3),'\u2103')
        #print('Minste returtemp: \t',round(np.min(bronntemp_retur),3),'\u2103')
        #print('Minste veggtemp: \t',round(np.min(bronntemp_vegg),3),'\u2103')
        #print('Gjennomsnittlig COP: \t',snittCOP)
        #print('Laveste COP-verdi: \t',round(np.min(COP),3))
        #print('')
        #print('Brønnlast gitt denne COP:')
        #print(BRONNLAST)
        #print('')
        #print(bronntemp_tur)
        #print(bronntemp_retur)




O_store_COP_beregning().kjor_hele()