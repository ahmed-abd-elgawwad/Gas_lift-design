from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *
import matplotlib
from PyQt5.uic import loadUiType
import sys
from os import path
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.pyplot  import Figure
import matplotlib.pyplot as plt
from pytube import Playlist
from pytube import YouTube
import pandas as pd
import numpy as np
import os
matplotlib.style.use("seaborn")
 #------------------------- functions used ----------------------------------------------
def get_equations(Pcs,j,BHSP,Q,Gs,D_total):
    """
    Pcs : surface injection pressure (psi) 
    j : productivity index 
    BHSP : bottomhole static pressure (psi)
    Q : oil flow rate ( bbl / day)
    Gs : static gradient ( psi/ft )
    D_total : total depth to the perforation (Ft)
    return 
    1- equation of reservoir line
    2- equation of caisng line
    """
        # so that we get injection point
    Gcs = Pcs/40000
    Pcs = Pcs-100
    Pwf= BHSP-(Q/j)
    # reservoir line
    reservoir_line = [1 ,-Gs, (Pwf - (Gs* D_total))]
    # casing line
    casing_line = [1,-Gcs,Pcs]
    
    return reservoir_line, casing_line ,Pwf

  # find the intersections
  
def find_intersection(line_1,line_2):
    """
    line_1 : coef of line_1 , in the form (a*P + b*Depth = c) , put list of [a1,b1,c1]
    line_2 : coef of line_2 , in the form (a*P + b*Depth = c) , put list of [a2,b2,c2]
    return : X,Y of the intersection
    """
    Z= np.array([ [line_1[0],line_1[1]],[line_2[0],line_2[1]] ])
    X= np.array([ [line_1[-1],line_1[-2]],[line_2[-1],line_2[-2]] ])
    Y= np.array([ [line_1[0],line_1[-1]],[line_2[0],line_2[-1]] ])
    try:
        D = np.linalg.det(Z)
        Dx= np.linalg.det(X)
        Dy= np.linalg.det(Y)
        P= Dx/D
        Depth= Dy/D 
        return np.round(Depth,2) , np.round(P,2)
    except:
        print("Error, there is no intersection")
        
    # get the gradeints 
    
def get_GF1_GF2(Tubing_head,Pcs,injection_depth, injection_pressure):
    """
    Tubing_head : tubing head pressure (Psi)
    Pcs : surface casing pressure (Psi)
    injection_depth : depth at the pint of injection (Ft)
    injection_pressure : pressure of the reservoir at the injection point (Psi)
    return : 
        Gf1 -> of the line from the wellhead to injection point
        Gf2 -> of the line from the modified wellhead to injection point
        Pwh2 -> modified wellhead pressure
    """
    Gf1 = (injection_pressure - Tubing_head) / injection_depth
    Pwh2= Tubing_head + 0.2*(Pcs - Tubing_head)
    Gf2=  (injection_pressure - Pwh2) / injection_depth
    return round(Gf1,3) , round(Gf2,3), round(Pwh2,3) 
 
    # the spacsing 
def get_spacings(Pwh1,Pwh2,Pcs,Pko,Glf,Gf2,injection_depth,injection_pressure):
    depths=[]
    P1=[]
    P2=[]  
    casing_ko=[1,(-Pko/40000),Pko]
    line_1= [1,-Glf,Pwh1]
    depth,p1= find_intersection(casing_ko,line_1)
    depths.append(depth)
    P1.append(p1)
    p2= Pwh2 + Gf2*depths[-1]
    P2.append(p2)
    # new_line
    Gcs = Pcs/40000
#     for i in range(30):
    while depths[-1] < injection_depth:
        Pc_0 = Pcs*(1+(depths[-1]/40000))
#         casing_ps=[1,-Pc_0/40000,Pc_0]
        casing_ps=[1,-Gcs,Pc_0]
        line_2= [1,-Glf,P2[-1]]
        depth,p1= find_intersection(casing_ps,line_2)
        depths.append(depth + depths[-1])
        P1.append(p1)
        p2= Pwh2 + Gf2*depths[-1]
        P2.append(p2)
    dics={
    "Depth":np.round(depths,2),
     "P1":np.round(P1,2),
     "P2":np.round(P2,2)
     }

    df =pd.DataFrame(dics)
    df.drop(axis=0,index=df.index[-1],inplace=True)
    final_depth=np.round(injection_depth,2)
    p1= Pcs * (1+ (final_depth / 40000))
    p2= injection_pressure
    final={
    "Depth":np.round(injection_depth,2),
     "P1":np.round(p1,2),
     "P2":np.round(p2,2)
     }
    
    df=df.append(final,ignore_index=True)
    # new lines 
    if (df.loc[df.index[-1],"Depth"] - df.loc[df.index[-2],"Depth"])  < 200:
        df.drop(axis=0,index=df.index[-2],inplace=True)
        
    # write the number of valves
    valves= [i+1 for i in range(len(df["Depth"]))]
    df["Valve .No"]= valves
    # df=df[["Valve .No","Depth","P1","P2"]]
    cols=["Valve .No","Depth","P1","P2"]
    df=df.reindex(columns=cols)
    return df 

 # all in one 
def all_in_one(total_depth,wellhead_pressure, Pcs ,Pko,Glf,Gs,Q,BHSP,J,Tre,Ts,R):
    """
    total_depth : total depth to the pay zone  ( ft )
    wellhead_pressure : tuning pressure at the surface ( psi )
    Pcs : surface casing pressure ( injection pressure) ( psi )
    Pko : kick_off pressure if it exists ( psi )
    Glf : load flowing gradient ( psi/ft )
    Gs : static gradient ( psi/ft )
    Q : oil flow rate ( bbl/day)
    BHSP : bottom_hole static pressure -> reservoir pressure
    J : productivity index 
    Tres : reservoir temperature ( F )
    Ts : surface following temperature ( F )
    R: port size 
    """
    reservoir_line,casing_line,Pwf = get_equations(Pcs,J,BHSP,Q,Gs,total_depth)
    injection_depth,injection_pressure = find_intersection(reservoir_line,casing_line)
    Gf1,Gf2,Pwh2 = get_GF1_GF2(wellhead_pressure,Pcs,injection_depth,injection_pressure)
    df=get_spacings(wellhead_pressure,Pwh2,Pcs,Pko,Glf,Gf2,injection_depth,injection_pressure) 
    
    # Tg= ((Tre-Ts) / total_depth)
    # df["Temp"] = np.round(Ts + Tg * df["Depth"],2)
    # df["Pdt"] = np.round((1-R)*df["P1"] + R*df["P2"],2)
    # make the Ct and get Pd and also Pvo
    Ct_data= pd.read_csv("Ct_Data.csv")
    Tg= ((Tre-Ts) / total_depth)
    df["Temp"] = np.round(Ts + Tg * df["Depth"],2)
    df["Pdt"] = np.round((1-R)*df["P1"] + R*df["P2"],2)
    df["Ct"]= df.Temp.apply(lambda T: (Ct_data[(Ct_data["temp"]== round(T))]["Ct"].values)[0] ) 
    df["Pd"]=round(df["Pdt"] * df["Ct"],2)
    df["Pvo"] = round(df["Pd"] / (1-R) ,2)
    return df,Pwf,wellhead_pressure,Pwh2,Pcs,Pko,total_depth,reservoir_line
#---------------------------------------- end functions ----------------------------------    

#----------------------------------------- canvas drawign --------------------------------

class MatplotlibWidget(QWidget):
     def __init__(self,parent =None):
            super(MatplotlibWidget,self).__init__(parent)
            self.figure = Figure()
            self.canvas = FigureCanvasQTAgg(self.figure)
            self.axis=self.figure.add_subplot(111)
            self.layoutvertical= QVBoxLayout(self)
            self.layoutvertical.addWidget(self.canvas)
            # initiate the grapfh
            self.figure.gca().invert_yaxis()
            self.axis.xaxis.tick_top()
            self.axis.xaxis.set_label_position("top")
            self.axis.set_xlabel("Pressure (psi)")
            self.axis.set_ylabel("Depth (ft)")
# ---------------------------------------- end canvas drawing ---------------------------- 
            
FORM_CLASS,_ =loadUiType(path.join(path.dirname(__file__),"final.ui"))
class MainApp(QMainWindow,FORM_CLASS):
    def __init__(self,parent =None):
        super(MainApp,self).__init__(parent)
        QMainWindow.__init__(self)
        self.setupUi(self)
        self.ui_settings()
        self.buttons()
        self.init_widget()
        self.dataframe =0
        # self.init_widget()
    # handling the ui settings
    def ui_settings(self):
        # self.setFixedSize(1130,830)
        self.setWindowTitle("GasLift")
    # init the drawing area 
    def init_widget(self):
        self.mat = MatplotlibWidget()
        self.layoutvertical= QVBoxLayout(self.widget)
        self.layoutvertical.addWidget(self.mat)
    # hindeling button and connecting them
    def buttons(self):
        self.pushButton.clicked.connect(self.save_file)
        self.pushButton_2.clicked.connect(self.get_location)
        self.pushButton_3.clicked.connect(self.make_design)
    def get_location(self):
        self.save_location=QFileDialog.getExistingDirectory(self,"select location to save file")
        self.lineEdit_6.setText(self.save_location)
    def save_file(self):
        file_name= self.lineEdit_10.text()
        if (file_name != None) & ( self.dataframe.empty == False):
            file_path=f"{self.save_location}/{file_name}.xlsx"
            with pd.ExcelWriter(file_path) as file:
                    self.dataframe.to_excel(file,sheet_name="GasLift_design",index=False)
            QMessageBox.information(self,"sucess","The file is made ,The file will open.")
            os.system(f'start "excel" "{file_path}"')
        else:
            QMessageBox.warning(self,"warning","please,file the file name and the required data")

            
     # initiate the canvas for the drawing 
     
    # def init_widget(self):
    #     self.matplotlib_w = mat()
    #     # self.layoutvertical =QBoxLayout(self)
    #     # self.layoutvertical.addWidget(self.matplotlib_w)
        
    def make_design(self):
        arguments=[]
        arguments.append( self.lineEdit.text() )# total_depth
        arguments.append(self.lineEdit_3.text() )# Pwf
        arguments.append(self.lineEdit_2.text() )# Pcs
        arguments.append(self.lineEdit_4.text()) # Pko
        arguments.append( self.lineEdit_17.text()) # Glf
        arguments.append( self.lineEdit_18.text()) # Gs 
        arguments.append(self.lineEdit_19.text() )# Q
        arguments.append( self.lineEdit_20.text()) # BHSP
        arguments.append( self.lineEdit_5.text() )# J
        arguments.append( self.lineEdit_7.text() )# Tres
        arguments.append( self.lineEdit_8.text() )# Ts
        arguments.append( self.lineEdit_9.text() )# R
        if len([x for x in arguments if x!=""]) <12 :
            QMessageBox.warning(self,"warning","please enter all the variables")
        else:
            # try:
                arguments= [float(x) for x in arguments]
                df,Pwf,wellhead_pressure,Pwh2,Pcs,Pko,total_depth,reservoir_line=all_in_one( arguments[0],arguments[1],arguments[2], arguments[3], arguments[4], arguments[5], arguments[6], arguments[7], arguments[8], arguments[9], arguments[10], arguments[11])
                # file_name= self.lineEdit_10.text()
                # try:
                # if (file_name != None) & ( df.empty == False):
                    # file_path=f"{self.save_location}/{file_name}.xlsx" 
                    # with pd.ExcelWriter(file_path) as file:
                    #     df.to_excel(file,sheet_name="GasLift_design",index=False)
                    # QMessageBox.information(self,"sucess","The file is made , it will open automaticly.")
                    # open the excel automatic
                    # os.system(f'start "excel" "{file_path}"')
                if df.empty == False:
                    # self.draw(df,Pwf,wellhead_pressure,Pwh2,Pcs,Pko,total_depth,reservoir_line)
                    self.createTable(df)
                    self.grapgh(df,Pwf,wellhead_pressure,Pwh2,Pcs,Pko,total_depth,reservoir_line)
                    self.dataframe=df
                else:
                    QMessageBox.warning(self,"warning","there some error in the data")
                # except:
                #     QMessageBox.warning(self,"warning","There is some error, make sure that all variables are entred right.")
    
                        
            # except:
            #     QMessageBox.warning(self,"warning","There is some error, make sure that all variables are entred right.")
    def draw(self,df,Pwf,Pwh1,Pwh2,Pcs,Pko,total_Depth,reservoir_line):
        plt.style.use("seaborn")
        fig,ax=plt.subplots(1,1,figsize=(8,6))
        fig.suptitle("Gas_Lift",size=13,color="blue")
        # draw the formation line:
        point1,point2=[df["P2"].iloc[-1],df["Depth"].iloc[-1]], [Pwf,total_Depth]
        ax.plot([point1[0],point2[0]],[point1[1],point2[1]],label="Formation",c="#080cdd")
        # draw the tubing line
        ax.plot([Pwh1,df["P2"].iloc[-1]],[0,df["Depth"].iloc[-1]],label="Tubing",c="#139100")
        # draw the modified tubing line
        ax.plot([Pwh2,df["P2"].iloc[-1]],[0,df["Depth"].iloc[-1]],label="Safty",c="#72f25c")
        # draw the casing line 
        #get_intersection
        d,p=find_intersection(reservoir_line,[1,-Pcs/40000,Pcs])
        ax.plot([Pcs,p],[0,d],label="Casing")
        # ploting the horizontal lines and inclindes lines
        ax.plot([Pwh1,df["P1"].iloc[0]],[0,df["Depth"].iloc[0]],c="#bbb") 
        # horisontal lines
        for index,row in df.iterrows():
            ax.plot([row["P1"],row["P2"]],[row["Depth"],row["Depth"]],c="#bbb") 
        # inclinded lines
        for index,row in df.iterrows():
            if index +1 >= len(df):
                break
            else:
                ax.plot([row["P2"],df["P1"].iloc[index+1]],[row["Depth"],df["Depth"].iloc[index+1]],c="#bbb")
        # draw the line of kick_off
        if Pko != Pcs:
            ax.plot([Pko,df["P1"].iloc[0]],[0,df["Depth"].iloc[0]],label="Kick_off",c="#1e647f")
        # adding the annotaion to injection point 
        ax.annotate(f'Injection point \n {df["Depth"].iloc[-1]} ft', xy=(df["P2"].iloc[-1], df["Depth"].iloc[-1]), 
            xytext=(df["P2"].iloc[-1]-400,df["Depth"].iloc[-1]), fontsize=10,
            arrowprops=dict(facecolor='#333',width=.5,frac=.4,headwidth=4),color="#0B5345")
        plt.gca().invert_yaxis()
        ax.xaxis.tick_top()
        ax.set_xlabel("Pressure (psi)")
        ax.set_ylabel("Depth (ft)")
        ax.xaxis.set_label_position("top")
        plt.legend() 
        plt.show()
    
    def createTable(self,df):
            self.df_rows = len(df)
            self.df_cols = len(df.columns)
            self.df = df
            columns= list(self.df.columns.values)
            self.tableWidget.setRowCount(self.df_rows)
            self.tableWidget.setColumnCount(self.df_cols)
            self.tableWidget.setHorizontalHeaderLabels(columns)
            for i in range(self.df_rows):
                for j in range(self.df_cols):
                    x = f'{self.df.iloc[i, j]}'
                    self.tableWidget.setItem(i, j, QTableWidgetItem(x))
    def grapgh(self,df,Pwf,Pwh1,Pwh2,Pcs,Pko,total_Depth,reservoir_line):
        self.mat.axis.clear()
        # legend= self.mat.figure.legend()
        # to remove the old legend
        #  
        # self.mat.figure.suptitle("Gas_Lift",size=13,color="blue")
        # draw the formation line:
        point1,point2=[df["P2"].iloc[-1],df["Depth"].iloc[-1]], [Pwf,total_Depth]
        self.mat.axis.plot([point1[0],point2[0]],[point1[1],point2[1]],label="Formation",c="#080cdd")
        # draw the tubing line
        self.mat.axis.plot([Pwh1,df["P2"].iloc[-1]],[0,df["Depth"].iloc[-1]],label="Tubing",c="#139100")
        # draw the modified tubing line
        self.mat.axis.plot([Pwh2,df["P2"].iloc[-1]],[0,df["Depth"].iloc[-1]],label="Safty",c="#72f25c")
        # draw the casing line 
        #get_intersection
        d,p=find_intersection(reservoir_line,[1,-Pcs/40000,Pcs])
        self.mat.axis.plot([Pcs,p],[0,d],label="Casing")
        # ploting the horizontal lines and inclindes lines
        self.mat.axis.plot([Pwh1,df["P1"].iloc[0]],[0,df["Depth"].iloc[0]],c="#bbb") 
        # horisontal lines
        for index,row in df.iterrows():
            self.mat.axis.plot([row["P1"],row["P2"]],[row["Depth"],row["Depth"]],c="#bbb") 
        # inclinded lines
        for index,row in df.iterrows():
            if index +1 >= len(df):
                break
            else:
                self.mat.axis.plot([row["P2"],df["P1"].iloc[index+1]],[row["Depth"],df["Depth"].iloc[index+1]],c="#bbb")
        # draw the line of kick_off
        if Pko != Pcs:
            self.mat.axis.plot([Pko,df["P1"].iloc[0]],[0,df["Depth"].iloc[0]],label="Kick_off",c="#1e647f")
        # adding the annotaion to injection point 
        self.mat.axis.annotate(f'Injection point \n {df["Depth"].iloc[-1]} ft', xy=(df["P2"].iloc[-1], df["Depth"].iloc[-1]), 
            xytext=(df["P2"].iloc[-1]-400,df["Depth"].iloc[-1]), fontsize=10,
            arrowprops=dict(facecolor='#333',width=.5,frac=.4,headwidth=4),color="#0B5345")
        self.mat.figure.gca().invert_yaxis()
        self.mat.axis.xaxis.tick_top()
        self.mat.axis.xaxis.set_label_position("top")
        self.mat.axis.set_xlabel("Pressure (psi)")
        self.mat.axis.set_ylabel("Depth (ft)")
        # self.mat.figure.legend(loc="center right") 
        self.mat.canvas.draw()
    
# run on the window
def main():
    app =QApplication(sys.argv)
    window = MainApp()
    window.show()
    app.exec_()
if __name__=="__main__":
    main()

