# -*- encoding:utf-8 -*-


from PyQt5 import QtWebEngineWidgets
import re
import sys, os
import pandas as pd
import re
from pykml.factory import KML_ElementMaker as KML
from datetime import datetime
from lxml import etree
from pykml.factory import KML_ElementMaker as KML
import folium
import folium.plugins
from folium.plugins import MarkerCluster
import numpy as np
import io
import chardet
from collections import defaultdict




def check(s):
    if s!= "":
        try:
            ss= float(s)
            return ss
        except:
            return None




def resource_path(relative_path):
    """ Get absolute path to resource, works for dev and for PyInstaller """
    base_path = getattr(sys, '_MEIPASS', os.path.dirname(os.path.abspath(__file__)))
    return os.path.join(base_path, relative_path)



class Mapper_Model:
    def __init__(self):
        d= defaultdict(lambda: None)
        self.result_files = defaultdict(lambda: None)

        self.input_files= defaultdict(lambda: None)
        self.run_files= {'user_input': d, 'file_input':d}


    # ======== WORKSPACES ============

    def filepath(self, tab, spart):
        self.input_files['tabfile']= tab
        self.input_files['spartfile']= spart


    def savepath(self, path):
        self.input_files['savefile']= path


    def user_input_data(self, plaintext):
        self.result_files['file_name']= 'input_data'
        x= plaintext
        output = io.StringIO(x)
        df1= pd.read_table(output, header= None)
        col= ['lat', 'lon', 'specimen']
        pp= int(df1.shape[-1])
        df1.columns= col[0:pp]
        if len(df1.columns)== 2: df1['specimen']= "unassigned"
        df1[['specimen']] = df1[['specimen']].fillna(value='unassigned')
        df1[['lat', 'lon']]= df1[['lat', 'lon']].applymap(lambda x: check(x))
        df1= df1.dropna()
        self.run_files['user_input']['table_file']= df1




    def user_input_html(self):
        df1= self.run_files['user_input']['table_file']
        center = [df1['lat'].values[0], df1['lon'].values[0]]
        map1 = folium.Map(location=center, zoom_start=8)
        for index, row in df1.iterrows():
            location = [row['lat'], row['lon']]
            folium.Marker(location, popup = f'Name:{row["specimen"]}').add_to(map1)
        self.result_files['html_file']= map1
        # map1.save(os.path.join(file3, "input.html"))
        # self.m_output.load(QUrl().fromLocalFile(os.path.join(file3, "input.html")))



    def user_input_kml(self):

        df1= self.run_files['user_input']['table_file']

        doc = KML.Document()
        for index, row in df1.iterrows():
            pm = KML.Placemark(
        KML.name("specimen={0}".format(row['specimen'])),


        KML.Point(
            KML.coordinates("{0},{1}".format(row['lon'],row['lat']))
        )
    )

            doc.append(pm)

        #return doc
        self.result_files['kml_file']= doc
            # outfile = open(os.path.join(file3, 'input.kml'),'w+')
            # outfile.write(etree.tostring(doc).decode('utf-8'))
            # outfile.close()


    def tab_file_data(self):
        tab_file = self.input_files['tabfile']
        spart_file = self.input_files['spartfile']
        print(self.input_files['spartfile'])
        buffer= open(tab_file, encoding="utf-8", errors="ignore")
        df1= pd.read_table(buffer)
        df1.columns= [x.lower() for x in df1.columns]
        accepted_1= ['specimen', 'specimen-voucher', 'specimenvoucher', 'specimen voucher', 'voucher', 'sample']
        accepted_2 = ['lat', 'lati', 'la']
        accepted_3= ['long', 'lon']
        df1.columns= ['specimen_voucher' if x in accepted_1 else x for x in df1.columns]
        df1.columns= ['latitude' if x in accepted_2 else x for x in df1.columns]
        df1.columns= ['longitude' if x in accepted_3 else x for x in df1.columns]
        df1['latitude']= df1['latitude'].map(lambda x: check(x))
        df1['longitude']= df1['longitude'].map(lambda x: check(x))
        df1= df1.dropna()

        if spart_file != '':
            xx= open(spart_file, mode= 'r')
            rr= xx.read()
            bb= re.findall(r'(?<=N_spartitions)[^A-Za-z]*(\w+\W+.*)(;)', rr)
            bb= bb[0][0]
            bb= bb.split(';')[0]
            aa= re.findall(r'(?sm)(?<=Individual_assignment)[^A-Za-z]*(\w+\W+.*)(;)', rr)
            aa= aa[0][0]
            aa= aa.strip().split(';', 1)
            aa= aa[0]
            aa= aa.split('\n')
            d= {'specimen_voucher': aa}
            dd= pd.DataFrame(d)
            dd['species_number']= dd['specimen_voucher'].map(lambda x: x.split(':')[-1])
            dd['specimen_voucher']= dd['specimen_voucher'].map(lambda x: x.split(':')[0])
            dd= dd.applymap(lambda x: x.strip(';'))
            dd= dd.applymap(lambda x: x.strip())
            df1= df1.merge(dd, on= 'specimen_voucher', how= 'left')
            df1= df1.fillna("0")

        elif spart_file == '':
            df1['species_number']= "0"
            bb= os.path.splitext(os.path.basename(tab_file))[0]

        buffer.close()
        self.run_files['file_input']['table_file']= df1
        self.run_files['file_input']['file_name']= bb



    def tab_kml(self):

        df1= self.run_files['file_input']['table_file']
        bb = self.run_files['file_input']['file_name']

        stylename = "earthquake-balloon-style"
        balloonstyle = KML.BalloonStyle(
          KML.text("""
        <table Border=1>
        <tr><th>Earthquake ID</th><td>$[Eqid_]</td></tr>
        <tr><th>Magnitude</th><td>$[Magnitude_]</td></tr>
        <tr><th>Depth</th><td>$[Depth_]</td></tr>
        <tr><th>Datetime</th><td>$[Datetime_]</td></tr>
        <tr><th>Coordinates</th><td>($[Lat_],$[Lat_])</td></tr>
        <tr><th>Region</th><td>$[Region_]</td></tr>
        </table>"""
          ),
        )


        doc = KML.Document()

        iconstyles = [
          [1,'ff000000'],
          [2,'ffff0000'],
          [3,'ff00ff55'],
          [4,'ffff00aa'],
          [5,'ff00ffff'],
          [6,'ff0000ff'],
        ]

        for threshold,color in iconstyles:
          doc.append(
              KML.Style(
                  KML.IconStyle(
                      KML.color(color),
                    KML.scale(1.1),
                    KML.Icon(
                        KML.href('http://www.gstatic.com/mapspro/images/stock/503-wht-blank_maps.png')
                    ),
                    KML.hotspot('', x='16', y='31', xunits='pixels', yunits='insetPixels'),
                ),

                KML.LabelStyle(
                    KML.scale(1.1)
                ),

                id='icon-503-{}-normal'.format(color)))


        def element(g, i):
            for k, row in g.iterrows():
                if len(iconstyles) > int(i):
                    jj= int(i)+1
                elif len(iconstyles) <= int(i):
                    jj= int(i)-(len(iconstyles)-1)
                pm = KML.Placemark(
                    KML.name("specimen={0}".format(row['specimen_voucher'])),
                     KML.styleUrl('#icon-503-{}-normal'.format(iconstyles[jj][-1])),

                    KML.Point(
                        KML.coordinates("{0},{1}".format(row['longitude'],row['latitude']))
                    )
                )

                fld.append(pm)
            return fld

        for i, g in df1.groupby('species_number'):

            fld= KML.Folder(KML.name(f'species number: {"unassigned" if i=="0" else i}'))

            gg= element(g, i)
            doc.append(gg)

        self.result_files['kml_file']= doc
        self.result_files['file_name']= bb


            # outfile = open(os.path.join(file3, bb+'.kml'),'w+')
            # outfile.write(etree.tostring(doc).decode('utf-8'))
            # outfile.close()

        # except Exception as e:
        #     QMessageBox.warning(self, "Warning", f"The spartmapper output not obtained, please check input file type because {e}")
        #     return
        # QMessageBox.information(self, "Information", "The spartmapper output data generated successfully")


    def tab_html(self):
        df1= self.run_files['file_input']['table_file']
        bb = self.run_files['file_input']['file_name']
        m = folium.Map(location=[df1['latitude'].values[0], df1['longitude'].values[0]], zoom_start=6)
        fg = folium.FeatureGroup(name= bb)
        m.add_child(fg)
        color= ['orange', 'purple', 'red', 'blue', 'green', 'pink', 'yellow']

        for i, g in df1.groupby('species_number'):
            if int(i) < len(color):
                cc= int(i)
            elif int(i)>= len(color):
                cc= int(i)- len(color)

            g1= folium.plugins.FeatureGroupSubGroup(fg, f'species number: {"unassigned" if i=="0" else i}')
            m.add_child(g1)
            for k, row in g.iterrows():

                folium.Marker(location=[float(row['latitude']), float(row['longitude'])], popup= row["specimen_voucher"], icon=folium.Icon(color=color[cc])).add_to(g1)

        folium.LayerControl(collapsed=False).add_to(m)
        #return m
        self.result_files['html_file']= m
        self.result_files['file_name']= bb

        # m.save(os.path.join(file3, bb+".html"))
        #
        # self.m_output.load(QUrl().fromLocalFile(os.path.join(file3, bb+".html")))
        # self.m_output.setZoomFactor(1.5)
        # self.m_output.resize(640, 400)
        # self.m_output.show()
