# -*- encoding:utf-8 -*-
from PyQt5 import QtWidgets, QtCore, QtGui
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
import re
import sys, os
from PyQt5.uic import loadUiType
import pandas as pd
from scipy import stats
from PyQt5.QtGui import *
from collections import defaultdict
import tempfile
from MolD_sDNCFASTA import *
import time, datetime
import asyncio
from PyQt5.QtCore import QThread


def resource_path(relative_path):
    """ Get absolute path to resource, works for dev and for PyInstaller """
    base_path = getattr(sys, '_MEIPASS', os.path.dirname(os.path.abspath(__file__)))
    return os.path.join(base_path, relative_path)



FORM_CLASS,_=loadUiType(resource_path("mold.ui"))

class MyAbstract(QThread):
    """Base export thread"""

    loop = asyncio.get_event_loop()
    def __init__(self, func, parent=None):
        super().__init__(parent)
        self.func= func


    def run(self):
        loop.run_until_complete(self.func())


class Main(QDialog, FORM_CLASS):
    def __init__(self,parent=None):
        super(Main, self).__init__(parent)
        self.setupUi(self)
        self.setWindowIcon(QIcon(resource_path('icon/stat.ico')))
        self.setWindowTitle("Mold")
        self.progressBar.setStyle(QStyleFactory.create("windows"))
        self.progressBar.setRange(0, 1)
        self.launcher= MyAbstract(self.run)
        self.gp.setVisible(False)
        self.reset_placement()
        self.filepath= defaultdict(lambda: None)
        self.outpath= defaultdict(lambda: None)
        self.toolButton_2.setEnabled(False)
        self.toolButton_3.setEnabled(False)
        self.plainTextEdit.setToolTip("newline separated list of taxa or ALL (if each taxon is to be diagnosed) if all taxa with more than N sequences available (where N is a natural number), enter >N")
        self.cb2.setToolTip("default 'no', 'yes' (if alignment gaps to be coded as a character), or 'no' if not.")
        self.cb1.setToolTip("1 for species, 2 for above species")
        self.m1.setToolTip("default 100 Number of informative positions for focus taxon to be considered, natural number")
        self.m2.setToolTip("default 25 Allowed number of ambiguously called nucleotides per sequence, natural number")
        self.m3.setToolTip("default 10000")
        self.m4.setToolTip("default 12")
        self.m5.setToolTip("default 7")
        self.r1.setToolTip("Percent difference between original and modified sequence (default 1 for species-level taxa, 3 for for supraspecific taxa)")
        self.r2.setToolTip("Max number of sequences per taxon to modify (default 10)")
        self.r3.setToolTip("To score each candidate rDNC, 100 simulated test datasets are created. If rDNC remains valid in a test dataset, it adds 1 to the score, so lowest possible score is 0 and highest is 100. If two consecutive scores are above the threshold value (defined here) the rDNC is output."+
                            "The threshold values are like: \n lousy: 66 moderate: 75 stringent: 90 very_stringent: 95 default: moderate")
        #self.toolButton_2.setEnabled(False)
        #self.toolButton_3.setEnabled(False)
        self.temporary= tempfile.TemporaryDirectory()
        quit = QAction("Quit", self)
        quit.triggered.connect(self.closeEvent)
        self.plainTextEdit.setPlainText('ALL')

        self.m1.setText('100')
        self.m2.setText('25')
        self.m3.setText('10000')
        self.m4.setText('12')
        self.m5.setText('7')
        self.r1.setText('1')
        self.r2.setText('10')
        self.r3.setText('75')
        self.m1.setValidator(QIntValidator())
        self.m2.setValidator(QIntValidator())
        self.m3.setValidator(QIntValidator())
        self.m4.setValidator(QIntValidator())
        self.m5.setValidator(QIntValidator())
        self.r1.setValidator(QDoubleValidator(0.99,99.99,2))
        self.r2.setValidator(QIntValidator())
        self.r3.setValidator(QIntValidator())

        self.Handel_Buttons()


    def closeEvent(self, event):
         close = QMessageBox.question(self, "QUIT", "Are you sure you want to close the program?",QMessageBox.Yes | QMessageBox.No)
         if close == QMessageBox.Yes:
             self.plainTextEdit.clear()
             self.temporary.cleanup()
             #os.removedirs(temp_directory)
             event.accept()
         else:
             event.ignore()


    #======= SETUP slots =================================



    def Handel_Buttons(self):
        self.toolButton.clicked.connect(self.open_file)
        self.toolButton_2.clicked.connect(self.download2)
        self.toolButton_3.clicked.connect(self.save_file)
        self.toolButton_4.clicked.connect(self.clear)
        self.pushButton.clicked.connect(self.BtnHandler)


    def download2(self):

        def started():

            self.progressBar.setRange(0, 0)


        def finished():
            self.progressBar.setRange(0, 1)



        self.launcher.started.connect(started)
        self.launcher.finished.connect(finished)

        self.launcher.start()



    def open_file(self):
        msg = 'The input file has each line corresponds to one sequenced specimen and contains three space‐separated records:\nspecies name\nfocus taxon name of the query level\n nucleotide sequence'
        QMessageBox.information(self, 'Add input file', msg)
        sel = 'Select text file'
        tab = self.file_dialog(sel, ".")
        if tab:
            print(tab)

            self.filepath['path']= tab
        self.toolButton_2.setEnabled(True)




    def file_dialog(self, msg, path):
        return QFileDialog.getOpenFileName(self, msg, path)[0]



    async def run(self):
        try:
            self.toolButton_2.setEnabled(True)
            self.toolButton_3.setEnabled(True)
            print(self.plainTextEdit.toPlainText())
            xx= self.plainTextEdit.toPlainText().split('\n')
            print(xx)
            gapsaschars=self.cb2.currentText()
            taxalist=self.plainTextEdit.toPlainText()
            taxonrank=int(self.cb1.currentText())
            cutoff= int(self.m1.text())
            numnucl= int(self.m2.text())
            numiter= int(self.m3.text())
            maxlenraw= int(self.m4.text())
            maxlenrefined= int(self.m5.text())
            pdiff= float(self.r1.text())
            nmax= int(self.r2.text())
            thresh= int(self.r3.text())
            tmpfname=self.filepath['path']
            origfname=None
            f= self.temporary
            _, filename= os.path.split(tmpfname)
            filename= filename.split(".")[0]
            output= os.path.join(f.name, filename + ".html")
            self.outpath['output']= output
            print(self.outpath['output'])
            mainprocessing(gapsaschars=gapsaschars, taxalist=taxalist, taxonrank=taxonrank, cutoff=cutoff, numnucl=numnucl, numiter=numiter, maxlenraw=maxlenraw, maxlenrefined=maxlenrefined, pdiff=pdiff, nmax=nmax, thresh=thresh, tmpfname=tmpfname, origfname=None, output=output)

        except Exception as e:
            print(e)
            QMessageBox.warning(self, "Warning", f"The output  is not obtained because {e}")




    def save_file(self):
        try:

            from pathlib import Path
            msg = 'Please browse to output file having .html extension'
            QMessageBox.information(self, 'Browse output file', msg)
            sel = 'Select html file'
            path1, ok = QFileDialog.getSaveFileName(self, caption="save_html", filter="HTML Files (*.HTML *.html)")
            if ok:
                print("PATH", path1)

            Path(self.outpath['output']).replace(path1)

        except Exception as e:
            print(e)
            QMessageBox.warning(self, "Warning", f"The output  is not saved because {e}")





    def clear(self):

        self.toolButton_2.setEnabled(False)
        self.toolButton_3.setEnabled(False)
        self.plainTextEdit.setPlainText('ALL')
        self.m1.setText('100')
        self.m2.setText('25')
        self.m3.setText('10000')
        self.m4.setText('12')
        self.m5.setText('7')
        self.r1.setText('1')
        self.r2.setText('10')
        self.r3.setText('75')
        self.cb1.setCurrentIndex(0)
        self.cb2.setCurrentIndex(0)







    def BtnHandler(self):
        if self.pushButton.isChecked():
            self.gp.setVisible(True)
            print(self.sizeHint())
        else:

            self.gp.setVisible(False)
            self.reset_placement()


    def reset_placement(self):

        g = QDesktopWidget().availableGeometry()
        self.resize(0.1 * g.width(), 0.4 * g.height())
        self.move(g.center().x() - self.width() / 2, g.center().y() - self.height() / 2)





def main1():

    app=QApplication(sys.argv)
    window=Main()
    window.show()
    app.exec_()


if __name__=='__main__':
    main1()
