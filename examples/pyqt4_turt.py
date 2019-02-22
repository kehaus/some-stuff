# -*- coding: utf-8 -*-
"""
Created on Tue Apr  3 23:15:39 2018

@author: kh
"""

import sys
from PyQt4 import QtGui, QtCore


#class Example(QtGui.QWidget):
#    
#    def __init__(self):
#        siper(Example, self).__init__()
#        self.initUI()
#        
#    def initUI(self):
#        
#        qbtn = QtGui.QPushButton('Quit', self)
##        qbtn.clicked.connect(Qtcore.QCoreApplication.instance().quit)
#        qbtn.resize(qbtn.sizeHint())
#        qbtn.moce(50, 50)
#        
#        self.setGeometry(300, 300, 250, 150)
#        self.setWindowTitle('Quit button')
#        self.show()
#
#
#
#
#
#
#
#def main():
#    app = QtGui.QApplication(sys.argv)
#

def main():
    
    app = QtGui.QApplication(sys.argv)
    
    w = QtGui.QWidget()
    w.resize(250, 150)
    w.move(300,300)
    w.setWindowTitle('Simple')
    
    sys.exit(app.exec_())
    
    

if __name__ == "__main__":
    main()
    pass