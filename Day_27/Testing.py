# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'Testing.ui'
#
# Created by: PyQt5 UI code generator 5.15.6
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(546, 515)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.tabWidget = QtWidgets.QTabWidget(self.centralwidget)
        self.tabWidget.setGeometry(QtCore.QRect(10, 10, 531, 461))
        self.tabWidget.setAutoFillBackground(False)
        self.tabWidget.setTabPosition(QtWidgets.QTabWidget.North)
        self.tabWidget.setElideMode(QtCore.Qt.ElideNone)
        self.tabWidget.setDocumentMode(False)
        self.tabWidget.setTabsClosable(False)
        self.tabWidget.setMovable(True)
        self.tabWidget.setTabBarAutoHide(False)
        self.tabWidget.setObjectName("tabWidget")

        self.Binomial = QtWidgets.QWidget()
        self.Binomial.setObjectName("Binomial")
        self.label = QtWidgets.QLabel(self.Binomial)
        self.label.setGeometry(QtCore.QRect(10, 20, 451, 31))
        font = QtGui.QFont()
        font.setPointSize(18)
        self.label.setFont(font)
        self.label.setObjectName("label")
        self.label_2 = QtWidgets.QLabel(self.Binomial)
        self.label_2.setGeometry(QtCore.QRect(10, 50, 311, 41))
        font = QtGui.QFont()
        font.setPointSize(18)
        self.label_2.setFont(font)
        self.label_2.setObjectName("label_2")
        self.label_3 = QtWidgets.QLabel(self.Binomial)
        self.label_3.setGeometry(QtCore.QRect(10, 90, 361, 31))
        font = QtGui.QFont()
        font.setPointSize(18)
        self.label_3.setFont(font)
        self.label_3.setObjectName("label_3")
        self.frame = QtWidgets.QFrame(self.Binomial)
        self.frame.setGeometry(QtCore.QRect(10, 160, 501, 261))
        self.frame.setFrameShape(QtWidgets.QFrame.Box)
        self.frame.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame.setLineWidth(1)
        self.frame.setMidLineWidth(0)
        self.frame.setObjectName("frame")
        self.label_4 = QtWidgets.QLabel(self.frame)
        self.label_4.setGeometry(QtCore.QRect(20, 20, 47, 31))
        font = QtGui.QFont()
        font.setPointSize(14)
        self.label_4.setFont(font)
        self.label_4.setObjectName("label_4")
        self.label_5 = QtWidgets.QLabel(self.frame)
        self.label_5.setGeometry(QtCore.QRect(20, 60, 47, 21))
        font = QtGui.QFont()
        font.setPointSize(14)
        self.label_5.setFont(font)
        self.label_5.setObjectName("label_5")
        self.label_6 = QtWidgets.QLabel(self.frame)
        self.label_6.setGeometry(QtCore.QRect(20, 90, 47, 31))
        font = QtGui.QFont()
        font.setPointSize(14)
        self.label_6.setFont(font)
        self.label_6.setObjectName("label_6")
        self.b_x = QtWidgets.QLineEdit(self.frame)
        self.b_x.setGeometry(QtCore.QRect(50, 20, 91, 21))
        self.b_x.setText("")
        self.b_x.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.b_x.setObjectName("b_x")
        self.b_x.setText("TESGING")
        self.b_x.get

        self.b_n = QtWidgets.QLineEdit(self.frame)
        self.b_n.setGeometry(QtCore.QRect(50, 60, 91, 21))
        self.b_n.setText("")
        self.b_n.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.b_n.setObjectName("b_n")
        self.b_p = QtWidgets.QLineEdit(self.frame)
        self.b_p.setGeometry(QtCore.QRect(50, 100, 91, 21))
        self.b_p.setText("")
        self.b_p.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.b_p.setObjectName("b_p")
        self.frame_2 = QtWidgets.QFrame(self.frame)
        self.frame_2.setGeometry(QtCore.QRect(170, 20, 311, 211))
        self.frame_2.setFrameShape(QtWidgets.QFrame.Box)
        self.frame_2.setFrameShadow(QtWidgets.QFrame.Plain)
        self.frame_2.setLineWidth(2)
        self.frame_2.setMidLineWidth(2)
        self.frame_2.setObjectName("frame_2")
        self.b_results = QtWidgets.QTextEdit(self.frame_2)
        self.b_results.setGeometry(QtCore.QRect(0, 0, 311, 211))
        self.b_results.setReadOnly(True)
        self.b_results.setObjectName("b_results")
        self.b_calculate = QtWidgets.QPushButton(self.frame)
        self.b_calculate.setGeometry(QtCore.QRect(20, 170, 131, 51))
        font = QtGui.QFont()
        font.setPointSize(15)
        font.setBold(False)
        font.setWeight(50)
        self.b_calculate.setFont(font)
        self.b_calculate.setCheckable(False)
        self.b_calculate.setFlat(False)
        self.b_calculate.setObjectName("b_calculate")
        self.tabWidget.addTab(self.Binomial, "")

        self.Geometric = QtWidgets.QWidget()
        self.Geometric.setObjectName("Geometric")
        self.label_7 = QtWidgets.QLabel(self.Geometric)
        self.label_7.setGeometry(QtCore.QRect(10, 20, 451, 31))
        font = QtGui.QFont()
        font.setPointSize(18)
        self.label_7.setFont(font)
        self.label_7.setObjectName("label_7")
        self.frame_3 = QtWidgets.QFrame(self.Geometric)
        self.frame_3.setGeometry(QtCore.QRect(10, 160, 501, 261))
        self.frame_3.setFrameShape(QtWidgets.QFrame.Box)
        self.frame_3.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_3.setLineWidth(1)
        self.frame_3.setMidLineWidth(0)
        self.frame_3.setObjectName("frame_3")
        self.label_8 = QtWidgets.QLabel(self.frame_3)
        self.label_8.setGeometry(QtCore.QRect(20, 20, 47, 31))
        font = QtGui.QFont()
        font.setPointSize(14)
        self.label_8.setFont(font)
        self.label_8.setObjectName("label_8")
        self.label_10 = QtWidgets.QLabel(self.frame_3)
        self.label_10.setGeometry(QtCore.QRect(20, 50, 47, 31))
        font = QtGui.QFont()
        font.setPointSize(14)
        self.label_10.setFont(font)
        self.label_10.setObjectName("label_10")
        self.g_x = QtWidgets.QLineEdit(self.frame_3)
        self.g_x.setGeometry(QtCore.QRect(50, 20, 91, 21))
        self.g_x.setText("")
        self.g_x.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.g_x.setObjectName("g_x")
        self.g_p = QtWidgets.QLineEdit(self.frame_3)
        self.g_p.setGeometry(QtCore.QRect(50, 60, 91, 21))
        self.g_p.setText("")
        self.g_p.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.g_p.setObjectName("g_p")
        self.frame_4 = QtWidgets.QFrame(self.frame_3)
        self.frame_4.setGeometry(QtCore.QRect(170, 20, 311, 211))
        self.frame_4.setFrameShape(QtWidgets.QFrame.Box)
        self.frame_4.setFrameShadow(QtWidgets.QFrame.Plain)
        self.frame_4.setLineWidth(2)
        self.frame_4.setMidLineWidth(2)
        self.frame_4.setObjectName("frame_4")
        self.g_results = QtWidgets.QTextEdit(self.frame_4)
        self.g_results.setGeometry(QtCore.QRect(0, 0, 311, 211))
        self.g_results.setReadOnly(True)
        self.g_results.setObjectName("g_results")
        self.g_calculate = QtWidgets.QPushButton(self.frame_3)
        self.g_calculate.setGeometry(QtCore.QRect(20, 170, 131, 51))
        font = QtGui.QFont()
        font.setPointSize(15)
        font.setBold(False)
        font.setWeight(50)
        self.g_calculate.setFont(font)
        self.g_calculate.setCheckable(False)
        self.g_calculate.setFlat(False)
        self.g_calculate.setObjectName("g_calculate")
        self.label_11 = QtWidgets.QLabel(self.Geometric)
        self.label_11.setGeometry(QtCore.QRect(10, 60, 361, 31))
        font = QtGui.QFont()
        font.setPointSize(18)
        self.label_11.setFont(font)
        self.label_11.setObjectName("label_11")
        self.tabWidget.addTab(self.Geometric, "")

        self.NegativBinomial = QtWidgets.QWidget()
        self.NegativBinomial.setObjectName("NegativBinomial")
        self.label_13 = QtWidgets.QLabel(self.NegativBinomial)
        self.label_13.setGeometry(QtCore.QRect(10, 20, 451, 31))
        font = QtGui.QFont()
        font.setPointSize(18)
        self.label_13.setFont(font)
        self.label_13.setObjectName("label_13")
        self.frame_5 = QtWidgets.QFrame(self.NegativBinomial)
        self.frame_5.setGeometry(QtCore.QRect(10, 160, 501, 261))
        self.frame_5.setFrameShape(QtWidgets.QFrame.Box)
        self.frame_5.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_5.setLineWidth(1)
        self.frame_5.setMidLineWidth(0)
        self.frame_5.setObjectName("frame_5")
        self.label_14 = QtWidgets.QLabel(self.frame_5)
        self.label_14.setGeometry(QtCore.QRect(20, 20, 47, 31))
        font = QtGui.QFont()
        font.setPointSize(14)
        self.label_14.setFont(font)
        self.label_14.setObjectName("label_14")
        self.label_15 = QtWidgets.QLabel(self.frame_5)
        self.label_15.setGeometry(QtCore.QRect(20, 60, 47, 21))
        font = QtGui.QFont()
        font.setPointSize(14)
        self.label_15.setFont(font)
        self.label_15.setObjectName("label_15")
        self.label_16 = QtWidgets.QLabel(self.frame_5)
        self.label_16.setGeometry(QtCore.QRect(20, 90, 47, 31))
        font = QtGui.QFont()
        font.setPointSize(14)
        self.label_16.setFont(font)
        self.label_16.setObjectName("label_16")
        self.n_x = QtWidgets.QLineEdit(self.frame_5)
        self.n_x.setGeometry(QtCore.QRect(50, 20, 91, 21))
        self.n_x.setText("")
        self.n_x.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.n_x.setObjectName("n_x")
        self.n_r = QtWidgets.QLineEdit(self.frame_5)
        self.n_r.setGeometry(QtCore.QRect(50, 60, 91, 21))
        self.n_r.setText("")
        self.n_r.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.n_r.setObjectName("n_r")
        self.n_p = QtWidgets.QLineEdit(self.frame_5)
        self.n_p.setGeometry(QtCore.QRect(50, 100, 91, 21))
        self.n_p.setText("")
        self.n_p.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.n_p.setObjectName("n_p")
        self.frame_6 = QtWidgets.QFrame(self.frame_5)
        self.frame_6.setGeometry(QtCore.QRect(170, 20, 311, 211))
        self.frame_6.setFrameShape(QtWidgets.QFrame.Box)
        self.frame_6.setFrameShadow(QtWidgets.QFrame.Plain)
        self.frame_6.setLineWidth(2)
        self.frame_6.setMidLineWidth(2)
        self.frame_6.setObjectName("frame_6")
        self.n_results = QtWidgets.QTextEdit(self.frame_6)
        self.n_results.setGeometry(QtCore.QRect(0, 0, 311, 211))
        self.n_results.setReadOnly(True)
        self.n_results.setObjectName("n_results")
        self.n_calculate = QtWidgets.QPushButton(self.frame_5)
        self.n_calculate.setGeometry(QtCore.QRect(20, 170, 131, 51))
        font = QtGui.QFont()
        font.setPointSize(15)
        font.setBold(False)
        font.setWeight(50)
        self.n_calculate.setFont(font)
        self.n_calculate.setCheckable(False)
        self.n_calculate.setFlat(False)
        self.n_calculate.setObjectName("n_calculate")
        self.label_17 = QtWidgets.QLabel(self.NegativBinomial)
        self.label_17.setGeometry(QtCore.QRect(10, 90, 361, 31))
        font = QtGui.QFont()
        font.setPointSize(18)
        self.label_17.setFont(font)
        self.label_17.setObjectName("label_17")
        self.label_18 = QtWidgets.QLabel(self.NegativBinomial)
        self.label_18.setGeometry(QtCore.QRect(10, 50, 301, 41))
        font = QtGui.QFont()
        font.setPointSize(18)
        self.label_18.setFont(font)
        self.label_18.setObjectName("label_18")

        self.tabWidget.addTab(self.NegativBinomial, "")
        self.Hypergeometric = QtWidgets.QWidget()
        self.Hypergeometric.setObjectName("Hypergeometric")
        self.label_19 = QtWidgets.QLabel(self.Hypergeometric)
        self.label_19.setGeometry(QtCore.QRect(10, 20, 451, 31))
        font = QtGui.QFont()
        font.setPointSize(18)
        self.label_19.setFont(font)
        self.label_19.setObjectName("label_19")
        self.frame_7 = QtWidgets.QFrame(self.Hypergeometric)
        self.frame_7.setGeometry(QtCore.QRect(10, 160, 501, 261))
        self.frame_7.setFrameShape(QtWidgets.QFrame.Box)
        self.frame_7.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame_7.setLineWidth(1)
        self.frame_7.setMidLineWidth(0)
        self.frame_7.setObjectName("frame_7")
        self.label_20 = QtWidgets.QLabel(self.frame_7)
        self.label_20.setGeometry(QtCore.QRect(20, 20, 47, 31))
        font = QtGui.QFont()
        font.setPointSize(14)
        self.label_20.setFont(font)
        self.label_20.setObjectName("label_20")
        self.label_21 = QtWidgets.QLabel(self.frame_7)
        self.label_21.setGeometry(QtCore.QRect(20, 60, 47, 21))
        font = QtGui.QFont()
        font.setPointSize(14)
        self.label_21.setFont(font)
        self.label_21.setObjectName("label_21")
        self.label_22 = QtWidgets.QLabel(self.frame_7)
        self.label_22.setGeometry(QtCore.QRect(20, 90, 47, 31))
        font = QtGui.QFont()
        font.setPointSize(14)
        self.label_22.setFont(font)
        self.label_22.setObjectName("label_22")
        self.h_x = QtWidgets.QLineEdit(self.frame_7)
        self.h_x.setGeometry(QtCore.QRect(50, 20, 91, 21))
        self.h_x.setText("")
        self.h_x.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.h_x.setObjectName("h_x")
        self.h_N = QtWidgets.QLineEdit(self.frame_7)
        self.h_N.setGeometry(QtCore.QRect(50, 60, 91, 21))
        self.h_N.setText("")
        self.h_N.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.h_N.setObjectName("h_N")
        self.h_n = QtWidgets.QLineEdit(self.frame_7)
        self.h_n.setGeometry(QtCore.QRect(50, 100, 91, 21))
        self.h_n.setText("")
        self.h_n.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.h_n.setObjectName("h_n")
        self.frame_8 = QtWidgets.QFrame(self.frame_7)
        self.frame_8.setGeometry(QtCore.QRect(170, 20, 311, 211))
        self.frame_8.setFrameShape(QtWidgets.QFrame.Box)
        self.frame_8.setFrameShadow(QtWidgets.QFrame.Plain)
        self.frame_8.setLineWidth(2)
        self.frame_8.setMidLineWidth(2)
        self.frame_8.setObjectName("frame_8")
        self.h_results = QtWidgets.QTextEdit(self.frame_8)
        self.h_results.setGeometry(QtCore.QRect(0, 0, 311, 211))
        self.h_results.setReadOnly(True)
        self.h_results.setObjectName("h_results")
        self.h_calculate = QtWidgets.QPushButton(self.frame_7)
        self.h_calculate.setGeometry(QtCore.QRect(20, 170, 131, 51))
        font = QtGui.QFont()
        font.setPointSize(15)
        font.setBold(False)
        font.setWeight(50)
        self.h_calculate.setFont(font)
        self.h_calculate.setCheckable(False)
        self.h_calculate.setFlat(False)
        self.h_calculate.setObjectName("h_calculate")
        self.h_r = QtWidgets.QLineEdit(self.frame_7)
        self.h_r.setGeometry(QtCore.QRect(50, 130, 91, 21))
        self.h_r.setText("")
        self.h_r.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.h_r.setObjectName("h_r")
        self.label_26 = QtWidgets.QLabel(self.frame_7)
        self.label_26.setGeometry(QtCore.QRect(20, 130, 47, 21))
        font = QtGui.QFont()
        font.setPointSize(14)
        self.label_26.setFont(font)
        self.label_26.setObjectName("label_26")
        self.label_23 = QtWidgets.QLabel(self.Hypergeometric)
        self.label_23.setGeometry(QtCore.QRect(10, 90, 351, 21))
        font = QtGui.QFont()
        font.setPointSize(18)
        self.label_23.setFont(font)
        self.label_23.setObjectName("label_23")
        self.label_24 = QtWidgets.QLabel(self.Hypergeometric)
        self.label_24.setGeometry(QtCore.QRect(10, 50, 311, 41))
        font = QtGui.QFont()
        font.setPointSize(18)
        self.label_24.setFont(font)
        self.label_24.setObjectName("label_24")
        self.label_25 = QtWidgets.QLabel(self.Hypergeometric)
        self.label_25.setGeometry(QtCore.QRect(10, 110, 361, 31))
        font = QtGui.QFont()
        font.setPointSize(18)
        self.label_25.setFont(font)
        self.label_25.setObjectName("label_25")
        self.tabWidget.addTab(self.Hypergeometric, "")
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 546, 22))
        self.menubar.setObjectName("menubar")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)

        self.retranslateUi(MainWindow)
        self.tabWidget.setCurrentIndex(3)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow"))
        self.label.setText(_translate("MainWindow", "x: Number of successes seeking"))
        self.label_2.setText(_translate("MainWindow", "n: Total population"))
        self.label_3.setText(_translate("MainWindow", "p: Probability of success"))
        self.label_4.setText(_translate("MainWindow", "x:"))
        self.label_5.setText(_translate("MainWindow", "n:"))
        self.label_6.setText(_translate("MainWindow", "p:"))
        self.b_results.setHtml(_translate("MainWindow", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'MS Shell Dlg 2\'; font-size:7.875pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:12px; margin-bottom:12px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p></body></html>"))
        self.b_calculate.setText(_translate("MainWindow", "Calculate"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.Binomial), _translate("MainWindow", "Binomial"))
        self.label_7.setText(_translate("MainWindow", "x: Index of first success"))
        self.label_8.setText(_translate("MainWindow", "x:"))
        self.label_10.setText(_translate("MainWindow", "p:"))
        self.g_results.setHtml(_translate("MainWindow", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'MS Shell Dlg 2\'; font-size:7.875pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:12px; margin-bottom:12px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p></body></html>"))
        self.g_calculate.setText(_translate("MainWindow", "Calculate"))
        self.label_11.setText(_translate("MainWindow", "p: Probability of success"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.Geometric), _translate("MainWindow", "Geometric"))
        self.label_13.setText(_translate("MainWindow", "x: Number of trials"))
        self.label_14.setText(_translate("MainWindow", "x:"))
        self.label_15.setText(_translate("MainWindow", "r:"))
        self.label_16.setText(_translate("MainWindow", "p:"))
        self.n_results.setHtml(_translate("MainWindow", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'MS Shell Dlg 2\'; font-size:7.875pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:12px; margin-bottom:12px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p></body></html>"))
        self.n_calculate.setText(_translate("MainWindow", "Calculate"))
        self.label_17.setText(_translate("MainWindow", "p: Probability of success"))
        self.label_18.setText(_translate("MainWindow", "r: Number of successes on x th trial"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.NegativBinomial), _translate("MainWindow", "Negative Binomial"))
        self.label_19.setText(_translate("MainWindow", "x: Number of successes seeking"))
        self.label_20.setText(_translate("MainWindow", "x:"))
        self.label_21.setText(_translate("MainWindow", "N:"))
        self.label_22.setText(_translate("MainWindow", "n:"))
        self.h_results.setHtml(_translate("MainWindow", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'MS Shell Dlg 2\'; font-size:7.875pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:12px; margin-bottom:12px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><br /></p></body></html>"))
        self.h_calculate.setText(_translate("MainWindow", "Calculate"))
        self.label_26.setText(_translate("MainWindow", "r:"))
        self.label_23.setText(_translate("MainWindow", "n: Number of selections made"))
        self.label_24.setText(_translate("MainWindow", "N: Total population size"))
        self.label_25.setText(_translate("MainWindow", "r: Total number of successes"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.Hypergeometric), _translate("MainWindow", "Hypergeometric"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec_())