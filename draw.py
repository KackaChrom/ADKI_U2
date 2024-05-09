from PyQt6.QtCore import *
from PyQt6.QtGui import *
from PyQt6.QtGui import QMouseEvent, QPaintEvent
from PyQt6.QtWidgets import *
from algorithms import *


class Draw(QWidget):
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        
        self.building = []
        self.mbr = QPolygonF()
        self.er = QPolygonF()
        self.wa = QPolygonF()
        self.le = QPolygonF()
        self.wb = QPolygonF()
        self.ch = QPolygonF()


    def mousePressEvent(self, e: QMouseEvent):
        
        # Get coordinates of q
        x = e.position().x()
        y = e.position().y()
        
        # Add new vertex
        p = QPointF(x, y)
        
        # Add to polygon
        if not self.building:  # if not new
            self.building.append(QPolygonF())
        self.building[-1].append(p)
            
        #Repaint screen
        self.repaint()


    def paintEvent(self, e: QPaintEvent):
        # Draw situation

        # Create new graphic object
        qp = QPainter(self)

        # Start drawing
        qp.begin(self)

        # Set graphical attributes for buildings
        qp.setPen(Qt.GlobalColor.black)
        qp.setBrush(Qt.GlobalColor.yellow)

        # Draw building
        # Ensure 'self.building' is a list even if it's a single building
        if isinstance(self.building, QPolygon):
            self.building = [self.building]  # Convert to list if it's a single QPolygon

        # Draw polygon for each building
        for building in self.building:
            qp.drawPolygon(building)
        
        # MBR ------
        # Set graphical attributes for MBR
        qp.setPen(Qt.GlobalColor.red)
        qp.setBrush(Qt.GlobalColor.transparent)
        
        # Draw MBRs
        for mbr in self.mbr:
            qp.drawPolygon(mbr)
        
        # PCA ------
        # Set graphical attributes for PCA
        qp.setPen(Qt.GlobalColor.blue)
        qp.setBrush(Qt.GlobalColor.transparent)
        
        # Draw PCAs
        for er in self.er:
            qp.drawPolygon(er)
        
        # WA ------
        # Set graphical attributes for WA
        qp.setPen(Qt.GlobalColor.green)
        qp.setBrush(Qt.GlobalColor.transparent)
        
        # Draw WAs
        for wa in self.wa:
            qp.drawPolygon(wa)
            
        # LE ------
        # Set graphical attributes for LE
        qp.setPen(Qt.GlobalColor.gray)
        qp.setBrush(Qt.GlobalColor.transparent)
        
        # Draw LEs
        for le in self.le:
            qp.drawPolygon(le)
            
        # WB ------
        # Set graphical attributes for WB
        qp.setPen(Qt.GlobalColor.darkGreen)
        qp.setBrush(Qt.GlobalColor.transparent)
        
        # Draw WBs
        for wb in self.wb:
            qp.drawPolygon(wb)
        
        # End drawing
        qp.end()

    def getBuilding(self):
        # Return building
        return self.building

    def setMBR(self, mbrs):
        # set mbr to input
        self.mbr = mbrs  

    def setData(self, build):
        self.clearData()
        
        self.building = build
        
        self.repaint()

    def setER(self, ers):
        self.er = ers
    
    def setWA(self, was):
        self.wa = was
        
    def setLE(self, les):
        self.le = les
    
    def setWB(self, wbs):
        self.wb = wbs

    def setCH(self,ch):
        #Set ch tp input
        self.ch = ch  

    def clearData(self):
        # Clear results
        self.mbr = QPolygonF()
        self.er = QPolygonF()
        self.ch = QPolygon()
        self.wa = QPolygonF()
        self.le = QPolygonF()
        self.wb = QPolygonF()
        
    def clearAllData(self):
        # Clear everything
        self.building.clear()
        self.mbr = QPolygonF()
        self.er = QPolygonF()
        self.ch = QPolygon()
        self.wa = QPolygonF()
        self.le = QPolygonF()
        self.wb = QPolygonF()