from PyQt6.QtCore import *
from PyQt6.QtGui import *
import fiona
from shapely.geometry import shape, Polygon
from PyQt6.QtWidgets import QFileDialog, QWidget

class pio(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.dia = QFileDialog(self)
        self.dia.setNameFilter("Shapefile (*.shp)")

    def loadData(self, w, h):
        # Load data
        polygons = []
        if self.dia.exec():
            fileNames = self.dia.selectedFiles()
            with fiona.open(fileNames[0], 'r') as shapefile:
                for feature in shapefile:
                    geom = shape(feature['geometry'])
                    if isinstance(geom, Polygon):
                        qpolygon = QPolygonF()
                        # Transform coordinates
                        for coords in geom.exterior.coords:
                            x, y = coords[:2]  # Read only the first two values to handle potential 3D coordinates
                            qpolygon.append(QPointF(x, -y))  # Assuming y needs to be inverted
                        polygons.append(qpolygon)

        return self.transformPolygons(polygons, w, h)

    def transformPolygons(self, polygons, w, h):
        if not polygons:
            return []

        # Find the bounding box of all polygons
        min_x = min(min(poly[i].x() for i in range(poly.count())) for poly in polygons)
        max_x = max(max(poly[i].x() for i in range(poly.count())) for poly in polygons)
        min_y = min(min(poly[i].y() for i in range(poly.count())) for poly in polygons)
        max_y = max(max(poly[i].y() for i in range(poly.count())) for poly in polygons)

        # Determine scale and offset
        scale_x = w / (max_x - min_x)
        scale_y = h / (max_y - min_y)
        scale = min(scale_x, scale_y)

        offset_x = -min_x * scale + (w - (max_x - min_x) * scale) / 2
        offset_y = -min_y * scale + (h - (max_y - min_y) * scale) / 2

        # Transform polygons
        transformed_polygons = []
        for poly in polygons:
            polygon = QPolygonF()
            points = []
            for i in range(poly.count()):
                x = poly[i].x() * scale + offset_x
                y = poly[i].y() * scale + offset_y
                
                if x not in points:
                    points.append(x)
                    point2 = QPointF(x,y)
                    polygon.append(point2)
                    
            transformed_polygons.append(polygon)

        return transformed_polygons