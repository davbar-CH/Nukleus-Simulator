
import pyvistaqt as pvqt
from PyQt5.QtWidgets import *

import re

def text_auslesen(input_text):
    """
    Die Funktion liest mit Regex den Text aus und berechnet die Neutronenzahl

    :param input_text: Text aus der Textbox, in der Form "Element-Massenzahl".
    :return: gibt den Namen des Elementes, die Massenzahl, die Neutronenzahl und Ordnungszahl zurück.
    """
    stereo_pattern = r"\(([EZRSezrs,\d]+)\)-"
    substituenten_pattern = (r"(?:(\d+(?:,\d+)*)-)?(?:(di|tri|tetra|penta|hexa|hepta|octa|nona|deca))?("
                             r"fluor|chlor|brom|iod|methyl|ethyl|propyl|butyl|pentyl|hexyl|heptyl|phenyl|nitro"
                             r"|hydroxy|amino|oxo)")

    stereo = re.findall(stereo_pattern, input_text,flags=re.IGNORECASE)

    substituent = re.findall(substituenten_pattern, input_text,flags=re.IGNORECASE)

    it_ohne_ster = re.sub(stereo_pattern, '', input_text, flags=re.IGNORECASE)

    it_ohne_ster_sub = re.sub(substituenten_pattern, '', it_ohne_ster, flags=re.IGNORECASE)

    isCyclo = False
    if re.search(r"cyclo", it_ohne_ster_sub,re.IGNORECASE):
        isCyclo = True

    stamm_pattern = r"(?:cyclo)?(meth|eth|prop|but|pent|hex|hept|oct|non|dec)"
    stamm = re.findall(stamm_pattern, it_ohne_ster_sub,flags=re.IGNORECASE)

    bindung_pattern = r"(?:(\d+(?:,\d+)*)-)?(?:(di|tri|tetra|penta|hexa|hepta|octa|nona|deca))?(en|in)(?!\w)"
    bindung = re.findall(bindung_pattern, it_ohne_ster_sub,flags=re.IGNORECASE)

    endung_pattern_saeure_al = r"(säure|al)"
    endung_saeure_al = re.findall(endung_pattern_saeure_al, it_ohne_ster_sub,flags=re.IGNORECASE)

    endung_pattern_rest = r"(?:(\d+(?:,\d+)*)-)?(?:(di|tri|tetra|penta|hexa|hepta|octa|nona|deca))?(ol|amin)"
    endung_rest = re.findall(endung_pattern_rest, it_ohne_ster_sub, flags=re.IGNORECASE)

    return stereo, substituent, isCyclo, stamm, bindung, endung_saeure_al, endung_rest

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("Nukleus Simulation")
        self.setGeometry(100, 100, 800, 600)

        layout_vertikal = QVBoxLayout()
        layout_horizontal = QHBoxLayout()

        self.infos = QLabel("")
        layout_horizontal.addWidget(self.infos)

        self.textbox = QTextEdit(self)
        layout_horizontal.addWidget(self.textbox)

        layout_horizontal.setSpacing(150)

        layout_vertikal.addLayout(layout_horizontal)

        self.plotter = pvqt.QtInteractor(self)
        layout_vertikal.addWidget(self.plotter.interactor)

        self.button = QPushButton("Start")
        layout_vertikal.addWidget(self.button)
        self.button.clicked.connect(self.update_plot)

        central_widget = QWidget()
        central_widget.setLayout(layout_vertikal)
        self.setCentralWidget(central_widget)

    def update_plot(self):
        input_text = self.textbox.toPlainText()
        print(input_text)
        resultat = text_auslesen(input_text)

        if not resultat:
            print("Keine Eingabe")
            return

        stereo, substituent, isCyclo, stamm, bindung, endung_saeure_al, endung_rest = resultat
        print(stereo, substituent, isCyclo, stamm, bindung, endung_saeure_al, endung_rest)

app = QApplication([])
window = MainWindow()
window.show()
app.exec_()