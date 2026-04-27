import pyvista as pv
import pyvistaqt as pvqt
from PyQt5.QtWidgets import *
import numpy as np
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

def darsteller(stereo, substituent, isCyclo, stamm, bindung, endung_saeure_al, endung_rest, plotter):
    plotter.clear()
    stamm_laenge = {
        "Eth": 2,
        "eth": 2,
        "Prop": 3,
        "prop": 3,
        "But": 4,
        "but": 4,
        "Pent": 5,
        "pent": 5,
        "Hex": 6,
        "hex": 6,
        "Hept": 7,
        "hept": 7,
        "Oct": 8,
        "oct": 8,
        "Non": 9,
        "non": 9,
        "Dec": 10,
        "dec": 10,
    }
    substituent_laenge = {
        "Methyl": 1,
        "methyl": 1,
        "Ethyl": 2,
        "ethyl": 2,
        "Propyl": 3,
        "propyl": 3,
        "Butyl": 4,
        "butyl": 4,
        "Pentyl": 5,
        "pentyl": 5,
        "Hexyl": 6,
        "hexyl": 6,
        "Heptyl": 7,
        "heptyl": 7,
        "Octyl": 8,
        "octyl": 8,
        "Nonyl": 9,
        "nonyl": 9,
        "Decyl": 10,
        "decyl": 10,
    }

    stamm_kette_punkte = np.array([[x * 0.5, (1 - (-1) ** x) / 2, 0]for x in range(0, stamm_laenge.get(stamm[0]))])
    print(stamm_kette_punkte)
    stamm_kette = pv.lines_from_points(stamm_kette_punkte)
    plotter.add_mesh(stamm_kette,line_width=4)

    try:
        alle_sub_pos = re.findall(r"\d", substituent[0][0])
        alle_sub_pos = [int(x) for x in alle_sub_pos]

        besetzt_liste = []

        for sub_pos in alle_sub_pos:
            if sub_pos in besetzt_liste:
                if sub_pos % 2 == 0:

                    sub_alkan_punkte = np.array([[stamm_kette_punkte[sub_pos-1][0] - ((-1 + (-1) ** x) / 8), (x + 2) * 0.5, 0] for x in
                                                 range(0, substituent_laenge.get(substituent[0][2])+1)])
                    besetzt_liste.append(sub_pos)
                    sub_alkan_kette = pv.lines_from_points(sub_alkan_punkte)
                    plotter.add_mesh(sub_alkan_kette, line_width=2)
                else:
                    sub_alkan_punkte = np.array(
                        [[stamm_kette_punkte[sub_pos - 1][0] - ((-1 + (-1) ** x) / 8), x * -0.5, 0] for x in
                         range(0, substituent_laenge.get(substituent[0][2]) + 1)])
                    besetzt_liste.append(sub_pos)
                    sub_alkan_kette = pv.lines_from_points(sub_alkan_punkte)
                    plotter.add_mesh(sub_alkan_kette, line_width=2)
                print(sub_alkan_punkte)
            else:
                if sub_pos % 2 == 0:

                    sub_alkan_punkte = np.array(
                        [[stamm_kette_punkte[sub_pos - 1][0] + ((-1 + (-1) ** x) / 8), (x + 2) * 0.5, 0] for x in
                         range(0, substituent_laenge.get(substituent[0][2]) + 1)])
                    besetzt_liste.append(sub_pos)
                    sub_alkan_kette = pv.lines_from_points(sub_alkan_punkte)
                    plotter.add_mesh(sub_alkan_kette, line_width=2)
                else:
                    sub_alkan_punkte = np.array(
                        [[stamm_kette_punkte[sub_pos - 1][0] + ((-1 + (-1) ** x) / 8), x * -0.5, 0] for x in
                         range(0, substituent_laenge.get(substituent[0][2]) + 1)])
                    besetzt_liste.append(sub_pos)
                    sub_alkan_kette = pv.lines_from_points(sub_alkan_punkte)
                    plotter.add_mesh(sub_alkan_kette, line_width=2)
                print(sub_alkan_punkte)


    except:
        print("kein Substituent")


    plotter.add_axes()
    plotter.camera_position = 'xy'
    plotter.render()

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
        darsteller(
            stereo,
            substituent,
            isCyclo,
            stamm,
            bindung,
            endung_saeure_al,
            endung_rest,
            self.plotter
        )


app = QApplication([])
window = MainWindow()
window.show()
app.exec_()